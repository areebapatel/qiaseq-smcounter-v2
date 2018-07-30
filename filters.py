import os
import sys
import scipy.stats
import fisher

#-------------------------------------------------------------------------------------
# filter "LM": low coverage
#-------------------------------------------------------------------------------------
def lm(fltrs, sMtCons):
   if sMtCons < 5:
      fltrs.add('LM') 
   # output variables
   return(fltrs)

#-------------------------------------------------------------------------------------
# Initial check of HP and LowC; NOTE: the results are not final filters and may be reversed by hp4indel() and rep4others()
#-------------------------------------------------------------------------------------
def isHPorLowComp(chrom, pos, length, refb, altb, refs, repTypeSet, hpInfo):
   # ref sequence of [pos-length, pos+length] interval
   chromLength = refs.get_reference_length(chrom)
   pos0 = int(pos) - 1   
   Lseq = refs.fetch(reference=chrom,start=max(0,pos0-length),end=pos0).upper()
   Rseq_ref = refs.fetch(reference=chrom,start=pos0+len(refb),end=min(pos0+len(refb)+length,chromLength)).upper()  
   Rseq_alt = refs.fetch(reference=chrom,start=min(pos0+len(altb),chromLength), end=min(pos0+len(altb)+length,chromLength)).upper()
   refSeq = Lseq + refb + Rseq_ref
   altSeq = Lseq + altb + Rseq_alt
   # check homopolymer
   homoA = refSeq.find('A'*length) >= 0 or altSeq.find('A'*length) >= 0
   homoT = refSeq.find('T'*length) >= 0 or altSeq.find('T'*length) >= 0
   homoG = refSeq.find('G'*length) >= 0 or altSeq.find('G'*length) >= 0
   homoC = refSeq.find('C'*length) >= 0 or altSeq.find('C'*length) >= 0
   homop = homoA or homoT or homoG or homoC

   # check low complexity -- window length is 2 * homopolymer region. If any 2 nucleotide >= 99% 
   len2 = 2 * length
   LseqLC = refs.fetch(reference=chrom,start=max(0,pos0-len2),end=pos0).upper()
   Rseq_refLC = refs.fetch(reference=chrom,start=pos0+len(refb),end=min(pos0+len(refb)+len2,chromLength)).upper()
   Rseq_altLC = refs.fetch(reference=chrom,start=min(pos0+len(altb),chromLength),end=min(pos0+len(altb)+len2,chromLength)).upper()
   # ref seq   
   refSeqLC = LseqLC + refb + Rseq_refLC
   # alt seq
   altSeqLC = LseqLC + altb + Rseq_altLC
   lowcomp = False

   # Ref seq
   totalLen = len(refSeqLC)
   for i in range(totalLen-len2):
      subseq = refSeqLC[i:(i+len2)]
      countA = subseq.count('A')
      countT = subseq.count('T')
      countG = subseq.count('G')
      countC = subseq.count('C')
      sortedCounts = sorted([countA, countT, countG, countC], reverse=True)
      top2Freq = 1.0 * (sortedCounts[0] + sortedCounts[1]) / len2
      if top2Freq >= 0.99:
         lowcomp = True
         break
      
   # If ref seq is not LC, check alt seq
   if not lowcomp:
      totalLen = len(altSeqLC)
      for i in range(totalLen-len2):
         subseq = altSeqLC[i:(i+len2)]
         countA = subseq.count('A')
         countT = subseq.count('T')
         countG = subseq.count('G')
         countC = subseq.count('C')
         sortedCounts = sorted([countA, countT, countG, countC], reverse=True)
         top2Freq = 1.0 * (sortedCounts[0] + sortedCounts[1]) / len2
         if top2Freq >= 0.99:
            lowcomp = True
            break
   
   
   if homop and hpInfo == '.':   # if REF is not HP but ALT is, count as HP and set length = 8
      repTypeSet.add('HP')
      hpInfo = 'chr0;100;108;8;100;108'
   if lowcomp:  
      repTypeSet.add('LowC')
   
   # output variables
   return(repTypeSet, hpInfo)
  
#-------------------------------------------------------------------------------------
# model-based homopolymer filter for indel
#-------------------------------------------------------------------------------------
def hp4indel(fltrs, repTypeSet, vtype, rpb, hpInfo, vafToVmfRatio, vmf_tmp, hqUmiEff, RppEffSize, altRppUmiMean):
   if 'HP' in repTypeSet and vtype == 'INDEL':
      if rpb > 1.8:
         hpChrom, hpStart, hpEnd, totalHpLen, realL, realR = hpInfo.strip().split(';')
         hpLen8 = 1 if int(totalHpLen) >= 8 else 0
         intercept = -1.65834
         b_vafToVmfRatio = -1.52957
         b_vmf = 0.04744
         b_hqUmiEff = 3.01795
         b_RppEffSize = -0.06622
         b_altRppUmiMean = 0.51300
         b_hpLen8 = -0.86837

         cutoffHP = 0.565415

         predHP = intercept + b_vafToVmfRatio * vafToVmfRatio + b_vmf * vmf_tmp + b_hqUmiEff * hqUmiEff + b_RppEffSize * RppEffSize + b_altRppUmiMean * altRppUmiMean + b_hpLen8 * hpLen8  
         isReal = True if predHP >= cutoffHP else False
         
      else:
         isReal = False

      if isReal and 'HP' in fltrs:
         fltrs.remove('HP')
      if not isReal:
         fltrs.add('HP')

   # output variables
   return(fltrs)

#-------------------------------------------------------------------------------------
# heuristic thresholds for other types of repetitive filters
#-------------------------------------------------------------------------------------
def rep4others(fltrs, repTypeSet, vtype, rpb, vafToVmfRatio, hqUmiEff, RppEffSize):
   if len(repTypeSet) > 0 and (vtype == 'SNP'  or (vtype == 'INDEL' and 'HP' not in repTypeSet)):
      if rpb >= 4:
         isReal = hqUmiEff > 0.1 and vafToVmfRatio < 3.0 and RppEffSize < 2.5
      elif rpb >= 1.8:
         isReal = hqUmiEff > 0.8 and vafToVmfRatio < 2.0 and RppEffSize < 1.5
      else:
         isReal = False

      if isReal:
         fltrs.difference_update(repTypeSet)
      else:
         fltrs.update(repTypeSet)

   # output variables
   return(fltrs)

#-------------------------------------------------------------------------------------
# filter "DP" and "SB": discordant read pairs and strand bias 
#-------------------------------------------------------------------------------------
def dp_sb(fltrs, origAlt, concordPairCnt, discordPairCnt, reverseCnt, forwardCnt, origRef, vaf_tmp):
   pairs = discordPairCnt[origAlt] + concordPairCnt[origAlt] # total number of paired reads covering the pos
   pDiscord = 1.0 * discordPairCnt[origAlt] / pairs if pairs > 0 else 0.0
   if pairs >= 1000 and pDiscord >= 0.5:
      fltrs.add('DP') 
   elif vaf_tmp <= 60.0:
      refR = reverseCnt[origRef]
      refF = forwardCnt[origRef]
      altR = reverseCnt[origAlt]
      altF = forwardCnt[origAlt]
      
      if (refF == 0 or altR == 0) and (refR == 0 or altF == 0):
         oddsRatio = float("Nan")
         return(fltrs)
      elif refF == 0 or altR == 0:
         oddsRatio = float("inf")
      else:
         oddsRatio = float(refR*altF)/(refF*altR)
         
      if oddsRatio < 50 and oddsRatio > 1.0/50:
         return(fltrs)
      
      temp = fisher.pvalue(refR,refF,altR,altF)
      pvalue = temp.two_tail
      if pvalue < 0.00001:
         fltrs.add('SB')
   # output variables
   return(fltrs)

#-------------------------------------------------------------------------------------
# filter "PB": primer direction bias
#-------------------------------------------------------------------------------------
def pb(fltrs, origAlt, sMtConsByDir, sMtConsByDirByBase):
   if sMtConsByDir['F'] >= 200 and sMtConsByDir['R'] >= 200:
      oddsRatioPB = ((sMtConsByDirByBase[origAlt]['F']+0.5)/(sMtConsByDir['F']+0.5)) / ((sMtConsByDirByBase[origAlt]['R']+.5)/(sMtConsByDir['R']+.5))
      oddsRatioPB = round(oddsRatioPB, 2)
      if oddsRatioPB > 10 or oddsRatioPB < 0.1:
         fltrs.add('PB')
      primerBiasOR = str(oddsRatioPB)
   else:
      primerBiasOR = 'NA'

   # output variables
   return(fltrs, primerBiasOR)

#-------------------------------------------------------------------------------------
# filter "LowQ": reads supporting the variant have low base quality 
#-------------------------------------------------------------------------------------
def lowq(fltrs, lowQReads, alleleCnt, origAlt, vafToVmfRatio, bqAlt):
   intercept = 7.652256
   bRatio = -1.254942
   bPLowQ = -6.602585
   cutoffLowQ = 1.068349

   if origAlt in alleleCnt and origAlt in lowQReads and alleleCnt[origAlt] > 0:
      predLowQ = intercept +  bRatio * vafToVmfRatio + bPLowQ * bqAlt
      isLowQ = True if predLowQ < cutoffLowQ else False
	  
      if bqAlt > 0.4 and vafToVmfRatio >= 0 and isLowQ:
         fltrs.add('LowQ')

   # output variables
   return(fltrs)

#-------------------------------------------------------------------------------------
# filter "RBCP": variant too close to the random (UMI) end on R1
#-------------------------------------------------------------------------------------
def rbcp(fltrs, endBase, mtSideBcEndPos, origRef, origAlt, vaf_tmp):
   refLeEnd = sum(d <= endBase for d in mtSideBcEndPos[origRef])  # number of REF R2 reads with distance <= endBase
   refGtEnd = len(mtSideBcEndPos[origRef]) - refLeEnd         # number of REF R2 reads with distance > endBase
   altLeEnd = sum(d <= endBase for d in mtSideBcEndPos[origAlt])  # number of ALT R2 reads with distance <= endBase
   altGtEnd = len(mtSideBcEndPos[origAlt]) - altLeEnd         # number of ALT R2 reads with distance > endBase
   fisher = scipy.stats.fisher_exact([[refLeEnd, refGtEnd], [altLeEnd, altGtEnd]])
   oddsRatio = fisher[0]
   pvalue = fisher[1]
   if pvalue < 0.001 and oddsRatio < 0.05 and vaf_tmp <= 60.0:
      fltrs.add('RBCP')

   # output variables
   return(fltrs)

#-------------------------------------------------------------------------------------
# filter "RPCP": variant too close to the random (UMI) end on R2
#-------------------------------------------------------------------------------------
def rpcp(fltrs, endBase, primerSideBcEndPos, origRef, origAlt, vaf_tmp):
   refLeEnd = sum(d <= endBase for d in primerSideBcEndPos[origRef])  # number of REF R2 reads with distance <= endBase
   refGtEnd = len(primerSideBcEndPos[origRef]) - refLeEnd         # number of REF R2 reads with distance > endBase
   altLeEnd = sum(d <= endBase for d in primerSideBcEndPos[origAlt])  # number of ALT R2 reads with distance <= endBase
   altGtEnd = len(primerSideBcEndPos[origAlt]) - altLeEnd         # number of ALT R2 reads with distance > endBase
   fisher = scipy.stats.fisher_exact([[refLeEnd, refGtEnd], [altLeEnd, altGtEnd]])
   oddsRatio = fisher[0]
   pvalue = fisher[1]
   if pvalue < 0.001 and oddsRatio < 0.05 and vaf_tmp <= 60.0:
      fltrs.add('RPCP')

   # output variables
   return(fltrs)

#-------------------------------------------------------------------------------------
# filter "PrimerCP": variant too close to the fixed (gene specific primer) end 
#-------------------------------------------------------------------------------------
def primercp(fltrs, primerDist, primerSidePrimerEndPos, origRef, origAlt, vmf_tmp, hqUmiEff, vafToVmfRatio, RppEffSize, rpb):
   # fixed end position filter
   refLeEnd = sum(d <= primerDist for d in primerSidePrimerEndPos[origRef])  # number of REF R2 reads with distance <= endBase
   refGtEnd = len(primerSidePrimerEndPos[origRef]) - refLeEnd         # number of REF R2 reads with distance > endBase
   altLeEnd = sum(d <= primerDist for d in primerSidePrimerEndPos[origAlt])  # number of ALT R2 reads with distance <= endBase
   altGtEnd = len(primerSidePrimerEndPos[origAlt]) - altLeEnd         # number of ALT R2 reads with distance > endBase
   fisher = scipy.stats.fisher_exact([[refLeEnd, refGtEnd], [altLeEnd, altGtEnd]])
   oddsRatio = fisher[0]
   pvalue = fisher[1]
   # updated PrimerCP -- depend on UMI efficiency
   if vmf_tmp < 40.0 and (altLeEnd + altGtEnd > 0) and (1.0 * altLeEnd / (altLeEnd + altGtEnd) >= 0.98 or (pvalue < 0.001 and oddsRatio < 0.05)):
      if rpb >= 4:
         isReal = hqUmiEff > 0.1 and vafToVmfRatio < 3.0 and RppEffSize < 2.5
      elif rpb >= 1.8:
         isReal = hqUmiEff > 0.8 and vafToVmfRatio < 2.0 and RppEffSize < 1.5
      else:
         isReal = False
      if not isReal:
         fltrs.add('PrimerCP')

   # output variables
   return(fltrs)

#-------------------------------------------------------------------------------------
# strict filters used in smCounter-v1 and earlier versions of smCounter2; for consensused BAM only
#-------------------------------------------------------------------------------------
def strict(fltrs, repTypeSet, vtype, vmf_tmp, primerDist, primerSidePrimerEndPos, origRef, origAlt):
   # homopolymer
   if 'HP' in repTypeSet and vmf_tmp < 99.0:
      fltrs.add('HP')
   # low complexity
   if 'LowC' in repTypeSet and vmf_tmp < 99.0:
      fltrs.add('LowC')
   # tandom repeat
   if 'RepT' in repTypeSet and vmf_tmp < 40.0:
      fltrs.add('RepT')
   # simple repeat
   if 'RepS' in repTypeSet and vmf_tmp < 40.0:
      fltrs.add('RepS')  
   # satellites
   if 'SL' in repTypeSet and vmf_tmp < 40.0:
      fltrs.add('SL')      	

   # PrimerCP (SNP only)
   if vtype == 'SNP':
      refLeEnd = sum(d <= primerDist for d in primerSidePrimerEndPos[origRef])  # number of REF R2 reads with distance <= endBase
      refGtEnd = len(primerSidePrimerEndPos[origRef]) - refLeEnd         # number of REF R2 reads with distance > endBase
      altLeEnd = sum(d <= primerDist for d in primerSidePrimerEndPos[origAlt])  # number of ALT R2 reads with distance <= endBase
      altGtEnd = len(primerSidePrimerEndPos[origAlt]) - altLeEnd         # number of ALT R2 reads with distance > endBase
      fisher = scipy.stats.fisher_exact([[refLeEnd, refGtEnd], [altLeEnd, altGtEnd]])
      oddsRatio = fisher[0]
      pvalue = fisher[1]   
      if altLeEnd + altGtEnd > 0 and 1.0 * altLeEnd / (altLeEnd + altGtEnd) >= 0.98 or (pvalue < 0.001 and oddsRatio < 0.05):
         fltrs.add('PrimerCP')
   
   # output variables
   return(fltrs)
