#!/usr/bin/python

# Chang Xu, 23Oct2017
import os
import sys
import datetime
import subprocess
import time
import operator
import multiprocessing
from collections import defaultdict
import random
import traceback

# 3rd party modules
import argparse
import pysam
import scipy.stats
import numpy

codePath = os.path.dirname(os.path.abspath(__file__))

homopolymerCode = os.path.join(codePath,'findhp.py')
pValCode = os.path.join(codePath,'get_pvalue_v2.R')
vcfCode = os.path.join(codePath,'make_vcf_v2.py')

atgc = ['A', 'T', 'G', 'C']
seed = 10262016
nsim = 5000
minTotalUMI = 5

mtTag = "Mi"
tagSeparator = "-"
primerTag = "pr"

_num_cols_ = 38 ## Number of columns in out_long returned by the vc() function of smCounter



# wrapper function for "vc()" - because Python multiprocessing module does not pass stack trace
# from runone/smcounter.py by John Dicarlo
#------------------------------------------------------------------------------------------------
def vc_wrapper(*args):
   try:
      output = vc(*args)
   except Exception as e:
      print("Exception thrown in vc() function at genome location:", args[1], args[2])
      output = ("Exception thrown!\n" + traceback.format_exc(),'no_bg')
      print output[0]
      raise Exception(e)
   return output
	
   
#-------------------------------------------------------------------------------------
# calculate mean rpb
#-------------------------------------------------------------------------------------
def getMeanRpb(bamName):
   samfile = pysam.AlignmentFile(bamName, 'rb')
   allFragSet = set()
   allBcSet = set()

   # fetch all reads
   for read in samfile.fetch():
      # read ID
      qname = read.query_name
      
      # barcode sequence
      BC = read.get_tag(mtTag)

      allFragSet.add(qname)
      allBcSet.add(BC)

   # total fragment count
   totalFrag = len(allFragSet)
   # total MT count
   totalMT = len(allBcSet)
   # mean rpb
   meanRpb = float(totalFrag) / totalMT
   samfile.close()
   return meanRpb

#-------------------------------------------------------------------------------------
# get reverse complement of bases
#-------------------------------------------------------------------------------------
def reverseBase(base):
   if base == 'A':
      revBase = 'T'
   elif base == 'T':
      revBase = 'A'
   elif base == 'G':
      revBase = 'C'
   elif base == 'C':
      revBase = 'G'
   else:
      revBase = 'N'
   return revBase
   
#-------------------------------------------------------------------------------------
# model-based homopolymer filter for indels
#-------------------------------------------------------------------------------------
def hpIndelPredict(vafToVmfRatio, vmf, hqUmiEff, rpbEffectSize, altMeanRpb, hpLen8):
   intercept = -1.65834
   b_vafToVmfRatio = -1.52957
   b_vmf = 0.04744
   b_hqUmiEff = 3.01795
   b_rpbEffectSize = -0.06622
   b_altMeanRpb = 0.51300
   b_hpLen8 = -0.86837

   cutoffHP = 0.565415
   
   predHP = intercept + b_vafToVmfRatio * vafToVmfRatio + b_vmf * vmf + b_hqUmiEff * hqUmiEff + b_rpbEffectSize * rpbEffectSize + b_altMeanRpb * altMeanRpb + b_hpLen8 * hpLen8  
   isRealIndelHP = True if predHP >= cutoffHP else False
   return isRealIndelHP

#-------------------------------------------------------------------------------------
# model-based LowQ filter for SNP
#-------------------------------------------------------------------------------------
def lowQPredict(bqAlt, vafToVmfRatio):
   intercept = 7.652256
   bRatio = -1.254942
   bPLowQ = -6.602585
   cutoffLowQ = 1.068349

   predLowQ = intercept +  bRatio * vafToVmfRatio + bPLowQ * bqAlt
   isLowQ = True if predLowQ < cutoffLowQ else False
   return isLowQ

#-------------------------------------------------------------------------------------
# find the consensus nucleotide (including indel) in a UMI family with high quality reads only
#-------------------------------------------------------------------------------------
def consHqMT(oneBC, mtThreshold):
   totalCnt = oneBC['all']
   cons = ''
   # find the majority base(s) whose proportion in the MT >= mtThreshold. NOTE: mtThreshold must be > 0.5 to ensure only one cons 
   for base in oneBC:
      if base == 'all':
         continue
      pCons = 1.0 * oneBC[base] / totalCnt if totalCnt > 0 else 0.0
      if pCons >= mtThreshold:
         cons = base
         break
   # report the consensus base. If no consensus or lack of read support, output ''. 
   return cons

#-------------------------------------------------------------------------------------
# find the consensus nucleotide (including indel) in a UMI family with all reads 
#-------------------------------------------------------------------------------------
def consAllMT(readList, mtThreshold):
   totalCnt = readList['all']
   cons = ''
   # find the majority base(s) whose proportion in the MT >= mtThreshold. NOTE: mtThreshold must be > 0.5 to ensure only one cons
   for base in readList:
      if base == 'all': ## just a counter
         continue
      pCons = 1.0 * readList[base] / totalCnt if totalCnt > 0 else 0.0
      if pCons >= mtThreshold:
         cons = base
         break
   # report the consensus base. If no consensus or lack of read support, output ''. 
   return cons

#-------------------------------------------------------------------------------------
# check if a locus is within or flanked by homopolymer region and/or low complexity region
#-------------------------------------------------------------------------------------
def isHPorLowComp(chrom, pos, length, refb, altb, refs):
   # get reference base
   #refs = pysam.FastaFile(refg)
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

   return [homop, lowcomp]

   
#-------------------------------------------------------------------------------------
# function to call variants
#-------------------------------------------------------------------------------------
def vc(bamName, chrom, pos, repType, hpInfo, srInfo, repInfo, minBQ, minMQ, hpLen, mismatchThr, primerDist, mtThreshold, rpb, primerSide, refg, minAltUMI, maxAltAllele, isRna, ds):   
   samfile = pysam.AlignmentFile(bamName, 'rb')

   mtSide = 'R1' if primerSide == 'R2' else 'R2'
   cvg = 0
   
   bcDictHq = defaultdict(lambda: defaultdict(list))
   bcDictHqBase = defaultdict(lambda:defaultdict(int))
   usedFrag = 0
   bcDictAll = defaultdict(lambda:defaultdict(int))
   allBcDict = defaultdict(set)
   allFrag = 0

   alleleCnt = defaultdict(int)
   mtSideBcEndPos = defaultdict(list)

   primerSideBcEndPos = defaultdict(list)
   primerSidePrimerEndPos = defaultdict(list)

   forwardCnt = defaultdict(int)
   reverseCnt = defaultdict(int)
   concordPairCnt = defaultdict(int)
   discordPairCnt = defaultdict(int)
   mismatchCnt = defaultdict(float)

   lowQReads = defaultdict(int)
   sMtCons = 0
   sMtConsByBase = defaultdict(int)
   sMtConsByDir  = defaultdict(int)
   sMtConsByDirByBase = defaultdict(lambda: defaultdict(int))
   pairedMTs = set()
   singleMTs = set()

   strands = defaultdict(int)
   subTypeCnt = defaultdict(int)
   smtSNP = 0

   repTypeSet0 = set() if repType == 'NA' else set(repType.strip().split(';'))

   hqAgree = defaultdict(int)
   hqDisagree = defaultdict(int)
   allAgree = defaultdict(int)
   allDisagree = defaultdict(int)
   rpbCnt = defaultdict(list)

   # output 
   out_long = ''
   # get reference base
   refseq = pysam.FastaFile(refg)
   origRef = refseq.fetch(reference=chrom, start=int(pos)-1, end=int(pos))
   origRef = origRef.upper()
   # splitting hpInfo here to avoid splitting inside the loop below
   hpInfoTmp = hpInfo.strip().split(';')

   sMtConsByBase['A'] = 0
   sMtConsByBase['T'] = 0
   sMtConsByBase['G'] = 0
   sMtConsByBase['C'] = 0
   
   downsamplePileupStackThr = 10 ** 5

   # pile up reads
   for read in samfile.pileup(region = chrom + ':' + pos + ':' + pos, truncate=True, max_depth=1000000000, stepper='nofilter'):
      # check pielup size, downsample if above threshold
      pileupStackSize = read.n
      downsamplePileup = True if pileupStackSize > downsamplePileupStackThr else False
      random.seed(pos)
   
      for pileupRead in read.pileups:
         # drop randomly
         if downsamplePileup and random.randint(1, pileupStackSize) > downsamplePileupStackThr:
            continue
            
         # check if position not on a gap (N or intron in RNAseq)
         dropRead = False
         cigar = pileupRead.alignment.cigar
         alignLen = int(pos) - pileupRead.alignment.pos
         # first case: perhaps outside the whole read
         if alignLen > sum([value if op in [0, 3] else 0 for (op, value) in cigar]):
            continue
         # second case: may lie on any segments
         for (op, value) in cigar:
            if op > 3:
               continue
            if alignLen <= value:
               if op == 3:
                  dropRead = True
               break
            elif op in [0, 3]:
               alignLen -= value
               
         # check if should drop read
         if dropRead:
            continue                       
         
         # read ID
         readid = pileupRead.alignment.query_name         
         BC = pileupRead.alignment.get_tag(mtTag)

         # read start and end coordinates in reference genome
         astart = pileupRead.alignment.reference_start
         aend = pileupRead.alignment.reference_end

         # get NM tag 
         NM = 0 
         allTags = pileupRead.alignment.tags
         for (tag, value) in allTags:
            if tag == 'NM':
               NM = value
               break
         # count number of INDELs in the read sequence
         nIndel = 0
         cigar = pileupRead.alignment.cigar
         cigarOrder = 1
         leftSP = 0  # soft clipped bases on the left
         rightSP = 0  # soft clipped bases on the right
         for (op, value) in cigar:
            # 1 for insertion
            if op == 1 or op == 2:
               nIndel += value 
            if cigarOrder == 1 and op == 4:
               leftSP = value
            if cigarOrder > 1 and op == 4:
               rightSP += value
            cigarOrder += 1

         # Number of mismatches except INDEL, including softcilpped sequences 
         mismatch = max(0, NM - nIndel)
         # read length, including softclip
         readLen = pileupRead.alignment.query_length
         # calculate mismatch per 100 bases
         mismatchPer100b = 100.0 * mismatch / readLen if readLen > 0 else 0.0

         # paired read
         if pileupRead.alignment.is_read1:
            pairOrder = 'R1'
         if pileupRead.alignment.is_read2:
            pairOrder = 'R2'

         # +/- strand
         strand = 'Reverse' if pileupRead.alignment.is_reverse else 'Forward'

         # mapping quality filter
         mq = pileupRead.alignment.mapping_quality
         minMQPass = True
         # get mapq of mate
         try:
            mateMq = pileupRead.alignment.get_tag("MQ")
            minFragMQ = min(mq,mateMq)
            if minFragMQ < minMQ:
               minMQPass = False
         except KeyError: 
            '''
            bam has not been tagged with the mate mapq,
            drop read pairs based on their respective mapqs only
            To note :
            warn user ? or make command line argument more descriptive
            settling on a more descriptive argument for now
            '''
            if mq < minMQ:
               minMQPass = False

         # repetitive region information
         if hpInfo == '.':
            hpCovered = True
         else:
            (hpChrom, hpStart, hpEnd, totalHpLen) = hpInfoTmp
            if astart < int(hpStart) - 1 and aend > int(hpEnd) + 1:
               hpCovered = True
            else:
               hpCovered = False

         # check if the site is the beginning of insertion
         if pileupRead.indel > 0:
            site = pileupRead.alignment.query_sequence[pileupRead.query_position]
            inserted = pileupRead.alignment.query_sequence[(pileupRead.query_position + 1) : (pileupRead.query_position + 1 +  pileupRead.indel)]
            base = 'INS|' + site + '|' + site + inserted
            bq = pileupRead.alignment.query_qualities[pileupRead.query_position]
            # if base quality not included in BAM
            if bq == None:
               bq = minBQ         
            
         # check if the site is the beginning of deletion
         elif pileupRead.indel < 0:
            site = pileupRead.alignment.query_sequence[pileupRead.query_position]
            deleted = refseq.fetch(reference=chrom, start=int(pos), end=int(pos)+abs(pileupRead.indel))
            deleted = deleted.upper()
            base = 'DEL|' + site + deleted + '|' + site
            bq = pileupRead.alignment.query_qualities[pileupRead.query_position]
            # if base quality not included in BAM
            if bq == None:
               bq = minBQ         

         # site is not beginning of any INDEL, but in the middle of a deletion
         elif  pileupRead.is_del:
            base = 'DEL'
            bq = minBQ

         # if the site is a regular locus, 
         else: 
            base = pileupRead.alignment.query_sequence[pileupRead.query_position] # note: query_sequence includes soft clipped bases
            bq = pileupRead.alignment.query_qualities[pileupRead.query_position]
            incCond = bq >= minBQ and minMQPass and mismatchPer100b <= mismatchThr and hpCovered
            # count the number of low quality reads (less than Q20 by default) for each base
            if bq < 20:   # why not minBQ???!!!
               lowQReads[base] += 1
            if pairOrder == mtSide:
               # distance to the barcode end on MT side read 
               if pileupRead.alignment.is_reverse:
                  distToBcEnd = pileupRead.alignment.query_alignment_length - (pileupRead.query_position - leftSP)      
               else:
                  distToBcEnd = pileupRead.query_position - leftSP
               if incCond:
                  mtSideBcEndPos[base].append(distToBcEnd)
            if pairOrder == primerSide:
               # distance to the barcode and/or primer end on primer side read. Different cases for forward and reverse strand
               if pileupRead.alignment.is_reverse:
                  distToBcEnd = pileupRead.query_position - leftSP
                  distToPrimerEnd = pileupRead.alignment.query_alignment_length - (pileupRead.query_position - leftSP)
               else:
                  distToBcEnd = pileupRead.alignment.query_alignment_length - (pileupRead.query_position - leftSP)
                  distToPrimerEnd = pileupRead.query_position - leftSP
               if incCond:
                  primerSideBcEndPos[base].append(distToBcEnd)
                  primerSidePrimerEndPos[base].append(distToPrimerEnd)

         # coverage -- read, not fragment
         cvg += 1

         if strand == 'Reverse':
            reverseCnt[base] += 1
         else:
            forwardCnt[base] += 1
         alleleCnt[base] += 1
         mismatchCnt[base] += mismatchPer100b

         # count total number of fragments and MTs
         if readid not in allBcDict[BC]:
            allFrag+=1 # total fragments
            allBcDict[BC].add(readid)

         # inclusion condition. NOTE: reads with duplex tag 'NN' are dropped from analysis
         incCond = bq >= minBQ and minMQPass and mismatchPer100b <= mismatchThr and hpCovered

         # constructing UMI family; this one with high quality reads only
         if incCond:
            if readid not in bcDictHq[BC]:
               readinfo = [base, pairOrder]
               bcDictHq[BC][readid] = readinfo
               # store base level information to avoid looping over read ids again
               bcDictHqBase[BC][base]+=1
               bcDictHqBase[BC]['all']+=1
               usedFrag+=1 # used fragments
               
            elif base == bcDictHq[BC][readid][0] or base in ['N', '*']:
               bcDictHq[BC][readid][1] = 'Paired'
               if base == bcDictHq[BC][readid][0]:
                  concordPairCnt[base] += 1
            else:
               # decrement fragment and base count
               usedFrag-=1
               bcDictHqBase[BC][bcDictHq[BC][readid][0]]-=1
               bcDictHqBase[BC]['all']-=1
               del bcDictHq[BC][readid]
               discordPairCnt[base] += 1

         # in non-HP region, include all reads for consensus. In HP region, including only the reads covering the HP. 
         if hpCovered:
            #bcDictAll[BC].append(base)
            bcDictAll[BC][base]+=1
            bcDictAll[BC]['all']+=1

   ##### end of looping through pileup reads ####

   # total number of MT, fragments, reads, including those dropped from analysis
   allMT = len(allBcDict)
   # gradually drop 1 read MTs
   bcToKeep = []
   
   if isRna:
      rbp = cvg / float(allFrag) if allFrag > 0 else 1.0
   
   # rpb < 2: no MT is dropped
   if rpb < 2.0: 
      bcToKeep = bcDictHq.keys()

   # 2 <= rpb < 3: gradually and randomly drop 1 read MTs 
   elif rpb >= 2.0 and rpb < 3.0:
      # set seed to be the genome position
      random.seed(pos)
      # count the numbers of paired and unpaired 1 read MTs; 
      pctToDrop = rpb - 2.0
      for bc in bcDictHq:
         readPairsInBc = len(bcDictHq[bc])
         if readPairsInBc == 1:
            readid = bcDictHq[bc].keys()[0]
            if bcDictHq[bc][readid][1] == 'Paired':
               pairedMTs.add(bc)
            else:
               singleMTs.add(bc)
      # total 1 read MTs
      pairedCnt = len(pairedMTs)
      singleCnt = len(singleMTs)
      oneReadMtCnt = pairedCnt + singleCnt
      # number of 1 read MTs to drop
      numToDrop = int(round(pctToDrop * oneReadMtCnt))
      # Decide which 1 read MTs to drop -- paired reads are kept with priority
      if numToDrop <= singleCnt:
         oneReadMtToDrop = set(random.sample(singleMTs, numToDrop))
      else:
         pairsToDrop = set(random.sample(pairedMTs, numToDrop - singleCnt))
         oneReadMtToDrop = singleMTs.union(pairsToDrop)
      # drop 1 read MTs
      bcToKeep = list(set(bcDictHq.keys()).difference(oneReadMtToDrop))      

   #rpb >= 3: drop all 1 read MTs;
   else:
      bcToKeep = [bc for bc in bcDictHq.iterkeys() if len(bcDictHq[bc]) >= 2]

   bcToKeepLen = len(bcToKeep)
   if bcToKeepLen <= minTotalUMI:
      out_long = '\t'.join([chrom, pos, origRef] + ['0'] * (_num_cols_ - 4) + ['LM']) + '\n'
      out_bkg = ''
   else:
      # down sample if too many UMIs
      if isRna and bcToKeepLen > ds:
         random.seed(pos)
         bcToKeep = random.sample(bcToKeep, ds)
         
      for bc in bcToKeep:
         # primer ID and direction
         bcSplit = bc.split(tagSeparator)
         primerDirCode = bcSplit[1]
         primerDirection = 'F' if primerDirCode == '0' else 'R' # 0 means the primer was priming the forward strand, 1 means priming the reverse strand

         # get consensus call of the UMI family
         consHq = consHqMT(bcDictHqBase[bc], mtThreshold)
         consAll = consAllMT(bcDictAll[bc], mtThreshold)
         cons = consHq if consHq == consAll else ''
         
         # count number of reads in concordant/discordant with consensus
         for base in bcDictHqBase[bc]:
            if base == 'all': ## just a counter
               continue
            if base == cons:
               hqAgree[base] += bcDictHqBase[bc][base]
            else:
               hqDisagree[base] += bcDictHqBase[bc][base]
   
         for base in bcDictAll[bc]:
            if base == 'all': ## just a counter
               continue
            if base == cons:
               allAgree[base] += bcDictAll[bc][base]
            else:
               allDisagree[base] += bcDictAll[bc][base]

         if cons != '':
            sMtCons += 1
            sMtConsByBase[cons] += 1
            # MT counts from + and - strands 
            sMtConsByDir[primerDirection] += 1
            sMtConsByDirByBase[cons][primerDirection] += 1
            # read pairs in the UMI            
            rpbCnt[cons].append(bcDictAll[bc]['all'])
            
            # base substitutions (snp only)
            # Note: smtSNP and strands are usually NOT equal to sMtCons and sMtConsByDir. The former include only base substitutions MTs, and the latter include indel MTs. 
            if len(cons) == 1:
               basePair = origRef + '/' + cons if primerDirCode == '0' else reverseBase(origRef) + '/' + reverseBase(cons)
               subTypeCnt[basePair] += 1
               smtSNP += 1
               strands[primerDirection] += 1

      # output the background error profile
      bkgErrList = [chrom, pos, origRef, str(subTypeCnt['A/G']), str(subTypeCnt['G/A']), str(subTypeCnt['C/T']), str(subTypeCnt['T/C']), str(subTypeCnt['A/C']), str(subTypeCnt['C/A']), str(subTypeCnt['A/T']), str(subTypeCnt['T/A']), str(subTypeCnt['C/G']), str(subTypeCnt['G/C']), str(subTypeCnt['G/T']), str(subTypeCnt['T/G']), str(strands['F']), str(strands['R']), str(smtSNP)]
      out_bkg = '\t'.join(bkgErrList) + '\n'

      sortedList = sorted(sMtConsByBase.items(), key=operator.itemgetter(1), reverse=True)
      firstAlt = True
      altCnt = 0
      # start multi-allelic loop
      for alleleInd in range(len(sortedList)):
         maxBase = sortedList[alleleInd][0]
         maxVMT = sortedList[alleleInd][1]

         if maxBase == origRef:
            continue
         if maxVMT < minAltUMI and not firstAlt:
            break

         # if the current allele has >= 3 UMIs and is not reference, treat as a candidate variant allele
         origAlt = maxBase

         # reset variant type, reference base, variant base 
         vtype = '.'
         ref = origRef
         alt = origAlt

         if len(origAlt) == 1:
            vtype = 'SNP'
         elif origAlt == 'DEL':
            vtype = 'SDEL'
         else:
            vals = origAlt.split('|')
            if vals[0] in ['DEL', 'INS']:
               vtype = 'INDEL'
               ref = vals[1]
               alt = vals[2]

         # initiate values for filters
         refForPrimer = sMtConsByDirByBase[origRef]['F']
         refRevPrimer = sMtConsByDirByBase[origRef]['R']
         altForPrimer = sMtConsByDirByBase[origAlt]['F']
         altRevPrimer = sMtConsByDirByBase[origAlt]['R']
         primerBiasOR, bqAlt, oddsRatio, pvalue, hqUmiEff, allUmiEff, refRppUmiMean, altRppUmiMean, RppEffSize = 'NA', -1.0, 1.0, 1.0, 0.0, 0.0, -1.0, -1.0, -1.0
         fltrs = set()
         repTypeSet = repTypeSet0

         # UMI efficiency metrics
         hqRcAgree = hqAgree[origAlt] 
         hqRcTotal = hqRcAgree + hqDisagree[origAlt] 
         hqUmiEff = round(1.0 * hqRcAgree / hqRcTotal, 2) if hqRcTotal > 0 else 0.0
         
         allRcAgree = allAgree[origAlt] 
         allRcTotal = allRcAgree + allDisagree[origAlt] 
         allUmiEff = round(1.0 * allRcAgree / allRcTotal, 3) if allRcTotal > 0 else 0.0

         if sMtConsByBase[origRef] >= 3 and sMtConsByBase[origAlt] >= 3:
            refRppUmiN = sMtConsByBase[origRef]
            refRppUmiMean = numpy.mean(rpbCnt[origRef])
            refRppUmiSd = numpy.std(rpbCnt[origRef])
            altRppUmiN = sMtConsByBase[origAlt]
            altRppUmiMean = numpy.mean(rpbCnt[origAlt])
            altRppUmiSd = numpy.std(rpbCnt[origAlt])
            sp = ( ((refRppUmiN-1) * refRppUmiSd**2 + (altRppUmiN-1) * altRppUmiSd**2) / (refRppUmiN + altRppUmiN-2) ) ** 0.5
            RppEffSize = (refRppUmiMean - altRppUmiMean) / (sp * (1.0/refRppUmiN + 1.0/altRppUmiN) ** 0.5) if sp > 0 else 1000.0
         else:
            refRppUmiMean = -1.0
            altRppUmiMean = -1.0
            RppEffSize = -1.0

         if vtype in ['SNP', 'INDEL'] and sMtCons > 0:
            vaf_tmp = 100.0 * alleleCnt[origAlt] / cvg if cvg > 0  else 0.0
            vmf_tmp = 100.0 * sMtConsByBase[origAlt] / sMtCons
            vafToVmfRatio = 1.0 * vaf_tmp / vmf_tmp if vmf_tmp > 0 else -1.0
            # low coverage filter
            if sMtCons < 5:
               fltrs.add('LM') 

            # repetative region filters
            hpIndelFilter = False
            (hp, lowc) = isHPorLowComp(chrom, pos, hpLen, ref, alt, refseq)
            if hp and hpInfo == '.':   # if REF is not HP but ALT is, count as HP and set length = 8
               repTypeSet.add('HP')
               hpInfo = 'chr0;100;108;8'
            if lowc:  
               repTypeSet.add('LowC')

            # update HP for indel
            if 'HP' in repTypeSet and vtype == 'INDEL':
               if rpb > 1.8:
                  (hpChrom, hpStart, hpEnd, totalHpLen) = hpInfo.strip().split(';')
                  hpLen8 = 1 if int(totalHpLen) >= 8 else 0
                  isReal = hpIndelPredict(vafToVmfRatio, vmf_tmp, hqUmiEff, RppEffSize, altRppUmiMean, hpLen8)
               else:
                  isReal = False

               if isReal and 'HP' in fltrs:
                  fltrs.remove('HP')
               if not isReal:
                  fltrs.add('HP')
               
            # update other repetitive region filters for SNP and indel, including HP for SNP
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

            # strand bias and discordant pairs filter
            pairs = discordPairCnt[origAlt] + concordPairCnt[origAlt] # total number of paired reads covering the pos
            pDiscord = 1.0 * discordPairCnt[origAlt] / pairs if pairs > 0 else 0.0
            if pairs >= 1000 and pDiscord >= 0.5:
               fltrs.add('DP') 
            elif vaf_tmp <= 60.0:
               refR = reverseCnt[origRef]
               refF = forwardCnt[origRef]
               altR = reverseCnt[origAlt]
               altF = forwardCnt[origAlt]
               fisher = scipy.stats.fisher_exact([[refR, refF], [altR, altF]])
               oddsRatio = fisher[0]
               pvalue = fisher[1]
               if pvalue < 0.00001 and (oddsRatio >= 50 or oddsRatio <= 1.0/50):
                  fltrs.add('SB')

            # primer bias filter
            if sMtConsByDir['F'] >= 200 and sMtConsByDir['R'] >= 200:
               oddsRatioPB = ((sMtConsByDirByBase[origAlt]['F']+0.5)/(sMtConsByDir['F']+0.5)) / ((sMtConsByDirByBase[origAlt]['R']+.5)/(sMtConsByDir['R']+.5))
               oddsRatioPB = round(oddsRatioPB, 2)
               if oddsRatioPB > 10 or oddsRatioPB < 0.1:
                  fltrs.add('PB')
               primerBiasOR = str(oddsRatioPB)

            # low base quality filter
            if vtype == 'SNP' and alt in alleleCnt and alt in lowQReads and alleleCnt[origAlt] > 0:
               bqAlt =  1.0 * lowQReads[origAlt] / alleleCnt[origAlt]
               if bqAlt > 0.4 and vafToVmfRatio >= 0:
                  if lowQPredict(bqAlt, vafToVmfRatio):
                     fltrs.add('LowQ')
            
            # random end position filter
            if vtype == 'SNP':
               # distance to barcode end of the read
               endBase = 20
               # MT side
               refLeEnd = sum(d <= endBase for d in mtSideBcEndPos[origRef])  # number of REF R2 reads with distance <= endBase
               refGtEnd = len(mtSideBcEndPos[origRef]) - refLeEnd         # number of REF R2 reads with distance > endBase
               altLeEnd = sum(d <= endBase for d in mtSideBcEndPos[origAlt])  # number of ALT R2 reads with distance <= endBase
               altGtEnd = len(mtSideBcEndPos[origAlt]) - altLeEnd         # number of ALT R2 reads with distance > endBase
               fisher = scipy.stats.fisher_exact([[refLeEnd, refGtEnd], [altLeEnd, altGtEnd]])
               oddsRatio = fisher[0]
               pvalue = fisher[1]
               if pvalue < 0.001 and oddsRatio < 0.05 and vaf_tmp <= 60.0:
                  fltrs.add('RBCP')

               # primer side
               refLeEnd = sum(d <= endBase for d in primerSideBcEndPos[origRef])  # number of REF R2 reads with distance <= endBase
               refGtEnd = len(primerSideBcEndPos[origRef]) - refLeEnd         # number of REF R2 reads with distance > endBase
               altLeEnd = sum(d <= endBase for d in primerSideBcEndPos[origAlt])  # number of ALT R2 reads with distance <= endBase
               altGtEnd = len(primerSideBcEndPos[origAlt]) - altLeEnd         # number of ALT R2 reads with distance > endBase
               fisher = scipy.stats.fisher_exact([[refLeEnd, refGtEnd], [altLeEnd, altGtEnd]])
               oddsRatio = fisher[0]
               pvalue = fisher[1]
               if pvalue < 0.001 and oddsRatio < 0.05 and vaf_tmp <= 60.0:
                  fltrs.add('RPCP')

               # fixed end position filter
               endBase = primerDist  # distance to primer end of the read
               refLeEnd = sum(d <= endBase for d in primerSidePrimerEndPos[origRef])  # number of REF R2 reads with distance <= endBase
               refGtEnd = len(primerSidePrimerEndPos[origRef]) - refLeEnd         # number of REF R2 reads with distance > endBase
               altLeEnd = sum(d <= endBase for d in primerSidePrimerEndPos[origAlt])  # number of ALT R2 reads with distance <= endBase
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


         firstAlt = False

         # output metrics for each non-reference allele with >= 3 UMIs; If none, output the one with most UMI
         # final FILTER to output
         fltrFinal = 'PASS' if len(fltrs) == 0 else ';'.join(list(fltrs)) 
         # read-based variant allele fraction (VAF)
         frac_alt = str(round((100.0 * alleleCnt[origAlt] / cvg),3)) if cvg > 0 else '0.0'  # based on all reads, including the excluded reads
         # UMI-based variant allele fraction (VMF)
         vmf = str(round((100.0 * sMtConsByBase[origAlt] / sMtCons),5)) if sMtCons > 0 else '.'  
         # UMI-based VMF for each strand
         vmfForward = str(round((100.0 * sMtConsByDirByBase[origAlt]['F'] / sMtConsByDir['F']),3)) if sMtConsByDir['F'] > 0 else '.'
         vmfReverse = str(round((100.0 * sMtConsByDirByBase[origAlt]['R'] / sMtConsByDir['R']),3)) if sMtConsByDir['R'] > 0 else '.'
         # UMI count for A,C,G,T
         sMTs = [str(sMtConsByBase['A']), str(sMtConsByBase['T']), str(sMtConsByBase['G']), str(sMtConsByBase['C'])]
         # proportion of <Q20 reads
         pLowQ = str(round(bqAlt,2)) if bqAlt >= 0 else 'NA'
         # type of repetitive region
         repTypeFinal = ';'.join(list(repTypeSet)) if len(repTypeSet) >= 1 else 'NA'

         out_long_list = [chrom, pos, ref, alt, vtype, str(sMtCons), str(sMtConsByDir['F']), str(sMtConsByDir['R']), str(sMtConsByBase[origAlt]), str(sMtConsByDirByBase[origAlt]['F']), str(sMtConsByDirByBase[origAlt]['R']),  vmf, vmfForward, vmfReverse, str(alleleCnt[origAlt]), frac_alt, str(refForPrimer), str(refRevPrimer), primerBiasOR, pLowQ, str(hqUmiEff), str(allUmiEff), str(refRppUmiMean), str(altRppUmiMean), str(RppEffSize), repTypeFinal, hpInfo, srInfo, repInfo, str(cvg), str(allFrag), str(allMT), str(usedFrag)] + sMTs + [fltrFinal]
         out_long_allele = '\t'.join(out_long_list) + '\n'
         out_long += out_long_allele

         altCnt += 1
         if altCnt >= maxAltAllele:
            break

   samfile.close()
   refseq.close()
   
   return (out_long, out_bkg)

#----------------------------------------------------------------------------------------------
# global for argument parsing (hack that works when calling from either command line or pipeline)
#------------------------------------------------------------------------------------------------
parser = None
def argParseInit():  # this is done inside a function because multiprocessing module imports the script
   global parser
   parser = argparse.ArgumentParser(description='Variant calling using molecular barcodes')
   parser.add_argument('--runPath', default=None, help='path to working directory')
   parser.add_argument('--bedTarget', default=None, help='BED file')
   parser.add_argument('--bamFile', default=None, help='BAM file')
   parser.add_argument('--outPrefix', default=None, help='file name prefix')
   parser.add_argument('--nCPU', type=int, default=1, help='number of CPU to use in parallel')
   parser.add_argument('--minBQ', type=int, default=25, help='minimum base quality allowed for analysis')
   parser.add_argument('--minMQ', type=int, default=50, help="minimum mapping quality allowed for analysis. If the bam is tagged with its mate's mapq, then the minimum of the R1 and R2 mapq will be used for comparison, if not each read is compared independently.")
   parser.add_argument('--hpLen', type=int, default=10, help='Minimum length for homopolymers')
   parser.add_argument('--mismatchThr', type=float, default=6.0, help='average number of mismatches per 100 bases allowed')
   parser.add_argument('--primerDist', type=int, default=2, help='filter variants that are within X bases to primer')
   parser.add_argument('--mtThreshold', type=float, default=0.8, help='threshold on read proportion to determine MT level consensus')
   parser.add_argument('--rpb', type=float, default=0.0, help='mean read pairs per UMI; default at 0 and will be calculated')
   parser.add_argument('--isRna', action = 'store_true', help='RNAseq varinat calling only; default is DNAseq')
   parser.add_argument('--primerSide', type=int, default=1, help='read end that includes the primer; default is 1')
   parser.add_argument('--minAltUMI', type=int, default=3, help='minimum requirement of ALT UMIs; default is 3')
   parser.add_argument('--maxAltAllele', type=int, default=2, help='maximum ALT alleles that meet minAltUMI to be reported; default is 2 (tri-allelic variants)')
   parser.add_argument('--refGenome',type=str,help='Path to the reference fasta file')
   parser.add_argument('--repBed',type=str,help='Path to the simpleRepeat bed file')
   parser.add_argument('--srBed',type=str,help='Path to the full repeat bed file')
   parser.add_argument('--ds', type=int, default=10000, help='down sample if number of UMIs greater than this value')

#--------------------------------------------------------------------------------------
# main function
#--------------------------------------------------------------------------------------
def main(args):
   # log run start
   timeStart = datetime.datetime.now()
   print("started at " + str(timeStart))
   
   # if argument parser global not assigned yet, initialize it
   if parser == None:
      argParseInit()
   
   # get arguments passed in via a lambda object (e.g. from upstream pipeline)
   if type(args) is not argparse.Namespace:
      argsList = []
      for argName, argVal in args.iteritems():
         argsList.append("--{0}={1}".format(argName, argVal))
      args = parser.parse_args(argsList)
   
   for argName, argVal in vars(args).iteritems():
      print(argName, argVal)
      
   # change working directory to runDir and make output directories
   if args.runPath != None:
      os.chdir(args.runPath)
   # make /intermediate directory to keep the long output
   if not os.path.exists('intermediate'):
      os.makedirs('intermediate')

   # intersect repeats and target regions
   if args.isRna: 
      seqType = 'rna'
      findHpLen = args.hpLen
   else:   
      seqType = 'dna'
      findHpLen = 6
   subprocess.check_call('python ' + homopolymerCode + ' ' + args.bedTarget + ' hp.roi.bed ' + str(findHpLen) + ' ' + args.refGenome + ' ' + seqType, shell=True)
   subprocess.check_call('bedtools intersect -a ' + args.repBed + ' -b ' + args.bedTarget + ' | bedtools sort -i > rep.roi.bed', shell=True)
   subprocess.check_call('bedtools intersect -a ' + args.srBed +  ' -b ' + args.bedTarget + ' | bedtools sort -i > sr.roi.bed', shell=True)

   # homopolymer 
   hpRegion = defaultdict(list)
   with open('hp.roi.bed','r') as IN:
      for line in IN:   
         [chrom, regionStart, regionEnd, repType, totalLen, unitLen, repLen, repBase] = line.strip().split()
         hpRegion[chrom].append([regionStart, regionEnd, repType, totalLen, unitLen, repLen])

   # tandem repeat 
   repRegion = defaultdict(list)
   with open('rep.roi.bed','r') as IN:
      for line in IN:
         lineList = line.strip().split()
         chrom = lineList[0]
         regionStart = lineList[1]
         regionEnd = lineList[2]
         unitLen = lineList[4]
         repLen = lineList[5]
         try:
            unitLen_num = float(unitLen)
         except ValueError:
            continue
         try:
            repLen_num = float(repLen)
         except ValueError:
            continue

         if args.isRna:
            totalLen = int(regionEnd) - int(regionStart)
            if totalLen < args.hpLen:
               continue
            repLen = str(totalLen / unitLen_num)
            totalLen = str(totalLen)
         else:
            totalLen = str(unitLen_num * repLen_num)
            
         repBase = lineList[-1]
         repType = 'RepT'
         repRegion[chrom].append([regionStart, regionEnd, repType, totalLen, unitLen, repLen])

   # simple repeat, low complexity, satelite 
   srRegion = defaultdict(list)
   with open('sr.roi.bed','r') as IN:
      for line in IN:
         [chrom, regionStart, regionEnd, repType, totalLen, unitLen, repLen, repBase] = line.strip().split()
         if repType == 'Simple_repeat':
            repType = 'RepS'
         elif repType == 'Low_complexity':
            repType = 'LowC'
         elif repType == 'Satellite':
            repType = 'SL'
         else:
            repType = 'Other_Repeat'
         
         if args.isRna:
            totalLen = int(regionEnd) - int(regionStart)
            if totalLen < args.hpLen:
               continue
            try:
               unitLen_num = float(unitLen)
               repLen = str(totalLen / unitLen_num)
            except ValueError:
               pass
            totalLen = str(totalLen)
         
         srRegion[chrom].append([regionStart, regionEnd, repType, totalLen, unitLen, repLen])

   # read in bed file and create a list of positions, annotated with repetitive region
   locList = []
   with open(args.bedTarget,'r') as IN:
      for line in IN:
         if line.startswith('track name='):
            continue
         lineList = line.strip().split('\t')
         chrom = lineList[0]
         regionStart = int(lineList[1]) + 1   # target region starts from 1-base after 
         regionEnd = lineList[2]

         pos = regionStart
         lineEnd = False

         while not lineEnd:
            (hpInfo, srInfo, repInfo) = ('.', '.', '.')
            repTypeSet = set()
            # check if the site is in homopolymer region (not including 1 base before) 
            for (regionStart_tmp, regionEnd_tmp, repType_tmp, totalLen_tmp, unitLen_tmp, repLen_tmp) in hpRegion[chrom]:
               if pos >= int(regionStart_tmp) - 0 and pos <= int(regionEnd_tmp):
                  repTypeSet.add(repType_tmp)
                  hpInfo = ';'.join([chrom, regionStart_tmp, regionEnd_tmp, totalLen_tmp])
                  break

            # check if the site is in other repeats region (including 1 base before) 
            for (regionStart_tmp, regionEnd_tmp, repType_tmp, totalLen_tmp, unitLen_tmp, repLen_tmp) in srRegion[chrom]:
               if pos >= int(regionStart_tmp) - 1 and pos <= int(regionEnd_tmp):
                  repTypeSet.add(repType_tmp)
                  srInfo = ';'.join([chrom, regionStart_tmp, regionEnd_tmp, totalLen_tmp, unitLen_tmp, repLen_tmp])
                  break

            for [regionStart_tmp, regionEnd_tmp, repType_tmp, totalLen_tmp, unitLen_tmp, repLen_tmp] in repRegion[chrom]:
               if pos >= int(regionStart_tmp) - 1 and pos <= int(regionEnd_tmp):
                  repTypeSet.add(repType_tmp)
                  repInfo = ';'.join([chrom, regionStart_tmp, regionEnd_tmp, totalLen_tmp, unitLen_tmp, repLen_tmp])
                  break

            repType = 'NA' if len(repTypeSet) == 0 else ';'.join(list(repTypeSet))
            locList.append((chrom, str(pos), repType, hpInfo, srInfo, repInfo))

            if str(pos) == regionEnd:
               lineEnd = True
            else:
               pos += 1

   # calculate rpb if args.rpb = 0
   if args.rpb == 0.0:
      rpb = getMeanRpb(args.bamFile) 
      print("rpb = " + str(round(rpb,1)) + ", computed by smCounter")
   else:
      rpb = args.rpb
      print("rpb = " + str(round(rpb,1)) + ", given by user")
      
   # set primer side
   primerSide = 'R1' if args.primerSide == 1 else 'R2'

   #----- loop over locs
   # prepare to save to disk
   outfile_long = open('intermediate/nopval.' + args.outPrefix + '.VariantList.long.txt', 'w')
   bkgFileName = 'intermediate/bkg.' + args.outPrefix + '.txt'
   outfile_bkg = open(bkgFileName, 'w')

   header_1 = ['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'sUMT', 'sForUMT', 'sRevUMT', 'sVMT', 'sForVMT', 'sRevVMT', 'sVMF', 'sForVMF', 'sRevVMF', 'VDP', 'VAF', 'RefForPrimer', 'RefRevPrimer', 'primerOR', 'pLowQ', 'hqUmiEff', 'allUmiEff', 'refMeanRpb', 'altMeanRpb', 'rpbEffectSize', 'repType', 'hpInfo', 'simpleRepeatInfo', 'tandemRepeatInfo', 'DP', 'FR', 'MT', 'UFR', 'sUMT_A', 'sUMT_T', 'sUMT_G', 'sUMT_C', 'FILTER']
   header_2 = ['CHROM', 'POS', 'REF', 'A/G', 'G/A', 'C/T', 'T/C', 'A/C', 'C/A', 'A/T', 'T/A', 'C/G', 'G/C', 'G/T', 'T/G', 'negStrand', 'posStrand', 'AllSMT' ]

   outfile_long.write('\t'.join(header_1) + '\n')
   outfile_bkg.write('\t'.join(header_2) + '\n')

   # process in chunks (for more granular logging)
   locChunkLen = 1000
   print "chunk size = " + str(locChunkLen)
   locLen = len(locList)
   locChunkCount = locLen / locChunkLen
   if locLen % locChunkLen:
      locChunkCount += 1
   
   for idx in range(locChunkCount):
      # this chunks
      idxStart = idx * locChunkLen
      idxEnd = min(idxStart + locChunkLen, locLen)
      locChunk = locList[idxStart:idxEnd]
      
      # run Python multiprocessing module
      pool = multiprocessing.Pool(processes=args.nCPU)
      results = [pool.apply_async(vc_wrapper, args=(args.bamFile, x[0], x[1], x[2], x[3], x[4], x[5], args.minBQ, args.minMQ, args.hpLen, args.mismatchThr, args.primerDist, args.mtThreshold, rpb, primerSide, args.refGenome, args.minAltUMI, args.maxAltAllele, args.isRna, args.ds)) for x in locChunk]
      
      # clear finished pool
      pool.close()
      pool.join()
      
      # get results - a list of tuples of 2 strings
      output = [p.get() for p in results]
      
      # check for exceptions thrown by vc()
      for idx in range(len(output)):
         line,bg = output[idx]
         if line.startswith("Exception thrown!"):
            print(line)
            raise Exception("Exception thrown in vc() at location: " + str(locChunk[idx]))
      
      # write output and bkg files to disk, for current chunk
      for (vcOutline, bkgOutline) in output:
         outfile_long.write(vcOutline)
         outfile_bkg.write(bkgOutline)
         
      # log some info
      memInfo = [x for x in subprocess.check_output('free -m', shell = True).split('\n')[1].split(' ') if x.isalnum()]
      memTotal = int(memInfo[0])
      memUsed = int(memInfo[1])
      memUsedPct = round(100.0 * memUsed / memTotal, 1)
      posStart = locChunk[0][0] + ":" + locChunk[0][1]
      posEnd = locChunk[-1][0] + ":" + locChunk[-1][1]
      print str(datetime.datetime.now()) + "\t" + posStart + "\t" + posEnd + "\t" + str(memUsed) + "\t" + str(memUsedPct)
      
      del pool, results, output

   outfile_long.close()
   outfile_bkg.close()

   # calculate p-value
   print("Calculating p-values at " + str(datetime.datetime.now()) + "\n")
   outfile1 = 'intermediate/nopval.' + args.outPrefix + '.VariantList.long.txt'

   outfile2 = 'intermediate/' + args.outPrefix + '.VariantList.long.txt'
   outfile_lod = 'intermediate/' + args.outPrefix + '.umi_depths.lod.bedgraph'
   pValCmd = ' '.join(['Rscript', pValCode, args.runPath, outfile1, bkgFileName, str(seed), str(nsim), outfile2, outfile_lod, args.outPrefix, str(rpb), str(args.minAltUMI)])
   subprocess.check_call(pValCmd, shell=True)
   print("completed p-values at " + str(datetime.datetime.now()) + "\n")

   # make VCFs
   vcfCmd = ' '.join(['python', vcfCode, args.runPath, outfile2, args.outPrefix])
   subprocess.check_call(vcfCmd, shell=True)

   # remove intermediate files
   os.remove('hp.roi.bed')
   os.remove('rep.roi.bed')
   os.remove('sr.roi.bed')
   os.remove(outfile1)

   # log run completion
   timeEnd = datetime.datetime.now()
   print("completed running at " + str(timeEnd) + "\n")
   print("total time: "+ str(timeEnd-timeStart) + "\n")  
   
#pythonism to run from the command line
#----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
if __name__ == "__main__":
   # init the argumet parser
   argParseInit()
   
   # get command line arguments
   args = parser.parse_args()
   
   # call main program   
   main(args)
