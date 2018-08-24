import os
import sys
import datetime
import subprocess
import operator
import multiprocessing
from collections import defaultdict
import random
import numpy
import string
import logging
import traceback
# modules from this project
import filters
# 3rd party modules
import pysam

# Set up some rudimentary logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
logger.addHandler(ch)

minTotalUMI = 5
lowBqThr = 20
endBase = 20
mtTag = "Mi"
mqTag = "MQ"
tagSeparator = "-"
primerTag = "pr"
_num_cols_ = 38 ## Number of columns in out_long returned by the vc() function of smCounter
maxDnaReadDepth = 1000000000
downsamplePileupStackThr = 10 ** 5
_base_complement_ = string.maketrans("ACTG","TGAC")

#-------------------------------------------------------------------------------------
# create variables
#-------------------------------------------------------------------------------------
def defineVariables():
   sMtCons, smtSNP = 0, 0
   sMtConsByBase = defaultdict(int)
   sMtConsByDir  = defaultdict(int)
   sMtConsByDirByBase = defaultdict(lambda: defaultdict(int))
   strands = defaultdict(int)
   subTypeCnt = defaultdict(int)
   hqAgree = defaultdict(int)
   hqDisagree = defaultdict(int)
   allAgree = defaultdict(int)
   allDisagree = defaultdict(int)
   rpbCnt = defaultdict(list)
   sMtConsByBase['A'] = 0
   sMtConsByBase['T'] = 0
   sMtConsByBase['G'] = 0
   sMtConsByBase['C'] = 0
   out_long = ''
   
   return(sMtCons, smtSNP, sMtConsByBase, sMtConsByDir, sMtConsByDirByBase, strands, subTypeCnt, hqAgree, hqDisagree, allAgree, allDisagree, rpbCnt, sMtConsByBase, out_long)

#-------------------------------------------------------------------------------------
# get the reference base
#-------------------------------------------------------------------------------------
def getRef(refseq, chrom, pos):
   origRef = refseq.fetch(reference=chrom, start=int(pos)-1, end=int(pos))
   origRef = origRef.upper()
   return(origRef)
  
#-------------------------------------------------------------------------------------
# condition to drop reads
#-------------------------------------------------------------------------------------
def dropRead(pileupRead, pos, cigar):
   # check if position not on a gap (N or intron in RNAseq)
   isDrop = False   
   alignLen = int(pos) - pileupRead.alignment.pos
   # first case: perhaps outside the whole read
   if alignLen > sum([value if op in [0, 3] else 0 for (op, value) in cigar]):
      isDrop = True
   # second case: may lie on any segments
   for (op, value) in cigar:
      if op > 3:
         continue
      if alignLen <= value:
         if op == 3:
            isDrop = True
         break
      elif op in [0, 3]:
         alignLen -= value
   
   return(isDrop)

#-------------------------------------------------------------------------------------
# get some basic information. NOTE: BC depends on the type of input BAM file
#-------------------------------------------------------------------------------------   
def getBasicInfo(pileupRead, bamType,readid):
   # if the input BAM is consensused, use read ID as UMI barcode
   BC = pileupRead.alignment.get_tag(mtTag) if bamType == 'raw' else readid
   # cigar
   cigar = pileupRead.alignment.cigar
   # alignment positions
   astart = pileupRead.alignment.reference_start
   aend = pileupRead.alignment.reference_end
   return(BC, cigar,astart,aend)

#-------------------------------------------------------------------------------------
# condition to drop reads
#-------------------------------------------------------------------------------------   
def getBaseAndBq(pileupRead, refseq, chrom, pos, minBQ):
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

   return(base, bq)

#-------------------------------------------------------------------------------------
# check if a read is high quality and can be included in the bcDictHq
#-------------------------------------------------------------------------------------   
def hqRead(pileupRead,cigar,minMQ,mismatchThr):
   # mapping quality filter - both R1 and R2 need to meet the minimum mapQ
   mq = pileupRead.alignment.mapping_quality
   minMQPass = True   
   try:   # get mapq of mate
      mateMq = pileupRead.alignment.get_tag(mqTag)
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
         
   # check if there are too many mismatches, excluding indel            
   NM = 0 # get NM tag
   allTags = pileupRead.alignment.tags
   for (tag, value) in allTags:
      if tag == 'NM':
         NM = value
         break
   nIndel = 0 # count number of INDELs in the read sequence
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
   
   # overall condition for high quality read
   incCond = minMQPass and mismatchPer100b <= mismatchThr

   return(leftSP,incCond)

#-------------------------------------------------------------------------------------
# check if the read covers the entire homopolymer stretch
#-------------------------------------------------------------------------------------
def isHPCovered(astart,aend,hpInfo):
   if hpInfo == '.':
      hpCovered = True
   else:
      hpChrom, hpStart, hpEnd, totalHpLen, realL, realR = hpInfo.strip().split(';')
      if astart < int(hpStart) - 1 and aend > int(hpEnd) + 1:
         hpCovered = True
      else:
         hpCovered = False
         
   return hpCovered

#-------------------------------------------------------------------------------------
# update read level metrics
#-------------------------------------------------------------------------------------  
def updateReadMetrics(pileupRead, base, bq, incCond, pairOrder, leftSP, mtSide, primerSide, alleleCnt, forwardCnt, reverseCnt, lowQReads, mtSideBcEndPos, primerSideBcEndPos, primerSidePrimerEndPos, cvg):
   # update metrics for all types of positions
   alleleCnt[base] += 1
   strand = 'Reverse' if pileupRead.alignment.is_reverse else 'Forward' # +/- strand and update read depth by strand
   if strand == 'Reverse':
      reverseCnt[base] += 1
   else:
      forwardCnt[base] += 1
   
   # update metrics for potential SNPs only; metrics will be used for filters
   if pileupRead.indel == 0 and not pileupRead.is_del:
      # count the number of low quality reads (less than Q20 by default) for each base
      if bq < lowBqThr:
         lowQReads[base] += 1
      if pairOrder == mtSide:
         # distance to the random (UMI) end
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
   
   return(alleleCnt, forwardCnt, reverseCnt, lowQReads, mtSideBcEndPos, primerSideBcEndPos, primerSidePrimerEndPos, cvg)

#-------------------------------------------------------------------------------------
# Group reads by UMI and update some metrics 
#-------------------------------------------------------------------------------------  
def groupByUMI(readid, BC, base, pairOrder, usedFrag, allFrag, incCond, hpCovered, allBcDict, bcDictHq, bcDictHqBase, bcDictAll, concordPairCnt, discordPairCnt):
   # count total number of fragments and UMIs
   if readid not in allBcDict[BC]:
      allFrag+=1 # total fragments
      allBcDict[BC].add(readid)

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
   
   return(allBcDict, bcDictHq, bcDictAll, bcDictHqBase, concordPairCnt, discordPairCnt, allFrag, usedFrag)

#-------------------------------------------------------------------------------------
# call pysam function to pile up an interval
#-------------------------------------------------------------------------------------
def pileup(bamName,chrom,start,end):
   samfile = pysam.AlignmentFile(bamName,"rb")
   current_pos = int(start)
   for p in samfile.pileup(region = chrom + ":" + start + ":" + end, truncate=True,max_depth = maxDnaReadDepth, stepper="nofilter"):
      ref_pos = p.pos+1
      while True:
         if ref_pos < current_pos:
            raise Exception("pysam returned out of bounds coordinate !")
         elif ref_pos > current_pos:
            yield None
            current_pos+=1
         else:
            yield p
            current_pos+=1
            break
      
   samfile.close()
   
#-------------------------------------------------------------------------------------
# pile up reads and group by UMI; some metrics are updated here
#-------------------------------------------------------------------------------------
def pileupAndGroupByUMI(bamName, bamType, chrom, pos, repType, hpInfo, minBQ, minMQ, hpLen, mismatchThr, primerDist, mtThreshold, rpb, primerSide, refseq, minAltUMI, maxAltAllele, isRna, read_pileup, hqCache, infoCache):

   # define variables
   cvg, usedFrag, allFrag = 0, 0, 0
   lowQReads = defaultdict(int)
   alleleCnt = defaultdict(int)
   forwardCnt = defaultdict(int)
   reverseCnt = defaultdict(int)
   concordPairCnt = defaultdict(int)
   discordPairCnt = defaultdict(int)
   mtSideBcEndPos = defaultdict(list)
   primerSideBcEndPos = defaultdict(list)
   primerSidePrimerEndPos = defaultdict(list)
   allBcDict = defaultdict(set)
   bcDictHqBase = defaultdict(lambda:defaultdict(int))
   bcDictAll = defaultdict(lambda:defaultdict(int))
   bcDictHq = defaultdict(lambda: defaultdict(list))

   mtSide = 'R1' if primerSide == 'R2' else 'R2'
  
   # check pielup size, downsample if above threshold; for RNA only
   pileupStackSize = read_pileup.n
   downsamplePileup = True if isRna and pileupStackSize > downsamplePileupStackThr else False
   random.seed(pos)

   # iterate over pileup reads and group by UMI
   for pileupRead in read_pileup.pileups:
      # drop reads randomly; for RNA only
      if downsamplePileup and random.randint(1, pileupStackSize) > downsamplePileupStackThr:
         continue

      # basic information that will be used in subsequent functions; NOTE: if the input BAM is consensused, use read ID as UMI barcode
      readid = pileupRead.alignment.query_name
      pairOrder = 'R1' if pileupRead.alignment.is_read1 else 'R2'
      key = readid + '_' + pairOrder
      if key in infoCache:
         BC,cigar,astart,aend = infoCache[key]
      else:
         BC,cigar,astart,aend = getBasicInfo(pileupRead, bamType,readid)
         infoCache[key] = (BC,cigar,astart,aend)

      # check if read should be dropped
      if dropRead(pileupRead, pos, cigar):
         continue

      # retrive base and base-quality, and re-format the base
      base, bq = getBaseAndBq(pileupRead, refseq, chrom, pos, minBQ)

      # check if the read is high quality
      hpCovered = isHPCovered(astart,aend,hpInfo)
      if key in hqCache:
         leftSP,incCondTemp = hqCache[key]
         incCond = incCondTemp and bq >= minBQ and hpCovered
      else:
         leftSP, incCondTemp = hqRead(pileupRead,cigar,minMQ,mismatchThr)
         incCond = incCondTemp and bq >= minBQ and hpCovered
         hqCache[key] = (leftSP,incCondTemp)

      # update read-level metrics
      alleleCnt, forwardCnt, reverseCnt, lowQReads, mtSideBcEndPos, primerSideBcEndPos, primerSidePrimerEndPos, cvg = updateReadMetrics(pileupRead, base, bq, incCond, pairOrder, leftSP, mtSide, primerSide, alleleCnt, forwardCnt, reverseCnt, lowQReads, mtSideBcEndPos, primerSideBcEndPos, primerSidePrimerEndPos, cvg)

      # group reads by UMIs
      allBcDict, bcDictHq, bcDictAll, bcDictHqBase, concordPairCnt, discordPairCnt, allFrag, usedFrag = groupByUMI(readid, BC, base, pairOrder, usedFrag, allFrag, incCond, hpCovered, allBcDict, bcDictHq, bcDictHqBase, bcDictAll, concordPairCnt, discordPairCnt)

   # output variables
   return(alleleCnt, forwardCnt, reverseCnt, lowQReads, mtSideBcEndPos, primerSideBcEndPos, primerSidePrimerEndPos, cvg, allBcDict, bcDictHq, bcDictAll, bcDictHqBase, concordPairCnt, discordPairCnt, allFrag, usedFrag,hqCache,infoCache)

#-------------------------------------------------------------------------------------
# gradually drop singleton UMIs, depending on rpb and input mode
#-------------------------------------------------------------------------------------
def dropSingleton(rpb, bcDictHq, pos, ds, cvg, allMT, isRna, bamType):
   singleUMIs = set()
   pairedUMIs = set()
   bcToKeep = []

   if isRna:
      rpb = cvg / float(allMT) if allMT > 0 else 1.0
   
   # rpb < 2 or consensused BAM input: no UMI is dropped
   if rpb < 2.0 or bamType == 'consensus':
      bcToKeep = bcDictHq.keys()

   # 2 <= rpb < 3: gradually and randomly drop singleton UMIs 
   elif rpb >= 2.0 and rpb < 3.0:
      # set seed to be the genome position
      random.seed(pos)
      # count the numbers of paired and unpaired singleton UMIs; 
      pctToDrop = rpb - 2.0
      for bc in bcDictHq:
         readPairsInBc = len(bcDictHq[bc])
         if readPairsInBc == 1:
            readid = bcDictHq[bc].keys()[0]
            if bcDictHq[bc][readid][1] == 'Paired':
               pairedUMIs.add(bc)
            else:
               singleUMIs.add(bc)
      # total number of singleton UMIs
      pairedCnt = len(pairedUMIs)
      singleCnt = len(singleUMIs)
      oneReadMtCnt = pairedCnt + singleCnt
      # number of singleton UMIs to drop
      numToDrop = int(round(pctToDrop * oneReadMtCnt))
      # Decide which singleton UMIs to drop -- paired reads are kept with priority
      if numToDrop <= singleCnt:
         oneReadMtToDrop = set(random.sample(singleUMIs, numToDrop))
      else:
         pairsToDrop = set(random.sample(pairedUMIs, numToDrop - singleCnt))
         oneReadMtToDrop = singleUMIs.union(pairsToDrop)
      # drop singleton UMIs
      bcToKeep = list(set(bcDictHq.keys()).difference(oneReadMtToDrop))      

   # rpb >= 3: drop all singleton UMIs;
   else:
      bcToKeep = [bc for bc in bcDictHq.iterkeys() if len(bcDictHq[bc]) >= 2]
   
   # additional downsample for RNA-seq data only
   if isRna and len(bcToKeep) > ds:
      random.seed(pos)
      bcToKeep = random.sample(bcToKeep, ds)
      
   # output variables
   return(bcToKeep)

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
# consensus for each UMI family; 
#-------------------------------------------------------------------------------------
def consensus(bcDictHqBase, bcDictAll, bc, mtThreshold, bamType):
   tmpHqBc = bcDictHqBase[bc]
   tmpAllBc = bcDictAll[bc]
   
   if bamType == 'raw':
      consHq = consHqMT(tmpHqBc, mtThreshold)
      consAll = consAllMT(tmpAllBc, mtThreshold)
      cons = consHq if consHq == consAll else ''
   else:
      if len(tmpHqBc) == 2 and 'all' in tmpHqBc: 
         del tmpHqBc['all']
         cons = tmpHqBc.keys()[0]
      else:
         cons = ''
       
   return(cons)

#-------------------------------------------------------------------------------------
# update the UMI metrics
#-------------------------------------------------------------------------------------
def updateUmiMetrics(bc, bcDictHqBase, cons, hqAgree, hqDisagree, bcDictAll, allAgree, allDisagree, origRef, sMtCons, sMtConsByBase, sMtConsByDir, sMtConsByDirByBase, rpbCnt, subTypeCnt, smtSNP, strands):
   # primer ID and direction
   bcSplit = bc.split(tagSeparator)
   primerDirCode = bcSplit[1]
   primerDirection = 'F' if primerDirCode == '0' else 'R' # 0 means the primer was priming the forward strand, 1 means priming the reverse strand

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
         basePair = origRef + '/' + cons if primerDirCode == '0' else origRef.translate(_base_complement_) + '/' + cons.translate(_base_complement_)
         subTypeCnt[basePair] += 1
         smtSNP += 1
         strands[primerDirection] += 1

   return(hqAgree, hqDisagree, allAgree, allDisagree, sMtCons, sMtConsByBase, sMtConsByDir, sMtConsByDirByBase, rpbCnt, subTypeCnt, smtSNP, strands)

#-------------------------------------------------------------------------------------
# save the background error profile
#-------------------------------------------------------------------------------------
def outbkg(chrom, pos, origRef, subTypeCnt, strands, smtSNP):
   bkgErrList = [chrom, pos, origRef, str(subTypeCnt['A/G']), str(subTypeCnt['G/A']), str(subTypeCnt['C/T']), str(subTypeCnt['T/C']), str(subTypeCnt['A/C']), str(subTypeCnt['C/A']), str(subTypeCnt['A/T']), str(subTypeCnt['T/A']), str(subTypeCnt['C/G']), str(subTypeCnt['G/C']), str(subTypeCnt['G/T']), str(subTypeCnt['T/G']), str(strands['F']), str(strands['R']), str(smtSNP)]
      
   out_bkg = '\t'.join(bkgErrList) + '\n'

   # output variables
   return(out_bkg)

#-------------------------------------------------------------------------------------
# reset variant type, reference base, variant base 
#-------------------------------------------------------------------------------------
def setRefAltType(origRef, origAlt):
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

   return(ref, alt, vtype)

#-------------------------------------------------------------------------------------
# compute UMI efficiency metrics
#-------------------------------------------------------------------------------------
def umiEfficiency(hqAgree, hqDisagree, allAgree, allDisagree, origRef, origAlt, rpbCnt, alleleCnt, sMtConsByBase, cvg, sMtCons, vaf_tmp, vmf_tmp):
   hqRcAgree = hqAgree[origAlt] 
   hqRcTotal = hqRcAgree + hqDisagree[origAlt] 
   hqUmiEff = 1.0 * hqRcAgree / hqRcTotal if hqRcTotal > 0 else 0.0
   
   allRcAgree = allAgree[origAlt] 
   allRcTotal = allRcAgree + allDisagree[origAlt] 
   allUmiEff = 1.0 * allRcAgree / allRcTotal if allRcTotal > 0 else 0.0

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

   vafToVmfRatio = 1.0 * vaf_tmp / vmf_tmp if vmf_tmp > 0 else -1.0

   return (vafToVmfRatio, hqUmiEff, allUmiEff, refRppUmiMean, altRppUmiMean, RppEffSize)


#-------------------------------------------------------------------------------------
# output to *.all.txt file
#-------------------------------------------------------------------------------------
def outlong(out_long, chrom, pos, ref, alt, vtype, origRef, origAlt, sMtCons, sMtConsByDir, sMtConsByBase, sMtConsByDirByBase, alleleCnt, primerBiasOR, hqUmiEff, allUmiEff, refRppUmiMean, altRppUmiMean, RppEffSize, repTypeSet, bqAlt, hpInfo, srInfo, repInfo, cvg, allFrag, allBcDict, usedFrag, fltrs):
   # total number of MT, fragments, reads, including those dropped from analysis
   allMT = len(allBcDict)
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
   # UMI counts by primer direction
   refForPrimer = sMtConsByDirByBase[origRef]['F']
   refRevPrimer = sMtConsByDirByBase[origRef]['R']
   altForPrimer = sMtConsByDirByBase[origAlt]['F']
   altRevPrimer = sMtConsByDirByBase[origAlt]['R']

   # round hqUmiEff and allUmiEff in the end
   hqUmiEff = round(hqUmiEff, 3)
   allUmiEff = round(allUmiEff, 3)

   out_long_list = [chrom, pos, ref, alt, vtype, str(sMtCons), str(sMtConsByDir['F']), str(sMtConsByDir['R']), str(sMtConsByBase[origAlt]), str(sMtConsByDirByBase[origAlt]['F']), str(sMtConsByDirByBase[origAlt]['R']),  vmf, vmfForward, vmfReverse, str(alleleCnt[origAlt]), frac_alt, str(refForPrimer), str(refRevPrimer), primerBiasOR, pLowQ, str(hqUmiEff), str(allUmiEff), str(refRppUmiMean), str(altRppUmiMean), str(RppEffSize), repTypeFinal, hpInfo, srInfo, repInfo, str(cvg), str(allFrag), str(allMT), str(usedFrag)] + sMTs + [fltrFinal]
   out_long_allele = '\t'.join(out_long_list) + '\n'
   out_long += out_long_allele

   return(out_long)
   
#-------------------------------------------------------------------------------------
# function to call variants
#-------------------------------------------------------------------------------------
def vc(bamName, chrom, pos, repType, hpInfo, srInfo, repInfo, minBQ, minMQ, hpLen, mismatchThr, primerDist, mtThreshold, rpb, primerSide, refseq, minAltUMI, maxAltAllele, isRna, ds, bamType, read_pileup, hqCache, infoCache, chromLength):

   # initiate variables
   sMtCons, smtSNP, sMtConsByBase, sMtConsByDir, sMtConsByDirByBase, strands, subTypeCnt, hqAgree, hqDisagree, allAgree, allDisagree, rpbCnt, sMtConsByBase, out_long = defineVariables()
   
   # find the reference base
   origRef = getRef(refseq,chrom,pos)
   
   # pile up and group reads by UMIs
   alleleCnt, forwardCnt, reverseCnt, lowQReads, mtSideBcEndPos, primerSideBcEndPos, primerSidePrimerEndPos, cvg, allBcDict, bcDictHq, bcDictAll, bcDictHqBase, concordPairCnt, discordPairCnt, allFrag, usedFrag, hqCache, infoCache = pileupAndGroupByUMI(bamName, bamType, chrom, pos, repType, hpInfo, minBQ, minMQ, hpLen, mismatchThr, primerDist, mtThreshold, rpb, primerSide, refseq, minAltUMI, maxAltAllele, isRna, read_pileup, hqCache, infoCache)
   
   # gradually drop singleton UMIs; keep all UMIs for consensused BAM
   allMT = len(allBcDict)
   bcToKeep = dropSingleton(rpb, bcDictHq, pos, ds, cvg, allMT, isRna, bamType)  
   bcToKeepLen = len(bcToKeep)

   # output zeros or blank if the remaining UMIs is lower than the minimum threshold 
   if bcToKeepLen <= minTotalUMI:
      out_long = '\t'.join([chrom, pos, origRef] + ['0'] * (_num_cols_ - 4) + ['LM']) + '\n'
      out_bkg = ''
   else:        
      for bc in bcToKeep:
         # generate consensus
         cons = consensus(bcDictHqBase, bcDictAll, bc, mtThreshold, bamType)
         # update UMI-level metrics
         hqAgree, hqDisagree, allAgree, allDisagree, sMtCons, sMtConsByBase, sMtConsByDir, sMtConsByDirByBase, rpbCnt, subTypeCnt, smtSNP, strands = updateUmiMetrics(bc, bcDictHqBase, cons, hqAgree, hqDisagree, bcDictAll, allAgree, allDisagree, origRef, sMtCons, sMtConsByBase, sMtConsByDir, sMtConsByDirByBase, rpbCnt, subTypeCnt, smtSNP, strands)

      # output the background error profile
      out_bkg = outbkg(chrom, pos, origRef, subTypeCnt, strands, smtSNP)

      sortedList = sorted(sMtConsByBase.items(), key=operator.itemgetter(1), reverse=True)
      firstAlt = True
      altCnt = 0
      repTypeSet0 = set() if repType == 'NA' else set(repType.strip().split(';'))

      # start multi-allelic loop
      for alleleInd in range(len(sortedList)):
         origAlt = sortedList[alleleInd][0]
         maxVMT = sortedList[alleleInd][1]

   # if the current allele has >= 3 UMIs and is not reference, treat as a candidate variant allele
         if origAlt == origRef:
            continue
         if maxVMT < minAltUMI and not firstAlt:
            break

         # reset variant type, reference base, variant base 
         ref, alt, vtype = setRefAltType(origRef, origAlt) 

         # initiate values for filters and output
         primerBiasOR, bqAlt, hqUmiEff, allUmiEff, refRppUmiMean, altRppUmiMean, RppEffSize = 'NA', -1.0, 0.0, 0.0, -1.0, -1.0, -1.0
         fltrs = set()
         repTypeSet = repTypeSet0
       
         if vtype in ['SNP', 'INDEL'] and sMtCons > 0:
            # compute read- and UMI-level allele frequency and proportion of low base quality reads supporting the variant
            vaf_tmp = 100.0 * alleleCnt[origAlt] / cvg if cvg > 0 else 0.0 
            vmf_tmp = 100.0 * sMtConsByBase[origAlt] / sMtCons

            ### common filters for both original and consensused BAM
            # low coverage filter
            fltrs = filters.lm(fltrs, sMtCons) 
            # strand bias and discordant pairs filter
            fltrs = filters.dp_sb(fltrs, origAlt, concordPairCnt, discordPairCnt, reverseCnt, forwardCnt, origRef, vaf_tmp)

            # Initial HP and LowC region filters
            repTypeSet, hpInfo = filters.isHPorLowComp(chrom, pos, hpLen, ref, alt, refseq, repTypeSet, hpInfo, chromLength)

            # SNP-only common filters
            if vtype == 'SNP':
               # random (UMI) end position filters
               fltrs = filters.rbcp(fltrs, endBase, mtSideBcEndPos, origRef, origAlt, vaf_tmp, isRna)
               fltrs = filters.rpcp(fltrs, endBase, primerSideBcEndPos, origRef, origAlt, vaf_tmp, isRna)
               # proportion of low base quality reads
               if origAlt in alleleCnt and origAlt in lowQReads and alleleCnt[origAlt] > 0:
                  bqAlt = 1.0 * lowQReads[origAlt] / alleleCnt[origAlt]
               
            ### original BAM only filters that use UMI efficiency metrics
            if bamType == 'raw':
               # UMI efficiency metrics; only for original BAM
               vafToVmfRatio, hqUmiEff, allUmiEff, refRppUmiMean, altRppUmiMean, RppEffSize = umiEfficiency(hqAgree, hqDisagree, allAgree, allDisagree, origRef, origAlt, rpbCnt, alleleCnt, sMtConsByBase, cvg, sMtCons, vaf_tmp, vmf_tmp)               
               if not isRna:
                  # update HP for indel
                  fltrs = filters.hp4indel(fltrs, repTypeSet, vtype, rpb, hpInfo, vafToVmfRatio, vmf_tmp, hqUmiEff, RppEffSize, altRppUmiMean)             
                  # update other repetitive region filters for SNP and indel, including HP for SNP
                  fltrs = filters.rep4others(fltrs, repTypeSet, vtype, rpb, vafToVmfRatio, hqUmiEff, RppEffSize)
               # primer bias filter
               fltrs, primerBiasOR = filters.pb(fltrs, origAlt, sMtConsByDir, sMtConsByDirByBase)
               if vtype == 'SNP':
                  # low base quality filter
                  fltrs = filters.lowq(fltrs, lowQReads, alleleCnt, origAlt, vafToVmfRatio, bqAlt, isRna)
                  # fixed end (gene specific primers) position filter
                  fltrs = filters.primercp(fltrs, primerDist, primerSidePrimerEndPos, origRef, origAlt, vmf_tmp, hqUmiEff, vafToVmfRatio, RppEffSize, rpb)
              
            ### consensus BAM only filters; use more strict repetitive and primerOR filters; Discard LowQ filter because the base qualities have different meanings
            else:
               fltrs = filters.strict(fltrs, repTypeSet, vtype, vmf_tmp, primerDist, primerSidePrimerEndPos, origRef, origAlt)

         firstAlt = False
         # output metrics for each non-reference allele with >= 3 UMIs; If none, output the one with most UMI
         out_long = outlong(out_long, chrom, pos, ref, alt, vtype, origRef, origAlt, sMtCons, sMtConsByDir, sMtConsByBase, sMtConsByDirByBase, alleleCnt, primerBiasOR, hqUmiEff, allUmiEff, refRppUmiMean, altRppUmiMean, RppEffSize, repTypeSet, bqAlt, hpInfo, srInfo, repInfo, cvg, allFrag, allBcDict, usedFrag, fltrs)
         
         altCnt += 1
         if altCnt >= maxAltAllele:
            break
   
   return (out_long, out_bkg, hqCache, infoCache)

#------------------------------------------------------------------------------------------------
# wrapper function for "vc()" - because Python multiprocessing module does not pass stack trace; from runone/smcounter.py by John Dicarlo
#------------------------------------------------------------------------------------------------
def vc_wrapper(general_args, interval):
   try:
      output = []
      hqCache = {}
      infoCache = {}
      bamName, minBQ, minMQ, hpLen, mismatchThr, primerDist, mtThreshold, rpb, primerSide, refg, minAltUMI, maxAltAllele, isRna, ds, bamType = general_args

      chrom = interval[0][0]
      intervalStartPos = interval[0][1]
      intervalEndPos = interval[-1][1]
      
      refseq = pysam.FastaFile(refg)
      chromLenghts = {}
      for idx in range(len(refseq.lengths)):
         chromLenghts[refseq.references[idx]] = refseq.lengths[idx]
         
      i = 0
      for read_pileup in pileup(bamName, chrom, intervalStartPos, intervalEndPos):
         site = interval[i]
         i+=1
         chrom,pos,repType,hpInfo,srInfo,repInfo = site
         
         if read_pileup is None: # site is not covered at all, pysam simply skips such sites
            origRef = getRef(refseq,chrom,pos)
            out_long = '\t'.join([chrom, pos, origRef] + ["0"] * (_num_cols_ - 4) + ["LM"]) + "\n"
            out = [out_long,""]
            output.append(out)
            continue
         
         temp = [bamName, chrom, pos, repType, hpInfo, srInfo, repInfo, minBQ, minMQ, hpLen, mismatchThr, primerDist, mtThreshold, rpb, primerSide, refseq, minAltUMI, maxAltAllele, isRna, ds, bamType, read_pileup, hqCache, infoCache, chromLenghts[chrom]]
         out_long,out_bkg,hqCache,infoCache = vc(*temp)
         out = [out_long, out_bkg]
         output.append(out)
         
   except Exception as e:
      out = ("Exception thrown!\n" + traceback.format_exc(), "no_bg")
      output.append(out)
      logger.info("Exception thrown in vc() function at genome location : {pos} in interval : {chrom}:{it1}-{it2}".format(pos=pos, chrom=chrom, it1=intervalStartPos, it2=intervalEndPos))
      logger.info(out[0])
      raise Exception(e)
   
   refseq.close()
   logger.info("Finished processing interval {chrom}:{it1}-{it2}".format(chrom=chrom,it1=intervalStartPos,it2=intervalEndPos))
   return output
