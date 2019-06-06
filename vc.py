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

minTotalUmi = 5
lowBqThr = 20
endBase = 20

# Number of columns in detailed output; normal DNA-seq
nColsSin = 38 
# Number of columns in detailed output; duplex-seq
nColsDup = 45
maxDnaReadDepth = 1000000000
downsamplePileupStackThr = 10 ** 5
_base_complement_ = string.maketrans("ACTG", "TGAC")

#-------------------------------------------------------------------------------------
# create variables
#-------------------------------------------------------------------------------------
def defVar():
   sUmiCons, sUmiSnp = 0, 0
   sUmiConsByBase = defaultdict(int)
   sUmiConsByDir  = defaultdict(int)
   sUmiConsByDirByBase = defaultdict(lambda: defaultdict(int))
   sStrands = defaultdict(int)
   sSubTypeCnt = defaultdict(int)
   hqAgree = defaultdict(int)
   hqDisagree = defaultdict(int)
   allAgree = defaultdict(int)
   allDisagree = defaultdict(int)
   rpuCnt = defaultdict(list)
   umiPairDict = defaultdict(set)
   sUmiConsByBase['A'] = 0
   sUmiConsByBase['T'] = 0
   sUmiConsByBase['G'] = 0
   sUmiConsByBase['C'] = 0
   outLineLong = ''

   # placeholder for duplex-specific variables
   dUmiCons = None 
   dUmiSnp = None
   dUmiConsByBase = None
   dUmiConsByDir = None
   dUmiConsByDirByBase = None
   dStrands = None
   dSubTypeCnt = None
   discordDupPairs = None

   return(sUmiCons, sUmiSnp, sUmiConsByBase, sUmiConsByDir, sUmiConsByDirByBase, sStrands, sSubTypeCnt, hqAgree, hqDisagree, allAgree, allDisagree, rpuCnt, umiPairDict, sUmiConsByBase, outLineLong, dUmiCons, dUmiSnp, dUmiConsByBase, dUmiConsByDir, dUmiConsByDirByBase, dStrands, dSubTypeCnt, discordDupPairs)

#-------------------------------------------------------------------------------------
# create additional variables; duplex-seq
#-------------------------------------------------------------------------------------
def dup_defVar():
   sUmiCons, sUmiSnp = 0, 0
   sUmiConsByBase = defaultdict(int)
   sUmiConsByDir  = defaultdict(int)
   sUmiConsByDirByBase = defaultdict(lambda: defaultdict(int))
   sStrands = defaultdict(int)
   sSubTypeCnt = defaultdict(int)
   hqAgree = defaultdict(int)
   hqDisagree = defaultdict(int)
   allAgree = defaultdict(int)
   allDisagree = defaultdict(int)
   rpuCnt = defaultdict(list)
   umiPairDict = defaultdict(set)
   sUmiConsByBase['A'] = 0
   sUmiConsByBase['T'] = 0
   sUmiConsByBase['G'] = 0
   sUmiConsByBase['C'] = 0
   outLineLong = ''
   
   # duplex-seq variables
   dUmiCons, dUmiSnp = 0, 0
   dUmiConsByBase = defaultdict(int)
   dUmiConsByDir  = defaultdict(int)
   dUmiConsByDirByBase = defaultdict(lambda: defaultdict(int))
   dStrands = defaultdict(int)
   dSubTypeCnt = defaultdict(int)
   discordDupPairs = defaultdict(int)
   dUmiConsByBase['A'] = 0
   dUmiConsByBase['T'] = 0
   dUmiConsByBase['G'] = 0
   dUmiConsByBase['C'] = 0
  
   return(sUmiCons, sUmiSnp, sUmiConsByBase, sUmiConsByDir, sUmiConsByDirByBase, sStrands, sSubTypeCnt, hqAgree, hqDisagree, allAgree, allDisagree, rpuCnt, umiPairDict, sUmiConsByBase, outLineLong, dUmiCons, dUmiSnp, dUmiConsByBase, dUmiConsByDir, dUmiConsByDirByBase, dStrands, dSubTypeCnt, discordDupPairs)
   
#-------------------------------------------------------------------------------------
# get the reference base
#-------------------------------------------------------------------------------------
def getRef(refseq, chrom, pos):
   origRef = refseq.fetch(reference = chrom, start = int(pos) - 1, end = int(pos))
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
# get UMI sequence from read
#-------------------------------------------------------------------------------------   
def getUmi(pileupRead, bamType, readid, umiTag, duplexTag = None):
   # UMI sequence not including duplex tags
   umiNoDupTag = pileupRead.alignment.get_tag(umiTag)
   # if the input BAM is consensused, use read ID as umi barcode
   umi = umiNoDupTag if bamType == 'raw' else readid
   dupTag = None

   return(umi, umiNoDupTag, dupTag)

#-------------------------------------------------------------------------------------
# get UMI sequence from read; duplex-seq
#-------------------------------------------------------------------------------------   
def dup_getUmi(pileupRead, bamType, readid, umiTag, duplexTag = None):
   # UMI sequence not including duplex tags
   umiNoDupTag = pileupRead.alignment.get_tag(umiTag)
   dupTag = pileupRead.alignment.get_tag(duplexTag)  
   # complete barcode sequence = duplex tag + UMI, so TT and CC are seperate strands
   umi = dupTag + ':' + umiNoDupTag
   
   return(umi, umiNoDupTag, dupTag)

#-------------------------------------------------------------------------------------
# get some basic information: cigar, start and end of alignment
#-------------------------------------------------------------------------------------   
def getBasicInfo(pileupRead):
   # cigar
   cigar = pileupRead.alignment.cigar
   # alignment positions
   astart = pileupRead.alignment.reference_start
   aend = pileupRead.alignment.reference_end
   
   return(cigar, astart, aend)

#-------------------------------------------------------------------------------------
# condition to drop reads
#-------------------------------------------------------------------------------------   
def getBaseAndBq(pileupRead, refseq, chrom, pos, minBq):
   # check if the site is the beginning of insertion
   if pileupRead.indel > 0:
      site = pileupRead.alignment.query_sequence[pileupRead.query_position]
      inserted = pileupRead.alignment.query_sequence[(pileupRead.query_position + 1) : (pileupRead.query_position + 1 +  pileupRead.indel)]
      base = 'INS|' + site + '|' + site + inserted
      bq = pileupRead.alignment.query_qualities[pileupRead.query_position]
      # if base quality not included in BAM
      if bq == None:
         bq = minBq         
   
   # check if the site is the beginning of deletion
   elif pileupRead.indel < 0:
      site = pileupRead.alignment.query_sequence[pileupRead.query_position]
      deleted = refseq.fetch(reference = chrom, start = int(pos), end = int(pos) + abs(pileupRead.indel))
      deleted = deleted.upper()
      base = 'DEL|' + site + deleted + '|' + site
      bq = pileupRead.alignment.query_qualities[pileupRead.query_position]
      # if base quality not included in BAM
      if bq == None:
         bq = minBq         

   # site is not beginning of any INDEL, but in the middle of a deletion
   elif  pileupRead.is_del:
      base = 'DEL'
      bq = minBq

   # if the site is a regular locus, 
   else: 
      base = pileupRead.alignment.query_sequence[pileupRead.query_position] # note: query_sequence includes soft clipped bases
      bq = pileupRead.alignment.query_qualities[pileupRead.query_position]

   return(base, bq)

#-------------------------------------------------------------------------------------
# check if a read is high quality and can be included in the umiDictHq
#-------------------------------------------------------------------------------------   
def hqRead(pileupRead, cigar, minMq, mismatchThr, mqTag):
   # mapping quality filter - both R1 and R2 need to meet the minimum mapQ
   mq = pileupRead.alignment.mapping_quality
   minMqPass = True   
   try:   # get mapq of mate
      mateMq = pileupRead.alignment.get_tag(mqTag)
      minFragMQ = min(mq,mateMq)
      if minFragMQ < minMq:
         minMqPass = False
   except KeyError: 
      '''
      bam has not been tagged with the mate mapq,
      drop read pairs based on their respective mapqs only
      To note :
      warn user ? or make command line argument more descriptive
      settling on a more descriptive argument for now
      '''
      if mq < minMq:
         minMqPass = False
         
   # check if there are too many mismatches, excluding indel            
   NM = 0 # get NM tag
   allTags = pileupRead.alignment.tags
   for (tag, value) in allTags:
      if tag == 'NM':
         NM = value
         break
   nIndel = 0 # count number of INDELs in the read sequence
   cigarOrder = 1
   leftSp = 0  # soft clipped bases on the left
   rightSp = 0  # soft clipped bases on the right
   for (op, value) in cigar:
      # 1 for insertion
      if op == 1 or op == 2:
         nIndel += value
      if cigarOrder == 1 and op == 4:
         leftSp = value
      if cigarOrder > 1 and op == 4:
         rightSp += value
      cigarOrder += 1

   # Number of mismatches except INDEL, including softcilpped sequences 
   mismatch = max(0, NM - nIndel)
   # read length, including softclip
   readLen = pileupRead.alignment.query_length
   # calculate mismatch per 100 bases
   mismatchPer100b = 100.0 * mismatch / readLen if readLen > 0 else 0.0
   
   # overall condition for high quality read
   incCond = minMqPass and mismatchPer100b <= mismatchThr

   return(leftSp, incCond)

#-------------------------------------------------------------------------------------
# check if the read covers the entire homopolymer stretch
#-------------------------------------------------------------------------------------
def isHPCovered(astart, aend, hpInfo):
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
def updateReadMetrics(pileupRead, base, bq, incCond, pairOrder, leftSp, umiSide, primerSide, alleleCnt, forwardCnt, reverseCnt, lowQReads, umiSideUmiEndPos, primerSideUmiEndPos, primerSidePrimerEndPos, cvg):
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
      if pairOrder == umiSide:
         # distance to the random (umi) end
         if pileupRead.alignment.is_reverse:
            distToUmiEnd = pileupRead.alignment.query_alignment_length - (pileupRead.query_position - leftSp)
         else:
            distToUmiEnd = pileupRead.query_position - leftSp
         if incCond:
            umiSideUmiEndPos[base].append(distToUmiEnd)
      if pairOrder == primerSide:
         # distance to the barcode and/or primer end on primer side read. Different cases for forward and reverse strand
         if pileupRead.alignment.is_reverse:
            distToUmiEnd = pileupRead.query_position - leftSp
            distToPrimerEnd = pileupRead.alignment.query_alignment_length - (pileupRead.query_position - leftSp)
         else:
            distToUmiEnd = pileupRead.alignment.query_alignment_length - (pileupRead.query_position - leftSp)
            distToPrimerEnd = pileupRead.query_position - leftSp
         if incCond:
            primerSideUmiEndPos[base].append(distToUmiEnd)
            primerSidePrimerEndPos[base].append(distToPrimerEnd)
   
   # coverage -- read, not fragment
   cvg += 1
   
   return(alleleCnt, forwardCnt, reverseCnt, lowQReads, umiSideUmiEndPos, primerSideUmiEndPos, primerSidePrimerEndPos, cvg)

#-------------------------------------------------------------------------------------
# Group reads by UMI and update some metrics 
#-------------------------------------------------------------------------------------  
def groupByUmi(readid, umi, base, pairOrder, usedFrag, allFrag, incCond, hpCovered, allUmiDict, umiDictHq, umiDictHqBase, umiDictAll, concordPairCnt, discordPairCnt, umiPairDict = None, umiNoDupTag = None, duplexTag = None):
   # count total number of fragments and umis
   if readid not in allUmiDict[umi]:
      allFrag += 1 # total fragments
      allUmiDict[umi].add(readid)

   # constructing umi family; this one with high quality reads only
   if incCond:
      if readid not in umiDictHq[umi]:
         readinfo = [base, pairOrder]
         umiDictHq[umi][readid] = readinfo
         # store base level information to avoid looping over read ids again
         umiDictHqBase[umi][base] += 1
         umiDictHqBase[umi]['all'] += 1
         usedFrag += 1     
      elif base == umiDictHq[umi][readid][0] or base in ['N', '*']:
         umiDictHq[umi][readid][1] = 'Paired'
         if base == umiDictHq[umi][readid][0]:
            concordPairCnt[base] += 1
      else:
         # decrement fragment and base count when R1 and R2 disagree
         usedFrag -= 1
         umiDictHqBase[umi][umiDictHq[umi][readid][0]] -= 1
         umiDictHqBase[umi]['all'] -= 1
         del umiDictHq[umi][readid]
         discordPairCnt[base] += 1
   
   # placeholder
   umiPairDict = None
      
   # in non-HP region, include all reads for consensus. In HP region, including only the reads covering the HP. 
   if hpCovered:
      #umiDictAll[umi].append(base)
      umiDictAll[umi][base] += 1
      umiDictAll[umi]['all'] += 1
   
   return(allUmiDict, umiDictHq, umiDictAll, umiDictHqBase, concordPairCnt, discordPairCnt, allFrag, usedFrag, umiPairDict)

#-------------------------------------------------------------------------------------
# Group reads by UMI and update some metrics; duplex-seq
#-------------------------------------------------------------------------------------  
def dup_groupByUmi(readid, umi, base, pairOrder, usedFrag, allFrag, incCond, hpCovered, allUmiDict, umiDictHq, umiDictHqBase, umiDictAll, concordPairCnt, discordPairCnt, umiPairDict = None, umiNoDupTag = None, duplexTag = None):
   # count total number of fragments and umis
   if readid not in allUmiDict[umi]:
      allFrag += 1 # total fragments
      allUmiDict[umi].add(readid)

   # constructing umi family; this one with high quality reads only
   if incCond:
      if readid not in umiDictHq[umi]:
         readinfo = [base, pairOrder]
         umiDictHq[umi][readid] = readinfo
         # store base level information to avoid looping over read ids again
         umiDictHqBase[umi][base] += 1
         umiDictHqBase[umi]['all'] += 1
         usedFrag += 1     
      elif base == umiDictHq[umi][readid][0] or base in ['N', '*']:
         umiDictHq[umi][readid][1] = 'Paired'
         if base == umiDictHq[umi][readid][0]:
            concordPairCnt[base] += 1
      else:
         # decrement fragment and base count when R1 and R2 disagree
         usedFrag -= 1
         umiDictHqBase[umi][umiDictHq[umi][readid][0]] -= 1
         umiDictHqBase[umi]['all'] -= 1
         del umiDictHq[umi][readid]
         discordPairCnt[base] += 1
      
      # keep track of UMIs with and without duplex tags
      umiPairDict[umiNoDupTag].add(duplexTag)
          
   # in non-HP region, include all reads for consensus. In HP region, including only the reads covering the HP. 
   if hpCovered:
      #umiDictAll[umi].append(base)
      umiDictAll[umi][base] += 1
      umiDictAll[umi]['all'] += 1
   
   return(allUmiDict, umiDictHq, umiDictAll, umiDictHqBase, concordPairCnt, discordPairCnt, allFrag, usedFrag, umiPairDict)
         
#-------------------------------------------------------------------------------------
# call pysam function to pile up an interval
#-------------------------------------------------------------------------------------
def pileup(bamName, chrom, start, end):
   samfile = pysam.AlignmentFile(bamName,"rb")
   current_pos = int(start)
   for p in samfile.pileup(region = chrom + ":" + start + ":" + end, truncate=True, max_depth = maxDnaReadDepth, stepper = "nofilter"):
      ref_pos = p.pos+1
      while True:
         if ref_pos < current_pos:
            raise Exception("pysam returned out of bounds coordinate !")
         elif ref_pos > current_pos:
            yield None
            current_pos += 1
         else:
            yield p
            current_pos += 1
            break
      
   samfile.close()
   
#-------------------------------------------------------------------------------------
# pile up reads and group by umi; some metrics are updated here
#-------------------------------------------------------------------------------------
def pileupAndGroupByUmi(bamName, bamType, chrom, pos, repType, hpInfo, minBq, minMq, hpLen, mismatchThr, primerDist, consThr, rpu, primerSide, refseq, minAltUmi, maxAltAllele, isRna, read_pileup, hqCache, infoCache, umiTag, primerTag, mqTag, tagSeparator, umiPairDict, duplexTag, getUmiFun, groupByUmiFun):
   # define variables
   cvg, usedFrag, allFrag = 0, 0, 0
   lowQReads = defaultdict(int)
   alleleCnt = defaultdict(int)
   forwardCnt = defaultdict(int)
   reverseCnt = defaultdict(int)
   concordPairCnt = defaultdict(int)
   discordPairCnt = defaultdict(int)
   umiSideUmiEndPos = defaultdict(list)
   primerSideUmiEndPos = defaultdict(list)
   primerSidePrimerEndPos = defaultdict(list)
   allUmiDict = defaultdict(set)
   umiDictHqBase = defaultdict(lambda:defaultdict(int))
   umiDictAll = defaultdict(lambda:defaultdict(int))
   umiDictHq = defaultdict(lambda:defaultdict(list))

   umiSide = 'R1' if primerSide == 'R2' else 'R2'
  
   # check pielup size, downsample if above threshold; for RNA only
   pileupStackSize = read_pileup.n
   downsamplePileup = True if isRna and pileupStackSize > downsamplePileupStackThr else False
   random.seed(pos)

   # iterate over pileup reads and group by umi
   for pileupRead in read_pileup.pileups:
      # drop reads randomly; for RNA only
      if downsamplePileup and random.randint(1, pileupStackSize) > downsamplePileupStackThr:
         continue

      # basic information that will be used in subsequent functions; NOTE: if the input BAM is consensused, use read ID as umi barcode
      readid = pileupRead.alignment.query_name
      pairOrder = 'R1' if pileupRead.alignment.is_read1 else 'R2'
      key = readid + '_' + pairOrder
      if key in infoCache:
         umi, umiNoDupTag, dupTag, cigar, astart, aend = infoCache[key]
      else:
         umi, umiNoDupTag, dupTag = getUmiFun(pileupRead, bamType, readid, umiTag, duplexTag)
         cigar, astart, aend = getBasicInfo(pileupRead)
         infoCache[key] = (umi, umiNoDupTag, dupTag, cigar, astart, aend)

      # check if read should be dropped
      if dropRead(pileupRead, pos, cigar):
         continue
      
      # drop all NN duplex tags
      if dupTag == 'NN':
         continue

      # retrive base and base-quality, and re-format the base
      base, bq = getBaseAndBq(pileupRead, refseq, chrom, pos, minBq)

      # check if the read is high quality
      hpCovered = isHPCovered(astart, aend, hpInfo)
      if key in hqCache:
         leftSp,incCondTemp = hqCache[key]
         incCond = incCondTemp and bq >= minBq and hpCovered
      else:
         leftSp, incCondTemp = hqRead(pileupRead, cigar, minMq, mismatchThr, mqTag)
         incCond = incCondTemp and bq >= minBq and hpCovered
         hqCache[key] = (leftSp, incCondTemp)

      # update read-level metrics
      alleleCnt, forwardCnt, reverseCnt, lowQReads, umiSideUmiEndPos, primerSideUmiEndPos, primerSidePrimerEndPos, cvg = updateReadMetrics(pileupRead, base, bq, incCond, pairOrder, leftSp, umiSide, primerSide, alleleCnt, forwardCnt, reverseCnt, lowQReads, umiSideUmiEndPos, primerSideUmiEndPos, primerSidePrimerEndPos, cvg)

      # group reads by umis
      allUmiDict, umiDictHq, umiDictAll, umiDictHqBase, concordPairCnt, discordPairCnt, allFrag, usedFrag, umiPairDict = groupByUmiFun(readid, umi, base, pairOrder, usedFrag, allFrag, incCond, hpCovered, allUmiDict, umiDictHq, umiDictHqBase, umiDictAll, concordPairCnt, discordPairCnt, umiPairDict, umiNoDupTag, dupTag)
      
   # output variables
   return(alleleCnt, forwardCnt, reverseCnt, lowQReads, umiSideUmiEndPos, primerSideUmiEndPos, primerSidePrimerEndPos, cvg, allUmiDict, umiDictHq, umiDictAll, umiDictHqBase, concordPairCnt, discordPairCnt, allFrag, usedFrag, hqCache, infoCache, umiPairDict)

#-------------------------------------------------------------------------------------
# gradually drop singleton UMIs, depending on rpu and input mode
#-------------------------------------------------------------------------------------
def dropSingleton(umiDictHq, minRpu, rpu = None, pos = None, ds = None, cvg = None, allUmiDict = None, isRna = None, bamType = None):
   singleUmis = set()
   pairedUmis = set()
   umiToKeep = []

   # drop singletons for duplex-seq runs
   if isRna:
      allUmi = len(allUmiDict)
      rpu = cvg / float(allUmi) if allUmi > 0 else 1.0  
      
   # rpu < 2 or consensused BAM input: no umi is dropped
   if rpu < 2.0 or bamType == 'consensus':
      umiToKeep = umiDictHq.keys()
      
   # 2 <= rpu < 3: gradually and randomly drop singleton umis 
   elif rpu >= 2.0 and rpu < 3.0:
      # set seed to be the genome position
      random.seed(pos)
      # count the numbers of paired and unpaired singleton umis; 
      pctToDrop = rpu - 2.0
      for bc in umiDictHq:
         readPairsInBc = len(umiDictHq[bc])
         if readPairsInBc == 1:
            readid = umiDictHq[bc].keys()[0]
            if umiDictHq[bc][readid][1] == 'Paired':
               pairedUmis.add(bc)
            else:
               singleUmis.add(bc)
      # total number of singleton umis
      pairedCnt = len(pairedUmis)
      singleCnt = len(singleUmis)
      oneReadUmiCnt = pairedCnt + singleCnt
      # number of singleton umis to drop
      numToDrop = int(round(pctToDrop * oneReadUmiCnt))
      # Decide which singleton umis to drop -- paired reads are kept with priority
      if numToDrop <= singleCnt:
         oneReadMtToDrop = set(random.sample(singleUmis, numToDrop))
      else:
         pairsToDrop = set(random.sample(pairedUmis, numToDrop - singleCnt))
         oneReadMtToDrop = singleUmis.union(pairsToDrop)
      # drop singleton umis
      umiToKeep = list(set(umiDictHq.keys()).difference(oneReadMtToDrop))      

   # rpu >= 3: drop UMIs with read fragments < minRpu;
   else:
      umiToKeep = [bc for bc in umiDictHq.iterkeys() if len(umiDictHq[bc]) >= minRpu]
      
   # additional downsample for RNA-seq data only
   if isRna and len(umiToKeep) > ds:
     random.seed(pos)
     umiToKeep = random.sample(umiToKeep, ds)
         
   # output variables
   return umiToKeep

#-------------------------------------------------------------------------------------
# drop UMIs with read fragments < minRpu; duplex-seq
#-------------------------------------------------------------------------------------
def dup_dropSingleton(umiDictHq, minRpu, rpu = None, pos = None, ds = None, cvg = None, allUmiDict = None, isRna = None, bamType = None):
   umiToKeep = [bc for bc in umiDictHq.iterkeys() if len(umiDictHq[bc]) >= minRpu]
   return umiToKeep

#-------------------------------------------------------------------------------------
# separate singleplex and duplex UMIs; duplex-seq 
#-------------------------------------------------------------------------------------
def dup_sepUmi(umiPairDict, umiDictHq, umiDictAll, umiToKeep):
   singleUmis = set()
   doubleUmiNoTags = set()
   for key in umiPairDict:
      tt = ':'.join(('TT',key))
      cc = ':'.join(('CC',key))
      umiPairDictVal = umiPairDict[key]

      if tt not in umiToKeep and 'TT' in umiPairDictVal:
         umiPairDictVal.discard('TT')
         if tt in umiDictHq:
            del umiDictHq[tt]
         if tt in umiDictAll:
            del umiDictAll[tt]

      if cc not in umiToKeep and 'CC' in umiPairDictVal:
         umiPairDictVal.discard('CC')
         if cc in umiDictHq:
            del umiDictHq[cc]
         if cc in umiDictAll:
            del umiDictAll[cc]
      
      nDupTag = len(umiPairDictVal)
      if nDupTag == 0:
         continue
      # single UMIs after dropping singleton UMIs
      elif nDupTag == 1:   
         fullUmi = list(umiPairDictVal)[0] + ':' + key
         singleUmis.add(fullUmi)
      # duplex UMIs after dropping singleton UMIs
      elif nDupTag == 2:   
         doubleUmiNoTags.add(key)
      else:
         # please modify the error message and triggering method to be consistent with the rest of code
         exit('UMI error: ' + key + ' has ' + str(nDupTag) + ' tags')  
   
   return(singleUmis, doubleUmiNoTags, umiDictHq, umiDictAll)
   
#-------------------------------------------------------------------------------------
# find the consensus nucleotide (including indel) in a UMI family with high quality reads only
#-------------------------------------------------------------------------------------
def consHqUmi(oneUmi, consThr):
   totalCnt = oneUmi['all']
   cons = ''
   # find the majority base(s) whose proportion >= consThr. NOTE: consThr must be > 0.5 to ensure only one cons
   for base in oneUmi:
      if base == 'all':
         continue
      pCons = 1.0 * oneUmi[base] / totalCnt if totalCnt > 0 else 0.0
      if pCons >= consThr:
         cons = base
         break
   # report the consensus base. If no consensus or lack of read support, output ''. 
   return cons

#-------------------------------------------------------------------------------------
# find the consensus nucleotide (including indel) in a UMI family with all reads 
#-------------------------------------------------------------------------------------
def consAllUmi(oneUmi, consThr):
   totalCnt = oneUmi['all']
   cons = ''
   # find the majority base(s) whose proportion >= consThr. NOTE: consThr must be > 0.5 to ensure only one cons
   for base in oneUmi:
      if base == 'all': ## just a counter
         continue
      pCons = 1.0 * oneUmi[base] / totalCnt if totalCnt > 0 else 0.0
      if pCons >= consThr:
         cons = base
         break
   # report the consensus base. If no consensus or lack of read support, output ''. 
   return cons
   
#-------------------------------------------------------------------------------------
# consensus for singleplex UMI  
#-------------------------------------------------------------------------------------
def consensus(umiDictHqBase, umiDictAll, umi, consThr, bamType):
   tmpHqUmi = umiDictHqBase[umi]
   tmpAllUmi = umiDictAll[umi]
   
   if bamType == 'raw':
      consHq = consHqUmi(tmpHqUmi, consThr)
      consAll = consAllUmi(tmpAllUmi, consThr)
      cons = consHq if consHq == consAll else ''
   else:
      if len(tmpHqUmi) == 2 and 'all' in tmpHqUmi: 
         del tmpHqUmi['all']
         cons = tmpHqUmi.keys()[0]
      else:
         cons = ''
       
   return cons

#-------------------------------------------------------------------------------------
# update single UMI metrics
#-------------------------------------------------------------------------------------
def updateUmiMetrics(umi, umiDictHqBase, cons, hqAgree, hqDisagree, umiDictAll, allAgree, allDisagree, origRef, sUmiCons, sUmiConsByBase, sUmiConsByDir, sUmiConsByDirByBase, rpuCnt, sSubTypeCnt, sUmiSnp, sStrands, tagSeparator):
   # primer ID and direction
   umiSplit = umi.split(tagSeparator)
   primerDirCode = umiSplit[1]
   primerDirection = 'F' if primerDirCode == '0' else 'R' # 0 means the primer was priming the forward strand, 1 means priming the reverse strand

   # count number of reads in concordant/discordant with consensus for UMI efficiency metrics
   for base in umiDictHqBase[umi]:
      if base == 'all': ## just a counter
         continue
      if base == cons:
         hqAgree[base] += umiDictHqBase[umi][base]
      else:
         hqDisagree[base] += umiDictHqBase[umi][base]

   for base in umiDictAll[umi]:
      if base == 'all': ## just a counter
         continue
      if base == cons:
         allAgree[base] += umiDictAll[umi][base]
      else:
         allDisagree[base] += umiDictAll[umi][base]

   if cons != '':
      sUmiCons += 1
      sUmiConsByBase[cons] += 1
      # UMI counts from + and - strands 
      sUmiConsByDir[primerDirection] += 1
      sUmiConsByDirByBase[cons][primerDirection] += 1
      # read pairs in the umi
      rpuCnt[cons].append(umiDictAll[umi]['all'])
      
      # base substitutions (snp only)
      # Note: sUmiSnp and strands are usually NOT equal to sUmiCons and sUmiConsByDir. The former include only base substitutions UMIs, and the latter include indel UMIs. 
      if len(cons) == 1:
         basePair = origRef + '/' + cons if primerDirCode == '0' else origRef.translate(_base_complement_) + '/' + cons.translate(_base_complement_)
         sSubTypeCnt[basePair] += 1
         sUmiSnp += 1
         sStrands[primerDirection] += 1

   return(hqAgree, hqDisagree, allAgree, allDisagree, sUmiCons, sUmiConsByBase, sUmiConsByDir, sUmiConsByDirByBase, rpuCnt, sSubTypeCnt, sUmiSnp, sStrands)

#-------------------------------------------------------------------------------------
# consensus for UMI and update UMI metrics; duplex-seq   
#-------------------------------------------------------------------------------------
def dup_consAndUpdateUmiMetrics(umiNoDupTag, umiDictHqBase, umiDictAll, consThr, hqAgree, hqDisagree, allAgree, allDisagree, origRef, dUmiCons, dUmiConsByBase, dUmiConsByDir, dUmiConsByDirByBase, rpuCnt, dSubTypeCnt, dUmiSnp, dStrands, tagSeparator, discordDupPairs):
   consHqDict = defaultdict(str)
   consAllDict = defaultdict(str)
   
   # primer ID and direction
   umiSplit = umiNoDupTag.split(tagSeparator)
   primerDirCode = umiSplit[1]
   primerDirection = 'F' if primerDirCode == '0' else 'R' # 0 means the primer was priming the forward strand, 1 means priming the reverse strand
   
   for dupTag in ['CC', 'TT']:
      # full UMI string
      umi = dupTag + ':' + umiNoDupTag
      
      tmpHqUmi = umiDictHqBase[umi]
      tmpAllUmi = umiDictAll[umi]
      
      # consensus by strand, both Hq and All
      hqCons = consHqUmi(tmpHqUmi, consThr)
      allCons = consAllUmi(tmpAllUmi, consThr)
      consHqDict[dupTag] = hqCons
      consAllDict[dupTag] = allCons
      strandCons = hqCons if hqCons == allCons else ''
      
      # count number of reads in concordant/discordant with consensus for UMI efficiency metrics
      for base in umiDictHqBase[umi]:
         if base == 'all': ## just a counter
            continue
         if base == strandCons:
            hqAgree[base] += umiDictHqBase[umi][base]
         else:
            hqDisagree[base] += umiDictHqBase[umi][base]

      for base in umiDictAll[umi]:
         if base == 'all': ## just a counter
            continue
         if base == strandCons:
            allAgree[base] += umiDictAll[umi][base]
         else:
            allDisagree[base] += umiDictAll[umi][base]
   
   consHq_CC = consHqDict['CC']
   consHq_TT = consHqDict['TT']
   consAll_CC = consAllDict['CC']
   consAll_TT = consAllDict['TT']

   # duplex consensus if all 4 (CC/TT, Hq/All) cons are identical and not equal to ''
   if consHq_CC == consHq_TT == consAll_CC == consAll_TT and consHq_CC != '':
      cons = consHq_CC
   else:
      cons = ''
   
   # update duplex-UMI metrics
   if cons != '':
      dUmiCons += 1
      dUmiConsByBase[cons] += 1
      # UMI counts from + and - strands 
      dUmiConsByDir[primerDirection] += 1
      dUmiConsByDirByBase[cons][primerDirection] += 1
      # read pairs in the duplex UMI
      rpuCnt[cons].append(umiDictAll['CC:' + umiNoDupTag]['all'])
      rpuCnt[cons].append(umiDictAll['TT:' + umiNoDupTag]['all'])
      
   # base substitutions (snp only) of duplex UMI
   if len(cons) == 1:
      basePair = origRef + '/' + cons if primerDirCode == '0' else origRef.translate(_base_complement_) + '/' + cons.translate(_base_complement_)
      dSubTypeCnt[basePair] += 1
      dUmiSnp += 1
      dStrands[primerDirection] += 1
         
   # discordant duplex UMIs
   if (consHq_CC != consHq_TT and consHq_CC != '' and consHq_TT != '') or (consAll_CC != consAll_TT and consAll_CC != '' and consAll_TT != ''):
      discordDupPairs[consHq_CC] += 1
      discordDupPairs[consHq_TT] += 1

   return (cons, hqAgree, hqDisagree, allAgree, allDisagree, dUmiCons, dUmiConsByBase, dUmiConsByDir, dUmiConsByDirByBase, rpuCnt, dSubTypeCnt, dUmiSnp, dStrands, discordDupPairs)

#-------------------------------------------------------------------------------------
# save the background error profile; normal DNA-seq 
#-------------------------------------------------------------------------------------
def outBkg(chrom, pos, origRef, sSubTypeCnt, sStrands, sUmiSnp, dSubTypeCnt = None, dStrands = None, dUmiSnp = None):
   # single UMI errors 
   sinBkgErrList = [chrom, pos, origRef, str(sSubTypeCnt['A/G']), str(sSubTypeCnt['G/A']), str(sSubTypeCnt['C/T']), str(sSubTypeCnt['T/C']), str(sSubTypeCnt['A/C']), str(sSubTypeCnt['C/A']), str(sSubTypeCnt['A/T']), str(sSubTypeCnt['T/A']), str(sSubTypeCnt['C/G']), str(sSubTypeCnt['G/C']), str(sSubTypeCnt['G/T']), str(sSubTypeCnt['T/G']), str(sStrands['F']), str(sStrands['R']), str(sUmiSnp), 'single']
   outLineBkg = '\t'.join(sinBkgErrList) + '\n'

   # output variables
   return outLineBkg 

#-------------------------------------------------------------------------------------
# save the background error profile; duplex-seq
#-------------------------------------------------------------------------------------
def dup_outBkg(chrom, pos, origRef, sSubTypeCnt, sStrands, sUmiSnp, dSubTypeCnt = None, dStrands = None, dUmiSnp = None):
   # single UMI errors 
   sinBkgErrList = [chrom, pos, origRef, str(sSubTypeCnt['A/G']), str(sSubTypeCnt['G/A']), str(sSubTypeCnt['C/T']), str(sSubTypeCnt['T/C']), str(sSubTypeCnt['A/C']), str(sSubTypeCnt['C/A']), str(sSubTypeCnt['A/T']), str(sSubTypeCnt['T/A']), str(sSubTypeCnt['C/G']), str(sSubTypeCnt['G/C']), str(sSubTypeCnt['G/T']), str(sSubTypeCnt['T/G']), str(sStrands['F']), str(sStrands['R']), str(sUmiSnp), 'single']
   outLineBkg = '\t'.join(sinBkgErrList) + '\n'

   # duplex UMI errors
   dupBkgErrList = [chrom, pos, origRef, str(dSubTypeCnt['A/G']), str(dSubTypeCnt['G/A']), str(dSubTypeCnt['C/T']), str(dSubTypeCnt['T/C']), str(dSubTypeCnt['A/C']), str(dSubTypeCnt['C/A']), str(dSubTypeCnt['A/T']), str(dSubTypeCnt['T/A']), str(dSubTypeCnt['C/G']), str(dSubTypeCnt['G/C']), str(dSubTypeCnt['G/T']), str(dSubTypeCnt['T/G']), str(dStrands['F']), str(dStrands['R']), str(dUmiSnp), 'duplex']
   
   # concatenate singleplex and duplex output lines
   outLineBkg += '\t'.join(dupBkgErrList) + '\n'

   # output variables
   return outLineBkg 
   
#-------------------------------------------------------------------------------------
# reset variant type, reference base, variant base 
#-------------------------------------------------------------------------------------
def setRefAltType(origRef, origAlt):
   vType = '.'
   ref = origRef
   alt = origAlt

   if len(origAlt) == 1:
      vType = 'SNP'
   elif origAlt == 'DEL':
      vType = 'SDEL'
   else:
      vals = origAlt.split('|')
      if vals[0] in ['DEL', 'INS']:
         vType = 'INDEL'
         ref = vals[1]
         alt = vals[2]

   return(ref, alt, vType)

#-------------------------------------------------------------------------------------
# compute UMI efficiency metrics
#-------------------------------------------------------------------------------------
def umiEfficiency(hqAgree, hqDisagree, allAgree, allDisagree, origRef, origAlt, rpuCnt, alleleCnt, sUmiConsByBase, cvg, sUmiCons, tmpVaf, tmpVmf):
   hqRcAgree = hqAgree[origAlt] 
   hqRcTotal = hqRcAgree + hqDisagree[origAlt] 
   hqUmiEff = 1.0 * hqRcAgree / hqRcTotal if hqRcTotal > 0 else 0.0
   
   allRcAgree = allAgree[origAlt] 
   allRcTotal = allRcAgree + allDisagree[origAlt] 
   allUmiEff = 1.0 * allRcAgree / allRcTotal if allRcTotal > 0 else 0.0

   if sUmiConsByBase[origRef] >= 3 and sUmiConsByBase[origAlt] >= 3:
      refRppUmiN = sUmiConsByBase[origRef]
      refRppUmiMean = numpy.mean(rpuCnt[origRef])
      refRppUmiSd = numpy.std(rpuCnt[origRef])
      altRppUmiN = sUmiConsByBase[origAlt]
      altRppUmiMean = numpy.mean(rpuCnt[origAlt])
      altRppUmiSd = numpy.std(rpuCnt[origAlt])
      sp = ( ((refRppUmiN - 1) * refRppUmiSd ** 2 + (altRppUmiN - 1) * altRppUmiSd ** 2) / (refRppUmiN + altRppUmiN - 2) ) ** 0.5
      RppEffSize = (refRppUmiMean - altRppUmiMean) / (sp * (1.0 / refRppUmiN + 1.0 / altRppUmiN) ** 0.5) if sp > 0 else 1000.0
   else:
      refRppUmiMean = -1.0
      altRppUmiMean = -1.0
      RppEffSize = -1.0

   vafToVmfRatio = 1.0 * tmpVaf / tmpVmf if tmpVmf > 0 else -1.0

   return (vafToVmfRatio, hqUmiEff, allUmiEff, refRppUmiMean, altRppUmiMean, RppEffSize)

#-------------------------------------------------------------------------------------
# detailed output file
#-------------------------------------------------------------------------------------
def outLong(outLine, chrom, pos, ref, alt, vType, origRef, origAlt, sUmiCons, sUmiConsByDir, sUmiConsByBase, sUmiConsByDirByBase, alleleCnt, primerBiasOR, hqUmiEff, allUmiEff, refRppUmiMean, altRppUmiMean, RppEffSize, repTypeSet, bqAlt, hpInfo, srInfo, repInfo, cvg, allFrag, allUmiDict, usedFrag, fltrs, dUmiCons = None, dUmiConsByBase = None, discordDupPairs = None, dUmiConsByDirByBase = None):
   # total number of UMIs, fragments, reads, including those dropped from analysis
   allUmi = len(allUmiDict)
   # final FILTER to output
   fltrFinal = 'PASS' if len(fltrs) == 0 else ';'.join(list(fltrs)) 
   # read-based variant allele fraction (VAF); based on all reads
   fracAlt = str(round((100.0 * alleleCnt[origAlt] / cvg), 3)) if cvg > 0 else '.'  
   # UMI-based variant allele fraction (VMF)
   sVmf = str(round((100.0 * sUmiConsByBase[origAlt] / sUmiCons), 5)) if sUmiCons > 0 else '.'  
   # UMI count for A,C,G,T
   sUmis = [str(sUmiConsByBase['A']), str(sUmiConsByBase['T']), str(sUmiConsByBase['G']), str(sUmiConsByBase['C'])]
   # proportion of <Q20 reads
   pLowQ = str(round(bqAlt, 2)) if bqAlt >= 0 else 'NA'
   # type of repetitive region
   repTypeFinal = ';'.join(list(repTypeSet)) if len(repTypeSet) >= 1 else 'NA'
   # UMI counts by primer direction
   refForPrimer = sUmiConsByDirByBase[origRef]['F']
   refRevPrimer = sUmiConsByDirByBase[origRef]['R']
   altForPrimer = sUmiConsByDirByBase[origAlt]['F']
   altRevPrimer = sUmiConsByDirByBase[origAlt]['R']
   # UMI-based VMF for each strand
   vmfForward = str(round((100.0 * sUmiConsByDirByBase[origAlt]['F'] / sUmiConsByDir['F']), 3)) if sUmiConsByDir['F'] > 0 else '.'
   vmfReverse = str(round((100.0 * sUmiConsByDirByBase[origAlt]['R'] / sUmiConsByDir['R']), 3)) if sUmiConsByDir['R'] > 0 else '.'
   # round hqUmiEff and allUmiEff in the end
   hqUmiEff = round(hqUmiEff, 3)
   allUmiEff = round(allUmiEff, 3)

   outList = [chrom, pos, ref, alt, vType, str(sUmiCons), str(sUmiConsByDir['F']), str(sUmiConsByDir['R']), str(sUmiConsByBase[origAlt]), str(sUmiConsByDirByBase[origAlt]['F']), str(sUmiConsByDirByBase[origAlt]['R']), sVmf, vmfForward, vmfReverse, str(alleleCnt[origAlt]), fracAlt, str(refForPrimer), str(refRevPrimer), primerBiasOR, pLowQ, str(hqUmiEff), str(allUmiEff), str(refRppUmiMean), str(altRppUmiMean), str(RppEffSize), repTypeFinal, hpInfo, srInfo, repInfo, str(cvg), str(allFrag), str(allUmi), str(usedFrag)] + sUmis + [fltrFinal]
   
   outLineAllele = '\t'.join(outList) + '\n'
   outLine += outLineAllele

   return outLine
   
#-------------------------------------------------------------------------------------
# detailed output file; duplex-seq
#-------------------------------------------------------------------------------------
def dup_outLong(outLine, chrom, pos, ref, alt, vType, origRef, origAlt, sUmiCons, sUmiConsByDir, sUmiConsByBase, sUmiConsByDirByBase, alleleCnt, primerBiasOR, hqUmiEff, allUmiEff, refRppUmiMean, altRppUmiMean, RppEffSize, repTypeSet, bqAlt, hpInfo, srInfo, repInfo, cvg, allFrag, allUmiDict, usedFrag, fltrs, dUmiCons = None, dUmiConsByBase = None, discordDupPairs = None, dUmiConsByDirByBase = None):
   # total number of UMIs, fragments, reads, including those dropped from analysis
   allUmi = len(allUmiDict)
   # final FILTER to output
   fltrFinal = 'PASS' if len(fltrs) == 0 else ';'.join(list(fltrs)) 
   # read-based variant allele fraction (VAF); based on all reads
   fracAlt = str(round((100.0 * alleleCnt[origAlt] / cvg), 3)) if cvg > 0 else '.'  
   # UMI-based variant allele fraction (VMF)
   sVmf = str(round((100.0 * sUmiConsByBase[origAlt] / sUmiCons), 3)) if sUmiCons > 0 else '.'  
   dVmf = str(round((100.0 * dUmiConsByBase[origAlt] / dUmiCons), 3)) if dUmiCons > 0 else '.'  
   # UMI count for A,C,G,T
   sUmis = [str(sUmiConsByBase['A']), str(sUmiConsByBase['T']), str(sUmiConsByBase['G']), str(sUmiConsByBase['C'])]
   dUmis = [str(dUmiConsByBase['A']), str(dUmiConsByBase['T']), str(dUmiConsByBase['G']), str(dUmiConsByBase['C'])]
   # discordant duplex UMI pairs
   vdm = discordDupPairs[origAlt]
   # proportion of <Q20 reads
   pLowQ = str(round(bqAlt, 3)) if bqAlt >= 0 else 'NA'
   # type of repetitive region
   repTypeFinal = ';'.join(list(repTypeSet)) if len(repTypeSet) >= 1 else 'NA'
   # single UMI counts by primer direction
   sForUmt = sUmiConsByDirByBase[origRef]['F']
   sRevUmt = sUmiConsByDirByBase[origRef]['R']
   sForVmt = sUmiConsByDirByBase[origAlt]['F']
   sRevVmt = sUmiConsByDirByBase[origAlt]['R']
   # duplex UMI counts by primer direction
   dForUmt = dUmiConsByDirByBase[origRef]['F']
   dRevUmt = dUmiConsByDirByBase[origRef]['R']
   dForVmt = dUmiConsByDirByBase[origAlt]['F']
   dRevVmt = dUmiConsByDirByBase[origAlt]['R']
   
   # round hqUmiEff and allUmiEff in the end
   hqUmiEff = round(hqUmiEff, 3)
   allUmiEff = round(allUmiEff, 3)

   outList = [chrom, pos, ref, alt, vType, str(sUmiCons), str(sUmiConsByBase[origAlt]), sVmf, str(dUmiCons), str(dUmiConsByBase[origAlt]), dVmf, str(cvg), str(alleleCnt[origAlt]), fracAlt, str(sForUmt), str(sForVmt), str(sRevUmt), str(sRevVmt), str(dForUmt), str(dForVmt), str(dRevUmt), str(dRevVmt), primerBiasOR, pLowQ, str(hqUmiEff), str(allUmiEff), str(refRppUmiMean), str(altRppUmiMean), str(RppEffSize), repTypeFinal, hpInfo, srInfo, repInfo, str(allFrag), str(allUmi), str(usedFrag)] + sUmis + dUmis + [fltrFinal]
   
   outLineAllele = '\t'.join(outList) + '\n'
   outLine += outLineAllele

   return outLine
   
#-------------------------------------------------------------------------------------
# function to call variants
#-------------------------------------------------------------------------------------
def vc(bamName, chrom, pos, repType, hpInfo, srInfo, repInfo, minBq, minMq, hpLen, mismatchThr, primerDist, consThr, rpu, primerSide, refseq, minAltUmi, maxAltAllele, isRna, ds, bamType, read_pileup, hqCache, infoCache, chromLength, umiTag, primerTag, mqTag, tagSeparator, nCols, defVarFun, getUmiFun, groupByUmiFun, dropSingletonFun, outBkgFun, outLongFun, isDuplex, duplexTag, minRpu):
   # initiate variables
   sUmiCons, sUmiSnp, sUmiConsByBase, sUmiConsByDir, sUmiConsByDirByBase, sStrands, sSubTypeCnt, hqAgree, hqDisagree, allAgree, allDisagree, rpuCnt, umiPairDict, sUmiConsByBase, outLineLong, dUmiCons, dUmiSnp, dUmiConsByBase, dUmiConsByDir, dUmiConsByDirByBase, dStrands, dSubTypeCnt, discordDupPairs = defVarFun()
   
   # find the reference base
   origRef = getRef(refseq, chrom, pos)
   
   # pile up and group reads by UMIs
   alleleCnt, forwardCnt, reverseCnt, lowQReads, umiSideUmiEndPos, primerSideUmiEndPos, primerSidePrimerEndPos, cvg, allUmiDict, umiDictHq, umiDictAll, umiDictHqBase, concordPairCnt, discordPairCnt, allFrag, usedFrag, hqCache, infoCache, umiPairDict = pileupAndGroupByUmi(bamName, bamType, chrom, pos, repType, hpInfo, minBq, minMq, hpLen, mismatchThr, primerDist, consThr, rpu, primerSide, refseq, minAltUmi, maxAltAllele, isRna, read_pileup, hqCache, infoCache, umiTag, primerTag, mqTag, tagSeparator, umiPairDict, duplexTag, getUmiFun, groupByUmiFun)
   
   # gradually drop singleton UMIs; keep all UMIs for consensused BAM; drop singletons for duplex-seq runs
   umiToKeep = dropSingletonFun(umiDictHq, minRpu, rpu, pos, ds, cvg, allUmiDict, isRna, bamType)  
 
   # output zeros or blank if the remaining UMIs is lower than the minimum threshold 
   if len(umiToKeep) <= minTotalUmi:
      outLineLong = '\t'.join([chrom, pos, origRef] + ['0'] * (nCols - 4) + ['LM']) + '\n'
      outLineBkg = ''
      return(outLineLong, outLineBkg, hqCache, infoCache)
   
   # run the following commands if UMI depth exceeds minTotalUmi 
   ## duplex-seq run
   if isDuplex:
      # separate single and duplex UMIs
      umiToKeep = set(umiToKeep)
      singleUmis, doubleUmiNoTags, umiDictHq, umiDictAll = dup_sepUmi(umiPairDict, umiDictHq, umiDictAll, umiToKeep)
      # process single UMIs
      for umi in singleUmis:
         # generate consensus
         cons = consensus(umiDictHqBase, umiDictAll, umi, consThr, bamType)
         
         # update umi-level metrics
         hqAgree, hqDisagree, allAgree, allDisagree, sUmiCons, sUmiConsByBase, sUmiConsByDir, sUmiConsByDirByBase, rpuCnt, sSubTypeCnt, sUmiSnp, sStrands = updateUmiMetrics(umi, umiDictHqBase, cons, hqAgree, hqDisagree, umiDictAll, allAgree, allDisagree, origRef, sUmiCons, sUmiConsByBase, sUmiConsByDir, sUmiConsByDirByBase, rpuCnt, sSubTypeCnt, sUmiSnp, sStrands, tagSeparator)    
      # process duplex UMIs
      for umiNoDupTag in doubleUmiNoTags:                  
         cons, hqAgree, hqDisagree, allAgree, allDisagree, dUmiCons, dUmiConsByBase, dUmiConsByDir, dUmiConsByDirByBase, rpuCnt, dSubTypeCnt, dUmiSnp, dStrands, discordDupPairs = dup_consAndUpdateUmiMetrics(umiNoDupTag, umiDictHqBase, umiDictAll, consThr, hqAgree, hqDisagree, allAgree, allDisagree, origRef, dUmiCons, dUmiConsByBase, dUmiConsByDir, dUmiConsByDirByBase, rpuCnt, dSubTypeCnt, dUmiSnp, dStrands, tagSeparator, discordDupPairs)  
   
   ## normal DNA-seq run
   else:
      for umi in umiToKeep:
         # generate consensus
         cons = consensus(umiDictHqBase, umiDictAll, umi, consThr, bamType)
         # update umi-level metrics
         hqAgree, hqDisagree, allAgree, allDisagree, sUmiCons, sUmiConsByBase, sUmiConsByDir, sUmiConsByDirByBase, rpuCnt, sSubTypeCnt, sUmiSnp, sStrands = updateUmiMetrics(umi, umiDictHqBase, cons, hqAgree, hqDisagree, umiDictAll, allAgree, allDisagree, origRef, sUmiCons, sUmiConsByBase, sUmiConsByDir, sUmiConsByDirByBase, rpuCnt, sSubTypeCnt, sUmiSnp, sStrands, tagSeparator)

   # output the background error profile
   outLineBkg = outBkgFun(chrom, pos, origRef, sSubTypeCnt, sStrands, sUmiSnp, dSubTypeCnt, dStrands, dUmiSnp)

   alleleList = sorted(sUmiConsByBase.items(), key = operator.itemgetter(1), reverse = True)
   firstAlt = True
   altCnt = 0
   repTypeSet0 = set() if repType == 'NA' else set(repType.strip().split(';'))

   # start multi-allelic loop
   for alleleInd in range(len(alleleList)):
      origAlt = alleleList[alleleInd][0]
      maxVMT = alleleList[alleleInd][1]

      # if a non-reference allele has >= 3 UMIs, consider it as a candidate variant 
      if origAlt == origRef:
         continue
      if maxVMT < minAltUmi and not firstAlt:
         break

      # reset variant type, reference base, variant base 
      ref, alt, vType = setRefAltType(origRef, origAlt) 

      # initiate values for filters and output
      primerBiasOR, bqAlt, hqUmiEff, allUmiEff, refRppUmiMean, altRppUmiMean, RppEffSize = 'NA', -1.0, 0.0, 0.0, -1.0, -1.0, -1.0
      fltrs = set()
      repTypeSet = repTypeSet0
    
      if vType in ['SNP', 'INDEL'] and sUmiCons > 0:
         # compute read- and umi-level allele frequency and proportion of low base quality reads supporting the variant
         tmpVaf = 100.0 * alleleCnt[origAlt] / cvg if cvg > 0 else 0.0 
         tmpVmf = 100.0 * sUmiConsByBase[origAlt] / sUmiCons

         ### common filters for both original and consensused BAM
         # low coverage filter
         fltrs = filters.lm(fltrs, sUmiCons) 
         # strand bias and discordant pairs filter
         fltrs = filters.dpSb(fltrs, origAlt, concordPairCnt, discordPairCnt, reverseCnt, forwardCnt, origRef, tmpVaf)

         # Initial HP and LowC region filters
         repTypeSet, hpInfo = filters.isHPorLowComp(chrom, pos, hpLen, ref, alt, refseq, repTypeSet, hpInfo, chromLength)

         # SNP-only common filters
         if vType == 'SNP':
            # random (umi) end position filters
            fltrs = filters.rbcp(fltrs, endBase, umiSideUmiEndPos, origRef, origAlt, tmpVaf, isRna)
            fltrs = filters.rpcp(fltrs, endBase, primerSideUmiEndPos, origRef, origAlt, tmpVaf, isRna)
            # proportion of low base quality reads
            if origAlt in alleleCnt and origAlt in lowQReads and alleleCnt[origAlt] > 0:
               bqAlt = 1.0 * lowQReads[origAlt] / alleleCnt[origAlt]
            
         ### original BAM only filters that use umi efficiency metrics
         if bamType == 'raw':
            # umi efficiency metrics; only for original BAM
            vafToVmfRatio, hqUmiEff, allUmiEff, refRppUmiMean, altRppUmiMean, RppEffSize = umiEfficiency(hqAgree, hqDisagree, allAgree, allDisagree, origRef, origAlt, rpuCnt, alleleCnt, sUmiConsByBase, cvg, sUmiCons, tmpVaf, tmpVmf)               
            if not isRna:
               # update HP for indel
               fltrs = filters.hp4indel(fltrs, repTypeSet, vType, rpu, hpInfo, vafToVmfRatio, tmpVmf, hqUmiEff, RppEffSize, altRppUmiMean)             
               # update other repetitive region filters for SNP and indel, including HP for SNP
               fltrs = filters.rep4others(fltrs, repTypeSet, vType, rpu, vafToVmfRatio, hqUmiEff, RppEffSize)
            # primer bias filter
            fltrs, primerBiasOR = filters.pb(fltrs, origAlt, sUmiConsByDir, sUmiConsByDirByBase)
            if vType == 'SNP':
               # low base quality filter - not for duplex-seq
               if not isDuplex:
                  fltrs = filters.lowq(fltrs, lowQReads, alleleCnt, origAlt, vafToVmfRatio, bqAlt, isRna)
               # fixed end (gene specific primers) position filter
               fltrs = filters.primercp(fltrs, primerDist, primerSidePrimerEndPos, origRef, origAlt, tmpVmf, hqUmiEff, vafToVmfRatio, RppEffSize, rpu)
           
         ### consensus BAM only filters; use more strict repetitive and primerOR filters; Discard LowQ filter because the base qualities have different meanings
         else:
            fltrs = filters.strict(fltrs, repTypeSet, vType, tmpVmf, primerDist, primerSidePrimerEndPos, origRef, origAlt)

      firstAlt = False
      # output metrics for each non-reference allele with >= 3 UMIs; If none, output the one with most UMIs
      outLineLong = outLongFun(outLineLong, chrom, pos, ref, alt, vType, origRef, origAlt, sUmiCons, sUmiConsByDir, sUmiConsByBase, sUmiConsByDirByBase, alleleCnt, primerBiasOR, hqUmiEff, allUmiEff, refRppUmiMean, altRppUmiMean, RppEffSize, repTypeSet, bqAlt, hpInfo, srInfo, repInfo, cvg, allFrag, allUmiDict, usedFrag, fltrs, dUmiCons, dUmiConsByBase, discordDupPairs, dUmiConsByDirByBase)
      
      altCnt += 1
      if altCnt >= maxAltAllele:
         break
   
   return (outLineLong, outLineBkg, hqCache, infoCache)

#------------------------------------------------------------------------------------------------
# wrapper function for "vc()" - because Python multiprocessing module does not pass stack trace
#------------------------------------------------------------------------------------------------
def vc_wrapper(general_args, interval):
   timeStart = datetime.datetime.now()
   
   try:
      output = []
      hqCache = {}
      infoCache = {}
      bamName, minBq, minMq, hpLen, mismatchThr, primerDist, consThr, rpu, primerSide, refg, minAltUmi, maxAltAllele, isRna, ds, bamType, umiTag, primerTag, mqTag, tagSeparator, isDuplex, duplexTag, minRpu = general_args

      chrom = interval[0][0]
      intervalStartPos = interval[0][1]
      intervalEndPos = interval[-1][1]
      
      # constants and functions to use depending on normal DNA-seq or duplex-seq
      if isDuplex:
         nCols = nColsDup
         defVarFun = dup_defVar
         outBkgFun = dup_outBkg
         outLongFun = dup_outLong
         getUmiFun = dup_getUmi
         dropSingletonFun = dup_dropSingleton
         groupByUmiFun = dup_groupByUmi
         
         if isRna or bamType != 'raw':
            exit("duplex-seq must be on DNA and with raw BAM input for now.")
            
      else:
         nCols = nColsSin
         defVarFun = defVar
         outBkgFun = outBkg
         outLongFun = outLong  
         getUmiFun = getUmi
         dropSingletonFun = dropSingleton
         groupByUmiFun = groupByUmi
         
         if duplexTag != None:
            exit("normal DNA-seq runs don't have duplex tag.")
      
      refseq = pysam.FastaFile(refg)
      chromLengths = {}
      for idx in range(len(refseq.lengths)):
         chromLengths[refseq.references[idx]] = refseq.lengths[idx]
         
      i = 0
      for read_pileup in pileup(bamName, chrom, intervalStartPos, intervalEndPos):
         site = interval[i]
         i += 1
         chrom, pos, repType, hpInfo, srInfo, repInfo = site
         
         if read_pileup is None: # site is not covered at all, pysam simply skips such sites
            origRef = getRef(refseq, chrom, pos)
            outLineLong = '\t'.join([chrom, pos, origRef] + ["0"] * (nCols - 4) + ["LM"]) + "\n"
            out = [outLineLong,""]
            output.append(out)
            continue

         temp = [bamName, chrom, pos, repType, hpInfo, srInfo, repInfo, minBq, minMq, hpLen, mismatchThr, primerDist, consThr, rpu, primerSide, refseq, minAltUmi, maxAltAllele, isRna, ds, bamType, read_pileup, hqCache, infoCache, chromLengths[chrom], umiTag, primerTag, mqTag, tagSeparator, nCols, defVarFun, getUmiFun, groupByUmiFun, dropSingletonFun, outBkgFun, outLongFun, isDuplex, duplexTag, minRpu]
                
         outLineLong, outLineBkg, hqCache, infoCache = vc(*temp)
         out = [outLineLong, outLineBkg]
         output.append(out)
         
   except Exception as e:
      out = ("Exception thrown!\n" + traceback.format_exc(), "no_bg")
      output.append(out)
      logger.info("Exception thrown in vc() function at genome location : {pos} in interval : {chrom}:{it1}-{it2}".format(pos = pos, chrom = chrom, it1 = intervalStartPos, it2 = intervalEndPos))
      logger.info(out[0])
      raise Exception(e)   
     
   refseq.close()
   timeEnd = datetime.datetime.now()   
   logger.info(str(timeEnd - timeStart) + "\t" + "{chrom}:{it1}-{it2}".format(chrom = chrom,it1 = intervalStartPos,it2 = intervalEndPos))
   return output
