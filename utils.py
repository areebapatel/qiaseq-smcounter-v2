import os
import sys
import subprocess
from collections import defaultdict
import pysam

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
# find homopolymer sequences
#-------------------------------------------------------------------------------------
def findhp(bedName, outName, minLength,refg,seqType='dna'):
   # how much to extend the roi to search for homopolymers
   extensionLen = 0 if seqType == 'rna' else 100
   
   # loop over roi BED
   outfile = open(outName, 'w')
   for line in open(bedName, 'r'):
      if line.startswith('track name='):
         continue
      lineList = line.strip().split('\t')
      chrom = lineList[0]
      start = int(lineList[1])
      end = int(lineList[2])

      # get reference base
      refseq = pysam.FastaFile(refg)
      
      start_coord = start - 1 - extensionLen
      if start_coord < 0:
         start_coord = start
      origRef = refseq.fetch(reference=chrom, start=start_coord, end=end + extensionLen)
      origRef = origRef.upper()

      hpL = 0
      for i in range(len(origRef))[1:]:
         if origRef[i] == origRef[i-1]:
            continue
         else:
            hpLen = i - hpL 
            realL = hpL - 1 + start - extensionLen
            realR = i - 1  + start - extensionLen
            if hpLen >= minLength and realL <= end and realR >= start:
               outline = '\t'.join([chrom, str(max(realL, start)), str(min(realR, end)), 'HP', str(hpLen), str(realL), str(realR), origRef[hpL]]) + '\n'
               outfile.write(outline)
            hpL = i 


   outfile.close()

#----------------------------------------------------------------------------------------------
# convert input VCF to BED file; assume the VCF contains only primitives, and follows the conventional VCF format
#------------------------------------------------------------------------------------------------
def vcf2bed(inputVCF):
   outBed = inputVCF + '.bed'
   outf = open(outBed, 'w')
   for line in open(inputVCF, 'r'):
      if line[0] == '#':
	continue
      lineList = line.strip().split('\t')
      # only the chrom and position matters
      chrom, pos = lineList[:2]

      try:
         ipos = int(pos)
      except ValueError:
         exit('Failed to convert position to integer')
      
      # create BED output line
      bed_chrom = chrom if chrom.startswith('chr') else 'chr' + chrom
      bed_start = str(max(0, ipos - 1))
      bed_end = pos
        
      # save to output
      outline = '\t'.join([bed_chrom, bed_start, bed_end]) + '\n'
      outf.write(outline)
	  
   outf.close()
   return(outBed)   
	     
#----------------------------------------------------------------------------------------------
# get homopolymer region information
#------------------------------------------------------------------------------------------------
def getHpInfo(bedTarget, refGenome, isRna, hpLen):
   # intersect repeats and target regions
   if isRna: 
      seqType = 'rna'
      findHpLen = hpLen
   else:   
      seqType = 'dna'
      findHpLen = 6

   findhp(bedTarget, 'hp.roi.bed', findHpLen, refGenome, seqType)
   
   # gather homopolymer region info
   hpRegion = defaultdict(list)
   with open('hp.roi.bed','r') as IN:
      for line in IN:   
         chrom, regionStart, regionEnd, repType, totalLen, realL, realR, repBase = line.strip().split()
         hpRegion[chrom].append([regionStart, regionEnd, repType, totalLen, realL, realR])
   
   # output variables
   return(hpRegion)
   
#----------------------------------------------------------------------------------------------
# get tandem region information
#------------------------------------------------------------------------------------------------
def getTrInfo(bedTarget, repBed, isRna, hpLen):
   # intersect repeats and target regions
   subprocess.check_call('bedtools intersect -a ' + repBed + ' -b ' + bedTarget + ' | bedtools sort -i > rep.roi.bed', shell=True)
   
   # gather tandem repeat region info
   repRegion = defaultdict(list)
   with open('rep.roi.bed','r') as IN:
      for line in IN:
         chrom, regionStart, regionEnd, repInfo = line.strip().split()[:4]
         unitLen, repLen = repInfo.split("|")[1:3]
         try:
            unitLen_num = float(unitLen)
         except ValueError:
            continue
         try:
            repLen_num = float(repLen)
         except ValueError:
            continue

         if isRna:
            totalLen = int(regionEnd) - int(regionStart)
            if totalLen < hpLen:
               continue
            repLen = str(totalLen / unitLen_num)
            totalLen = str(totalLen)
         else:
            totalLen = str(unitLen_num * repLen_num)
            
         repBase = repInfo[-1]
         repType = 'RepT'
         repRegion[chrom].append([regionStart, regionEnd, repType, totalLen, unitLen, repLen])
   
   # output variables
   return(repRegion)
   
#----------------------------------------------------------------------------------------------
# get other repeats region (simple repeats, low complexity, micro-satelites) information
#------------------------------------------------------------------------------------------------
def getOtherRepInfo(bedTarget, srBed, isRna, hpLen):
   # intersect repeats and target regions
   subprocess.check_call('bedtools intersect -a ' + srBed +  ' -b ' + bedTarget + ' | bedtools sort -i > sr.roi.bed', shell=True)   
   
   # gather other repeat region info
   srRegion = defaultdict(list)
   with open('sr.roi.bed','r') as IN:
      for line in IN:
         chrom, regionStart, regionEnd, repInfo = line.strip().split()
         repType, totalLen, unitLen, repLen, repBase = repInfo.strip().split("|")
         if repType == 'Simple_repeat':
            repType = 'RepS'
         elif repType == 'Low_complexity':
            repType = 'LowC'
         elif repType == 'Satellite':
            repType = 'SL'
         else:
            repType = 'Other_Repeat'
         
         if isRna:
            totalLen = int(regionEnd) - int(regionStart)
            if totalLen < hpLen:
               continue
            try:
               unitLen_num = float(unitLen)
               repLen = str(totalLen / unitLen_num)
            except ValueError:
               pass
            totalLen = str(totalLen)
         
         srRegion[chrom].append([regionStart, regionEnd, repType, totalLen, unitLen, repLen])
   
   # output variables
   return(srRegion)

#----------------------------------------------------------------------------------------------
# generate locList, where each member is a target site
#------------------------------------------------------------------------------------------------
def getLocList(bedTarget, hpRegion, repRegion, srRegion):
   locList = []
   with open(bedTarget,'r') as IN:
      for line in IN:
         if line.startswith('track name='):
            continue
         lineList = line.strip().split('\t')
         chrom = lineList[0]
         regionStart = int(lineList[1]) + 1   # target region starts from 1-base after 
         regionEnd = lineList[2]
         interval = [] # information for an interval

         pos = regionStart
         lineEnd = False

         while not lineEnd:
            (hpInfo, srInfo, repInfo) = ('.', '.', '.')
            repTypeSet = set()
            # check if the site is in homopolymer region (not including 1 base before) 
            for (regionStart_tmp, regionEnd_tmp, repType_tmp, totalLen_tmp, realL, realR) in hpRegion[chrom]:
               if pos >= int(regionStart_tmp) - 0 and pos <= int(regionEnd_tmp):
                  repTypeSet.add(repType_tmp)
                  hpInfo = ';'.join([chrom, regionStart_tmp, regionEnd_tmp, totalLen_tmp, realL, realR])
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
            interval.append((chrom, str(pos), repType, hpInfo, srInfo, repInfo))

            if str(pos) == regionEnd:
               lineEnd = True
            else:
               pos += 1
         locList.append(interval)
   
   # output variables
   return(locList)
   
   
