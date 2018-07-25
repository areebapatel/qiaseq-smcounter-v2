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
def findhp(bedName, outName, minLength,refg):
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
      
      if (start - 1 - 100) < 0:
         start_coord = start
      else:
         start_coord = start - 1 - 100
         
      origRef = refseq.fetch(reference=chrom, start=start_coord, end=end + 100)
      origRef = origRef.upper()

      hpL = 0
      for i in range(len(origRef))[1:]:
         if origRef[i] == origRef[i-1]:
            continue
         else:
            hpLen = i - hpL 
            realL = hpL - 1 + start - 100
            realR = i - 1  + start - 100
            if hpLen >= minLength and realL <= end and realR >= start:
               outline = '\t'.join([chrom, str(max(realL, start)), str(min(realR, end)), 'HP', str(hpLen), str(realL), str(realR), origRef[hpL]]) + '\n'
               outfile.write(outline)
            hpL = i 


   outfile.close()

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
   
def assign_gt(alt,chrom,vmf):
   ''' Function for faking the Genotype i.e. GT field
   for downstream tools
   :param alts (str) alternative allele(s)
   :param chrom (str) chromosome the variant is on
   :param vmf (str) variant minor allele frequency (comma seperated for multi-allelic sites)
   '''
   alts = alt.split(",")
   if len(alts) >= 2: ## Treat all multiallelic sites as heterozygotes for the first 2 variant alleles
      genotype = '1/2'
   elif chrom == "chrY" or chrom == "chrM":
      genotype = '1'
   elif float(vmf) > 0.95 : ## Treat as Heterozygous
      genotype = '1/1'
   else:
      genotype = '0/1'

   return genotype

def assign_ad(uumi,vumi):
   ''' Function for faking the Allele Depth i.e. AD field
   for downstream tools
   :param uumi (str) total umis at the variant site
   :param vumi (str) umis corresponding to the non-reference allele(s) at the variant site (comma seperated for multi-allelic sites)
   '''    
   vumis = vumi.split(',')   
   refumi = int(uumi)
   for umi in vumis:
      refumi = refumi - int(umi)
   refumi = str(refumi)
   ad = refumi + ',' + ','.join(vumis)
   return ad 
      
#--------------------------------------------------------------------------------------
# function to handle normal variants
#--------------------------------------------------------------------------------------
def biAllelicVar(alleles, RepRegion, outVcf, outVariants):
   ID = '.'
   chrom, pos, ref, alt, typ, dp, vdp, vaf, umt, vmt, vmf, qual, fqual, fltr = alleles[0]

   INFO = ';'.join(
      ['TYPE=' +typ,'RepRegion=' + RepRegion,'DP='+dp,'UMT='+umt,'VMT='+vmt,
      'VMF='+vmf]
      ) 
   FORMAT = 'GT:AD:VF'
   gt = assign_gt(alt,chrom,vmf)      
   ad = assign_ad(umt,vmt)         
   SAMPLE = ':'.join([gt,ad,vmf])
   vcfLine = '\t'.join([chrom, pos, ID, ref, alt, qual, fltr, INFO, FORMAT, SAMPLE]) + '\n'
   outVcf.write(vcfLine)
   cutVarLine = '\t'.join([chrom, pos, ref, alt, typ, dp, vdp, vaf, umt, vmt, vmf, qual, fltr]) + '\n'
   outVariants.write(cutVarLine)

#--------------------------------------------------------------------------------------
# function to handle multi-allelic variants
#--------------------------------------------------------------------------------------
def multiAllelicVar(alleles, RepRegion, outVcf, outVariants):
  ID = '.'
  tmpAlleles = [x for x in alleles if x[-1] == 'PASS']
  lenTmpAlleles = len(tmpAlleles)
  if lenTmpAlleles == 0:
     pass
  elif lenTmpAlleles == 1:
     chrom, pos, ref, alt, typ, dp, vdp, vaf, umt, vmt, vmf, qual, fqual, fltr = tmpAlleles[0]
     INFO = ';'.join(
            ['TYPE=' +typ,'RepRegion=' + RepRegion,'DP='+dp,'UMT='+umt,'VMT='+vmt,
             'VMF='+vmf]
         )
     FORMAT = 'GT:AD:VF'
     gt = assign_gt(alt,chrom,vmf)
     ad = assign_ad(umt,vmt)
     SAMPLE = ':'.join([gt,ad,vmf])     
     vcfLine = '\t'.join([chrom, pos, ID, ref, alt, qual, fltr, INFO, FORMAT, SAMPLE]) + '\n'
     outVcf.write(vcfLine)
     cutVarLine = '\t'.join([chrom, pos, ref, alt, typ, dp, vdp, vaf, umt, vmt, vmf, qual, fltr]) + '\n'
     outVariants.write(cutVarLine)
  else:
     VDPs, VAFs, VMTs, UMTs, VMFs, QUALs, fQUALs, TYPEs, REFs, ALTs, DPs = [], [], [], [], [], [], [], [], [], [], []
     for allele in tmpAlleles:
        chrom, pos, ref, alt, typ, dp, vdp, vaf, umt, vmt, vmf, qual, fqual, fltr = allele
        VDPs.append(vdp)
        VAFs.append(vaf)
        VMTs.append(vmt)
        UMTs.append(umt)
        VMFs.append(vmf)
        QUALs.append(qual)
        TYPEs.append(typ)
        REFs.append(ref)
        ALTs.append(alt)
        fQUALs.append(fqual)
        DPs.append(dp)

     # debug check
     assert len(set(UMTs)) == 1, "The number of used UMIs at a site should be the same across all alleles"

     # align multiple alleles to the same REF if necessary
     if all(x==REFs[0] for x in REFs):
        finalRef = REFs[0]
        finalAlt = ','.join(ALTs)
     else:
        # Assumption: the first bases are the same 
        finalRef = max(REFs, key=len)
        for j in range(len(ALTs)):
           ALTs[j] = ALTs[j] if REFs[j] == finalRef else ALTs[j] + finalRef[len(REFs[j]):]
        finalAlt = ','.join(ALTs)

     newQual = str(min(fQUALs))
     allTypes = ','.join(TYPEs)
     allVDPs = ','.join(VDPs)
     allVAFs = ','.join(VAFs)
     allVMTs = ','.join(VMTs)
     allVMFs = ','.join(VMFs)
     allDPs = ','.join(DPs)

     INFO = ';'.join(
            ['TYPE=' +allTypes,'RepRegion=' + RepRegion,'DP='+allDPs,'UMT='+umt,'VMT='+allVMTs,
             'VMF='+allVMFs]
         )          
     FORMAT = 'GT:AD:VF' 
     gt = assign_gt(finalAlt,chrom,allVMFs)
     ad = assign_ad(umt,allVMTs)   
     SAMPLE = ':'.join([gt,ad,allVMFs])
     vcfLine = '\t'.join([chrom, pos, ID, finalRef, finalAlt, newQual, 'PASS', INFO, FORMAT, SAMPLE]) + '\n'     
     outVcf.write(vcfLine)
     cutVarLine = '\t'.join([chrom, pos, finalRef, finalAlt, allTypes, allDPs, allVDPs, allVAFs, umt, allVMTs, allVMFs, newQual,'PASS']) + '\n'
     outVariants.write(cutVarLine)

#--------------------------------------------------------------------------------------
# main function
#--------------------------------------------------------------------------------------
def makeVcf(runPath, outlong, sampleName):
   # change working directory to runDir
   os.chdir(runPath)
   outAll = open(sampleName + '.smCounter.all.txt', 'w')
   outVariants = open(sampleName + '.smCounter.cut.txt','w')
   outVcf = open(sampleName + '.smCounter.cut.vcf','w')
   outLowPi = open(sampleName + '.smCounter.lowQ.txt','w')
   
   cutoff = 6
   minCutoff = {'INDEL': 2,'SNP':2} ## Cutoff for the low-PI file
   
   ID = '.'
   headerAll = ['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'sUMT', 'sForUMT', 'sRevUMT', 'sVMT', 'sForVMT', 'sRevVMT', 'sVMF', 'sForVMF', 'sRevVMF', 'VDP', 'VAF', 'RefForPrimer', 'RefRevPrimer', 'primerOR', 'pLowQ', 'hqUmiEff', 'allUmiEff', 'refMeanRpb', 'altMeanRpb', 'rpbEffectSize', 'repType', 'hpInfo', 'simpleRepeatInfo', 'tandemRepeatInfo', 'DP', 'FR', 'MT', 'UFR', 'sUMT_A', 'sUMT_T', 'sUMT_G', 'sUMT_C', 'logpval', 'FILTER']
   headerVariants = ['CHROM','POS','REF','ALT','TYPE','DP','VDP','VAF','sUMT','sVMT','sVMF','QUAL','FILTER']
   headerLowPi = [sampleName] + headerVariants   
   headerVcf = '##fileformat=VCFv4.2' + '\n' + \
         '##reference=GRCh37' + '\n' + \
         '##FILTER=<ID=LM,Description="Low coverage (fewer than 5 barcodes)">' + '\n' + \
         '##FILTER=<ID=RepT,Description="Variant in tandem repeat (TFR) regions">' + '\n' + \
         '##FILTER=<ID=RepS,Description="Variant in simple repeats (RepeatMasker) regions">' + '\n' + \
         '##FILTER=<ID=HP,Description="Inside or flanked by homopolymer regions">' + '\n' + \
         '##FILTER=<ID=LowC,Description="Variant in Low complexity regions, as defined in RepeatMasker">' + '\n' + \
         '##FILTER=<ID=SL,Description="Variant in micro-satelite regions, as defined in RepeatMasker">' + '\n' + \
         '##FILTER=<ID=SB,Description="Strand Bias">' + '\n' + \
         '##FILTER=<ID=DP,Description="Too many discordant pairs">' + '\n' + \
         '##FILTER=<ID=MM,Description="Too many mismatches in a read. Default threshold is 6.5 per 100 bases">' + '\n' + \
         '##FILTER=<ID=LowQ,Description="Low base quality">' + '\n' + \
         '##FILTER=<ID=RBCP,Description="Variant are clustered at the end of barcode-side reads">' + '\n' + \
         '##FILTER=<ID=RPCP,Description="Variant are clustered at the end of primer-side reads">' + '\n' + \
         '##FILTER=<ID=PB,Description="Primer bias filter. odds ratio > 10 or < 0.1">' + '\n' + \
         '##FILTER=<ID=PrimerCP,Description="variant is clustered within 2 bases from primer sequence due to possible primer dimers">' + '\n' + \
         '##INFO=<ID=TYPE,Number=.,Type=String,Description="Variant type: SNP/INDEL/COMPLEX">' + '\n' + \
         '##INFO=<ID=RepRegion,Number=.,Type=String,Description="Repetitive region">' + '\n' + \
         '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">' + '\n' + \
         '##INFO=<ID=UMT,Number=1,Type=Integer,Description="Total used UMI depth">' + '\n' + \
         '##INFO=<ID=VMT,Number=.,Type=Integer,Description="Variant UMI depth">' + '\n' + \
         '##INFO=<ID=VMF,Number=.,Type=Float,Description="Variant UMI allele frequency">' + '\n' + \
         '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' + '\n' + \
         '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Filtered allelic MT depths for the ref and alt alleles">' + '\n' + \
         '##FORMAT=<ID=VF,Number=.,Type=Float,Description="Variant UMI allele frequency, same as VMF">' + '\n' + \
         '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + [sampleName]) + '\n'

   alleles = []
   lastCHROM, lastPOS = '', ''

   outAll.write('\t'.join(headerAll)+'\n')
   outVariants.write('\t'.join(headerVariants)+'\n')
   outLowPi.write('\t'.join(headerLowPi)+'\n')
   outVcf.write(headerVcf)

   cnt = 1
   with open(outlong, 'r') as f:
      next(f)
      for line in f:
         outAll.write(line)
         cnt += 1
         CHROM, POS, REF, ALT, TYPE, sUMT, sForUMT, sRevUMT, sVMT, sForVMT, sRevVMT, sVMF, sForVMF, sRevVMF, VDP, VAF, RefForPrimer, RefRevPrimer, primerOR, pLowQ, hqUmiEff, allUmiEff, refMeanRpb, altMeanRpb, rpbEffectSize, repType, hpInfo, simpleRepeatInfo, tandemRepeatInfo, DP, FR, MT, UFR, sUMT_A, sUMT_T, sUMT_G, sUMT_C, logpval, FILTER = line.strip().split('\t')
         
         if TYPE == '0':
            continue
         
         if ALT == 'DEL': 
            continue

         QUAL = logpval if logpval != 'NA' else '0.00'

         try:
            fQUAL = float(QUAL)
         except ValueError:
            fQUAL = 0.00

         if fQUAL < minCutoff[TYPE.upper()]:
            lastCHROM = '.'
            continue
         try:
            VAF = str(float(VAF)/100)
         except ValueError:
            VAF = '-1'
         try:
            sVMF = str(float(sVMF)/100)
         except ValueError:
            sVMF = '-1'

         # rep types are separeted by ";" in the long output. Replace to "," to comply with VCF format
         RepRegion = repType.replace(';', ',')
         
         currentAllele = (CHROM, POS, REF, ALT, TYPE, DP, VDP, VAF, sUMT, sVMT, sVMF, QUAL, fQUAL, FILTER)
         tempVar = (CHROM, POS, REF, ALT, TYPE, DP, VDP, VAF, sUMT, sVMT, sVMF, QUAL, FILTER)
         lenAlleles = len(alleles)

         if fQUAL < cutoff: ## Write to low-PI file
            outLowPi.write(sampleName+'\t'+'\t'.join(tempVar)+'\n')
            continue
            
         # if current chrom and position equal to last line, append it for potential multi-allelic output
         if lenAlleles == 0 or (CHROM == lastCHROM and POS == lastPOS):
            alleles.append(currentAllele)

         # for new chrom or position, if last variant is not multi-allelic, write to vcf directly
         elif lenAlleles == 1:
            biAllelicVar(alleles, RepRegion, outVcf, outVariants)
            alleles = [currentAllele]

         # if last variant is possible multi-allelic, combine and write as one 
         else:            
            multiAllelicVar(alleles, RepRegion, outVcf, outVariants)
            alleles = [currentAllele]

         lastCHROM, lastPOS = CHROM, POS

   # take care of the last line
   lenAlleles = len(alleles)
   if lenAlleles == 1:
      biAllelicVar(alleles, RepRegion, outVcf, outVariants)
   elif lenAlleles >= 2:
      multiAllelicVar(alleles, RepRegion, outVcf, outVariants)
   else:
      pass

   # close all output file handles
   outAll.close()
   outVariants.close()
   outVcf.close()
   outLowPi.close()

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
   findhp(bedTarget, 'hp.roi.bed', str(findHpLen), refGenome, seqType)
   
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
            locList.append((chrom, str(pos), repType, hpInfo, srInfo, repInfo))

            if str(pos) == regionEnd:
               lineEnd = True
            else:
               pos += 1
   
   # output variables
   return(locList)
   
   