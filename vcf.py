import os

# NOTE: temporarily output a VCF for each threshold
tmp_dup_cutoff = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]
rangeDupCutoff = range(len(tmp_dup_cutoff))

#-------------------------------------------------------------------------------------
# assign genotype
#-------------------------------------------------------------------------------------
def assign_gt(alt, chrom, vmf, fltr):
   ''' Function for faking the Genotype i.e. GT field
   for downstream tools
   :param alts (str) alternative allele(s)
   :param chrom (str) chromosome the variant is on
   :param vmf (str) variant minor allele frequency (comma seperated for multi-allelic sites)
   :param fltr (str) FILTER for this variant
   '''
   if fltr == 'LOH_HomRef': # TumorNormal case
      genotype = '0/0'
   else:
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
def biAllelicVar(alleles, RepRegion, outVcf, isDuplex, tumorNormal = False, outVariants = None):
   ID = '.'
   FORMAT = 'GT:AD:VF'

   # duplex-seq runs; temp: multiple VCFs for duplex-seq runs and do not output outVariants
   if isDuplex:
      chrom, pos, ref, alt, typ, dp, vdp, vaf, sumt, svmt, svmf, dumt, dvmt, dvmf, qual, fqual, fltr = alleles[0]
      INFO = ';'.join(['TYPE=' + typ, 'RepRegion=' + RepRegion, 'DP=' + dp, 'sUMT=' + sumt, 'sVMT=' + svmt, 'sVMF=' + svmf, 'dUMT=' + dumt, 'dVMT=' + dvmt, 'dVMF=' + dvmf])
      gt = assign_gt(alt, chrom, svmf, fltr)
      ad = assign_ad(sumt, svmt)         
      SAMPLE = ':'.join([gt, ad, svmf])
      vcfLine = '\t'.join([chrom, pos, ID, ref, alt, qual, fltr, INFO, FORMAT, SAMPLE]) + '\n'
      
      for i in rangeDupCutoff:
         if fqual >= tmp_dup_cutoff[i]:
            outVcf[i].write(vcfLine)
               
   # normal DNA-seq runs; single VCF output
   else:
      # temporarily assume tumor-normal mode for non-duplex runs only
      if tumorNormal:
         chrom, pos, ref, alt, typ, dp, vdp, vaf, umt, vmt, vmf, qual, fqual, fltr, fetpval = alleles[0]
         INFO = ';'.join(
            ['TYPE=' + typ, 'RepRegion=' + RepRegion, 'TNB=' + fetpval,
             'DP=' + dp, 'UMT=' + umt, 'VMT=' + vmt, 'VMF=' + vmf])
        
      # single VCF for normal DNA-seq runs
      else:
         chrom, pos, ref, alt, typ, dp, vdp, vaf, umt, vmt, vmf, qual, fqual, fltr = alleles[0]
         INFO = ';'.join(
            ['TYPE=' + typ, 'RepRegion=' + RepRegion,
             'DP=' + dp, 'UMT=' + umt, 'VMT=' + vmt, 'VMF=' + vmf])

      gt = assign_gt(alt, chrom, vmf, fltr)
      ad = assign_ad(umt, vmt)         
      SAMPLE = ':'.join([gt, ad, vmf])

      if fltr.find('LOH_HomRef') != -1:
         alt = '.'
         fltr = fltr.replace('_HomRef','')

      cutVarLine = '\t'.join([chrom, pos, ref, alt, typ, dp, vdp, vaf, umt, vmt, vmf, qual, fltr]) + '\n'
      vcfLine = '\t'.join([chrom, pos, ID, ref, alt, qual, fltr, INFO, FORMAT, SAMPLE]) + '\n'
      outVcf.write(vcfLine)
      outVariants.write(cutVarLine)

#--------------------------------------------------------------------------------------
# function to handle multi-allelic variants
#--------------------------------------------------------------------------------------
def multiAllelicVar(alleles, RepRegion, outVcf, isDuplex, tumorNormal = False, outVariants = None):
   ID = '.'
   FORMAT = 'GT:AD:VF'

   if tumorNormal:
      fltr_col = -2
   else:
      fltr_col = -1

   tmpAlleles = [x for x in alleles if x[fltr_col] == 'PASS']
   lenTmpAlleles = len(tmpAlleles)

   if lenTmpAlleles == 0:
      pass

   elif lenTmpAlleles == 1:

      if isDuplex:
         chrom, pos, ref, alt, typ, dp, vdp, vaf, sumt, svmt, svmf, dumt, dvmt, dvmf, qual, fqual, fltr = alleles[0]
         INFO = ';'.join(['TYPE=' + typ, 'RepRegion=' + RepRegion, 'DP=' + dp, 'sUMT=' + sumt, 'sVMT=' + svmt, 'sVMF=' + svmf, 'dUMT=' + dumt, 'dVMT=' + dvmt, 'dVMF=' + dvmf])
         gt = assign_gt(alt, chrom, svmf, fltr)
         ad = assign_ad(sumt, svmt)         
         SAMPLE = ':'.join([gt, ad, svmf])
         vcfLine = '\t'.join([chrom, pos, ID, ref, alt, qual, fltr, INFO, FORMAT, SAMPLE]) + '\n'
         
         for i in rangeDupCutoff:
            if fqual >= tmp_dup_cutoff[i]:
               outVcf[i].write(vcfLine)
      else:      
         if tumorNormal:
            chrom, pos, ref, alt, typ, dp, vdp, vaf, umt, vmt, vmf, qual, fqual, fltr, fetpval = alleles[0]
            INFO = ';'.join(
               ['TYPE=' + typ, 'RepRegion=' + RepRegion, 'TNB=' + fetpval,
                'DP=' + dp, 'UMT=' + umt, 'VMT=' + vmt, 'VMF=' + vmf])
         else:
            chrom, pos, ref, alt, typ, dp, vdp, vaf, umt, vmt, vmf, qual, fqual, fltr = tmpAlleles[0]
            INFO = ';'.join(
               ['TYPE=' + typ, 'RepRegion=' + RepRegion,
                'DP=' + dp, 'UMT=' + umt, 'VMT=' + vmt, 'VMF=' + vmf])

         gt = assign_gt(alt, chrom, vmf, fltr)
         ad = assign_ad(umt, vmt)
         SAMPLE = ':'.join([gt, ad, vmf])

         if fltr.find('LOH_HomRef') != -1:
            alt = '.'
            fltr = fltr.replace('_HomRef','')

         vcfLine = '\t'.join([chrom, pos, ID, ref, alt, qual, fltr, INFO, FORMAT, SAMPLE]) + '\n'
         outVcf.write(vcfLine)
         cutVarLine = '\t'.join([chrom, pos, ref, alt, typ, dp, vdp, vaf, umt, vmt, vmf, qual, fltr]) + '\n'
         outVariants.write(cutVarLine)

   else:
      VDPs, VAFs, VMTs, UMTs, VMFs, dVMTs, dUMTs, dVMFs, QUALs, fQUALs, TYPEs, REFs, ALTs, DPs, Pvals = [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
      for allele in tmpAlleles:
         if isDuplex: # duplex-seq runs
            chrom, pos, ref, alt, typ, dp, vdp, vaf, umt, vmt, vmf, dumt, dvmt, dvmf, qual, fqual, fltr = allele
            dUMTs.append(dumt)
            dVMTs.append(dvmt)
            dVMFs.append(dvmf)

         else: # regular dna-seq
            if tumorNormal:
               chrom, pos, ref, alt, typ, dp, vdp, vaf, umt, vmt, vmf, qual, fqual, fltr, fetpval = allele
               Pvals.append(fetpval)
            else:
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
      if all(x == REFs[0] for x in REFs):
         finalRef = REFs[0]
         finalAlt = ','.join(ALTs)
      else:
         # Assumption: the first bases are the same
         finalRef = max(REFs, key = len)
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

      if isDuplex:
         # debug check
         assert len(set(dUMTs)) == 1, "The number of used duplex UMIs at a site should be the same across all alleles"

         alldVMTs = ','.join(dVMTs)
         alldVMFs = ','.join(dVMFs)
         
         INFO = ';'.join(['TYPE=' + typ, 'RepRegion=' + RepRegion, 'DP=' + dp, 'sUMT=' + umt, 'sVMT=' + vmt, 'sVMF=' + vmf, 'dUMT=' + dumt, 'dVMT=' + dvmt, 'dVMF=' + dvmf])
         gt = assign_gt(alt, chrom, vmf, fltr)
         ad = assign_ad(umt, vmt)
         SAMPLE = ':'.join([gt, ad, vmf])
         vcfLine = '\t'.join([chrom, pos, ID, ref, alt, qual, fltr, INFO, FORMAT, SAMPLE]) + '\n'
         
         for i in rangeDupCutoff:
            if fqual >= tmp_dup_cutoff[i]:
               outVcf[i].write(vcfLine)
         
      else:
         if tumorNormal:
            allPvals = ','.join(Pvals)
            INFO = ';'.join(['TYPE=' + allTypes, 'RepRegion=' + RepRegion,
                             'TNB=' + allPvals,
                             'DP='+ allDPs,'UMT=' + umt,'VMT=' + allVMTs,'VMF=' + allVMFs])
         else:
            INFO = ';'.join(['TYPE=' + allTypes, 'RepRegion=' + RepRegion,
                             'DP=' + allDPs, 'UMT=' + umt, 'VMT=' + allVMTs, 'VMF=' + allVMFs])

         # debug check
         assert fltr == 'PASS', "Error in assinging FILTER at Multi-Allelic Site"

         FORMAT = 'GT:AD:VF'
         gt = assign_gt(finalAlt, chrom, allVMFs, fltr)
         ad = assign_ad(umt, allVMTs)
         SAMPLE = ':'.join([gt, ad, allVMFs])
         vcfLine = '\t'.join([chrom, pos, ID, finalRef, finalAlt, newQual, 'PASS', INFO, FORMAT, SAMPLE]) + '\n'
         outVcf.write(vcfLine)
         cutVarLine = '\t'.join([chrom, pos, finalRef, finalAlt, allTypes, allDPs, allVDPs, allVAFs, umt, allVMTs, allVMFs, newQual,'PASS']) + '\n'
         outVariants.write(cutVarLine)

#--------------------------------------------------------------------------------------
# main function
#--------------------------------------------------------------------------------------
def makeVcf(runPath, outlong, sampleName, refGenome, isDuplex, tumorNormal = False):
   # change working directory to runDir
   os.chdir(runPath)

   ID = '.'   
   headerVariants = ['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'DP', 'VDP', 'VAF', 'sUMT', 'sVMT', 'sVMF', 'QUAL', 'FILTER']
   
   if isDuplex:
      if not os.path.exists('vcf'):
         os.makedirs('vcf')
         
      outVariants = None
      outAll = open(sampleName + '.smCounter-duplex.all.txt', 'w')
      outVcf = []
      for i in rangeDupCutoff:
         vcfName = 'vcf/' + sampleName + '.smCounter-duplex.cut.' + str(tmp_dup_cutoff[i]) + '.vcf'
         outVcf.append(open(vcfName, 'w'))
      
      # header for .all.txt output
      headerAll = ['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'sUMT', 'sVMT', 'sVMF', 'dUMT', 'dVMT', 'dVMF', 'DP', 'VDP', 'VAF', 'RefForPrimer', 'RefRevPrimer', 'primerOR', 'pLowQ', 'hqUmiEff', 'allUmiEff', 'refMeanRpb', 'altMeanRpb', 'rpbEffectSize', 'repType', 'hpInfo', 'simpleRepeatInfo', 'tandemRepeatInfo', 'FR', 'MT', 'UFR', 'sUMT_A', 'sUMT_T', 'sUMT_G', 'sUMT_C', 'dUMT_A', 'dUMT_T', 'dUMT_G', 'dUMT_C', 'plowDupVMF', 'logLH1', 'logLR', 'FILTER']
      
   else:
      outAll = open(sampleName + '.smCounter.all.txt', 'w')
      outVcf = open(sampleName + '.smCounter.cut.vcf','w')
      outVariants = open(sampleName + '.smCounter.cut.txt','w')
      outLowPi = open(sampleName + '.smCounter.lowQ.txt','w')
      
      # header for .all.txt output
      headerAll = ['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'sUMT', 'sForUMT', 'sRevUMT', 'sVMT', 'sForVMT', 'sRevVMT', 'sVMF', 'sForVMF', 'sRevVMF', 'VDP', 'VAF', 'RefForPrimer', 'RefRevPrimer', 'primerOR', 'pLowQ', 'hqUmiEff', 'allUmiEff', 'refMeanRpb', 'altMeanRpb', 'rpbEffectSize', 'repType', 'hpInfo', 'simpleRepeatInfo', 'tandemRepeatInfo', 'DP', 'FR', 'MT', 'UFR', 'sUMT_A', 'sUMT_T', 'sUMT_G', 'sUMT_C', 'logpval', 'FILTER']
      
      # default(6.0) and relaxed(2.5) cutoff
      cutoff = 6.0
      minCutoff = {'INDEL':2.5, 'SNP':2.5} ## Cutoff for the low-threshold file

      if tumorNormal:
         headerAll.append('TNB')
         headerLowPi = ['READ_SET', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'TYPE', 'RepRegion', 'TNB', 'DP', 'UMT', 'VMT', 'VMF']
      else:
         headerLowPi = ['READ_SET', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'TYPE', 'RepRegion', 'DP', 'UMT', 'VMT', 'VMF']

   
   headerVcf = '##fileformat=VCFv4.2' + '\n' + \
      '##reference={ref}'.format(ref=refGenome) + '\n' + \
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
      '##FILTER=<ID=PrimerCP,Description="variant is clustered within 2 bases from primer sequence due to possible primer dimers">' + '\n'
   
   if isDuplex:
      headerVcf = headerVcf + \
         '##FILTER=<ID=LowDupVMF,Description="Duplex VMF is significantly lower than singleplex VMF">' + '\n' 
   
   if tumorNormal:
      headerVcf = headerVcf + \
         '##FILTER=<ID=LOH,Description="Loss of Heterozygocity">' + '\n' \
         '##FILTER=<ID=Germline_Risk,Description="Not a significant difference in UMI counts for the variant allele between Normal and Tumor Samples">' + '\n'
   
   # common INFO fields
   headerVcf = headerVcf + \
      '##INFO=<ID=TYPE,Number=.,Type=String,Description="Variant type: SNP/INDEL/COMPLEX">' + '\n' + \
      '##INFO=<ID=RepRegion,Number=.,Type=String,Description="Repetitive region">' + '\n' + \
      '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">' + '\n'

   if isDuplex:
      headerVcf = headerVcf + \
         '##INFO=<ID=sUMT,Number=1,Type=Integer,Description="Total used singleplex UMI depth">' + '\n' + \
         '##INFO=<ID=sVMT,Number=.,Type=Integer,Description="Variant singleplex UMI depth">' + '\n' + \
         '##INFO=<ID=sVMF,Number=.,Type=Float,Description="Variant singleplex UMI allele frequency">' + '\n' \
         '##INFO=<ID=dUMT,Number=1,Type=Integer,Description="Total used duplex UMI depth">' + '\n' + \
         '##INFO=<ID=dVMT,Number=.,Type=Integer,Description="Variant duplex UMI depth">' + '\n' + \
         '##INFO=<ID=dVMF,Number=.,Type=Float,Description="Variant duplex UMI allele frequency">' + '\n' 
   else:
      headerVcf = headerVcf + \
         '##INFO=<ID=UMT,Number=1,Type=Integer,Description="Total used UMI depth">' + '\n' + \
         '##INFO=<ID=VMT,Number=.,Type=Integer,Description="Variant UMI depth">' + '\n' + \
         '##INFO=<ID=VMF,Number=.,Type=Float,Description="Variant UMI allele frequency">' + '\n' 

   if tumorNormal:
      headerVcf = headerVcf + \
         '##INFO=<ID=TNB,Number=.,Type=Float,Description="FDR Corrected Phred-scaled p-value using Fisher\'s exact test to detect Tumor Normal Bias">' + '\n'

   # common GT tag
   headerVcf = headerVcf + \
      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' + '\n'
   
   if isDuplex:
      headerVcf = headerVcf + \
         '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Filtered allelic sUMI depths for the ref and alt alleles">' + '\n' + \
         '##FORMAT=<ID=VF,Number=.,Type=Float,Description="Variant sUMI allele frequency, same as sVMF">' + '\n' + \
         '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + [sampleName]) + '\n'
   else:
      headerVcf = headerVcf + \
         '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Filtered allelic UMI depths for the ref and alt alleles">' + '\n' + \
         '##FORMAT=<ID=VF,Number=.,Type=Float,Description="Variant UMI allele frequency, same as VMF">' + '\n' + \
         '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + [sampleName]) + '\n'\
         
   alleles = []
   lastCHROM, lastPOS = '', ''

   # write headers
   outAll.write('\t'.join(headerAll) + '\n')
   if isDuplex:
      for i in rangeDupCutoff:
         outVcf[i].write(headerVcf)
   else:
      outVariants.write('\t'.join(headerVariants) + '\n')
      outLowPi.write('\t'.join(headerLowPi) + '\n')
      outVcf.write(headerVcf)

   cnt = 1
   with open(outlong, 'r') as f:
      next(f)
      for line in f:
         tempLine = line
         if tempLine.find('LOH_HomRef') != -1:
            tempLine = tempLine.replace('_HomRef','')
         outAll.write(tempLine)
         cnt += 1
         
         if isDuplex:
            CHROM, POS, REF, ALT, TYPE, sUMT, sVMT, sVMF, dUMT, dVMT, dVMF, DP, VDP, VAF, RefForPrimer, RefRevPrimer, primerOR, pLowQ, hqUmiEff, allUmiEff, refMeanRpb, altMeanRpb, rpbEffectSize, repType, hpInfo, simpleRepeatInfo, tandemRepeatInfo, FR, MT, UFR, sUMT_A, sUMT_T, sUMT_G, sUMT_C, dUMT_A, dUMT_T, dUMT_G, dUMT_C, plowDupVMF, logLH1, logLR, FILTER = line.strip().split()
         elif tumorNormal:
            CHROM, POS, REF, ALT, TYPE, sUMT, sForUMT, sRevUMT, sVMT, sForVMT, sRevVMT, sVMF, sForVMF, sRevVMF, VDP, VAF, RefForPrimer, RefRevPrimer, primerOR, pLowQ, hqUmiEff, allUmiEff, refMeanRpb, altMeanRpb, rpbEffectSize, repType, hpInfo, simpleRepeatInfo, tandemRepeatInfo, DP, FR, MT, UFR, sUMT_A, sUMT_T, sUMT_G, sUMT_C, logpval, FILTER, TNFetPval = line.strip().split('\t')
         else:
            CHROM, POS, REF, ALT, TYPE, sUMT, sForUMT, sRevUMT, sVMT, sForVMT, sRevVMT, sVMF, sForVMF, sRevVMF, VDP, VAF, RefForPrimer, RefRevPrimer, primerOR, pLowQ, hqUmiEff, allUmiEff, refMeanRpb, altMeanRpb, rpbEffectSize, repType, hpInfo, simpleRepeatInfo, tandemRepeatInfo, DP, FR, MT, UFR, sUMT_A, sUMT_T, sUMT_G, sUMT_C, logpval, FILTER = line.strip().split('\t')
         
         if TYPE == '0':
            continue
         
         if ALT == 'DEL': 
            continue
         
         if isDuplex:
            QUAL = logLR if logLR != 'NA' else '0.00'
         else:
            QUAL = logpval if logpval != 'NA' else '0.00'

         try:
            fQUAL = float(QUAL)
         except ValueError:
            fQUAL = 0.00

         minThr = tmp_dup_cutoff[0] if isDuplex else minCutoff[TYPE.upper()]
         if fQUAL < minThr and FILTER.find('LOH_HomRef') == -1: # let LOH_HomRef variants in tumor through
            lastCHROM = '.'
            continue
         try:
            VAF = str(float(VAF) / 100)
         except ValueError:
            VAF = '-1'
         try:
            sVMF = str(float(sVMF) / 100)
         except ValueError:
            sVMF = '-1'
            
         if isDuplex:
            try:
               dVMF = str(float(dVMF) / 100)
            except ValueError:
               dVMF = '-1'

         # rep types are separeted by ";" in the long output. Replace to "," to comply with VCF format
         RepRegion = repType.replace(';', ',')
         
         lenAlleles = len(alleles)
         
         if isDuplex:
            currentAllele = (CHROM, POS, REF, ALT, TYPE, DP, VDP, VAF, sUMT, sVMT, sVMF, dUMT, dVMT, dVMF, QUAL, fQUAL, FILTER)
         elif tumorNormal:
            currentAllele = (CHROM, POS, REF, ALT, TYPE, DP, VDP, VAF, sUMT, sVMT, sVMF, QUAL, fQUAL, FILTER, TNFetPval)
            # tempVar = (CHROM, POS, REF, ALT, TYPE, DP, VDP, VAF, sUMT, sVMT, sVMF, QUAL, fQUAL, FILTER, TNFetPval)  
         else:
            currentAllele = (CHROM, POS, REF, ALT, TYPE, DP, VDP, VAF, sUMT, sVMT, sVMF, QUAL, fQUAL, FILTER)
            # tempVar = (CHROM, POS, REF, ALT, TYPE, DP, VDP, VAF, sUMT, sVMT, sVMF, QUAL, FILTER)
         
         ## Write to low-PI file, want the LOH_HomRef variants to appear in the cut files; temp: for non-duplex runs only
         if not isDuplex and fQUAL < cutoff and FILTER.find('LOH_HomRef') == -1: 
            if tumorNormal:
               outLowPi.write('\t'.join([sampleName,CHROM,POS,".", REF, ALT, QUAL, FILTER, TYPE, RepRegion, TNFetPval, DP, sUMT, sVMT, sVMF]))
               outLowPi.write('\n')
            else:
               outLowPi.write('\t'.join([sampleName, CHROM, POS, ".", REF, ALT, QUAL, FILTER, TYPE, RepRegion, DP, sUMT, sVMT, sVMF]))
               outLowPi.write('\n')
            
            continue
            
         # if current chrom and position equal to last line, append it for potential multi-allelic output
         if lenAlleles == 0 or (CHROM == lastCHROM and POS == lastPOS):
            alleles.append(currentAllele)

         # for new chrom or position, if last variant is not multi-allelic, write to vcf directly
         elif lenAlleles == 1:
            biAllelicVar(alleles, RepRegion, outVcf, isDuplex, tumorNormal, outVariants)
            alleles = [currentAllele]

         # if last variant is possible multi-allelic, combine and write as one 
         else:            
            multiAllelicVar(alleles, RepRegion, outVcf, isDuplex, tumorNormal, outVariants)
            alleles = [currentAllele]

         lastCHROM, lastPOS = CHROM, POS

   # take care of the last line
   lenAlleles = len(alleles)
   if lenAlleles == 1:
      biAllelicVar(alleles, RepRegion, outVcf, isDuplex, tumorNormal, outVariants)
   elif lenAlleles >= 2:
      multiAllelicVar(alleles, RepRegion, outVcf, isDuplex, tumorNormal, outVariants)
   else:
      pass

   # close all output file handles
   outAll.close()
   if isDuplex:
      for vcf in outVcf:
         vcf.close()
   else:
      outVariants.close()
      outVcf.close()
      outLowPi.close()
