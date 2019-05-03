#!/usr/bin/python

import os
import subprocess
import multiprocessing
import datetime
import argparse
import functools

# our modules
import utils
import vcf
from vc import vc_wrapper

# global constants
codePath = os.path.dirname(os.path.abspath(__file__))
pValCode_sin = os.path.join(codePath,'get_pvalue.R')
pValCode_dup = os.path.join(codePath,'get_lr.duplex.R')
locChunkLen = 1000
seed = 10262016
nsim = 5000
parser = None
header_sin = ['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'sUMT', 'sForUMT', 'sRevUMT', 'sVMT', 'sForVMT', 'sRevVMT', 'sVMF', 'sForVMF', 'sRevVMF', 'VDP', 'VAF', 'RefForPrimer', 'RefRevPrimer', 'primerOR', 'pLowQ', 'hqUmiEff', 'allUmiEff', 'refMeanRpb', 'altMeanRpb', 'rpbEffectSize', 'repType', 'hpInfo', 'simpleRepeatInfo', 'tandemRepeatInfo', 'DP', 'FR', 'MT', 'UFR', 'sUMT_A', 'sUMT_T', 'sUMT_G', 'sUMT_C', 'FILTER']
header_dup = ['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'sUMT', 'sVMT', 'sVMF', 'dUMT', 'dVMT', 'dVMF', 'DP', 'VDP', 'VAF', 'RefForPrimer', 'RefRevPrimer', 'primerOR', 'pLowQ', 'hqUmiEff', 'allUmiEff', 'refMeanRpb', 'altMeanRpb', 'rpbEffectSize', 'repType', 'hpInfo', 'simpleRepeatInfo', 'tandemRepeatInfo', 'FR', 'MT', 'sForUMT', 'sRevUMT', 'sForVMT', 'sRevVMT', 'UFR', 'sUMT_A', 'sUMT_T', 'sUMT_G', 'sUMT_C', 'dUMT_A', 'dUMT_T', 'dUMT_G', 'dUMT_C', 'FILTER']
header_bkg = ['CHROM', 'POS', 'REF', 'A/G', 'G/A', 'C/T', 'T/C', 'A/C', 'C/A', 'A/T', 'T/A', 'C/G', 'G/C', 'G/T', 'T/G', 'negStrand', 'posStrand', 'AllSMT', 'workflow' ]

#----------------------------------------------------------------------------------------------
# global for argument parsing (hack that works when calling from either command line or pipeline)
# this is done inside a function because multiprocessing module imports the script
#------------------------------------------------------------------------------------------------
def argParseInit():  
   global parser
   parser = argparse.ArgumentParser(description = 'smCounter2: variant calling using Unique Molecular Identifiers')
   parser.add_argument('--runPath', default = None, help = 'path to working directory')
   parser.add_argument('--bedTarget', default = None, help = 'BED file')
   parser.add_argument('--bamFile', default = None, help = 'BAM file')
   parser.add_argument('--outPrefix', default = None, help = 'file name prefix')
   parser.add_argument('--nCPU', type = int, default = 1, help = 'number of CPU to use in parallel')
   parser.add_argument('--minBQ', type = int, default = 25, help = 'minimum base quality allowed for analysis')
   parser.add_argument('--minMQ', type = int, default = 50, help = "minimum mapping quality allowed for analysis. If the bam is tagged with its mate's mapq, then the minimum of the R1 and R2 mapq will be used for comparison, if not each read is compared independently.")
   parser.add_argument('--hpLen', type = int, default = 10, help = 'minimum length for homopolymers')
   parser.add_argument('--mismatchThr', type = float, default = 6.0, help = 'average number of mismatches per 100 bases allowed')
   parser.add_argument('--primerDist', type = int, default = 2, help = 'filter variants that are within X bases to primer')
   parser.add_argument('--consThr', type = float, default = 0.8, help = 'threshold on read proportion to determine UMI level consensus')
   parser.add_argument('--rpu', type = float, default = 0.0, help = 'mean read pairs per UMI; default at 0 and will be calculated')
   parser.add_argument('--isRna', action = 'store_true', help = 'RNAseq varinat calling only; default is DNAseq')
   parser.add_argument('--primerSide', type = int, default = 1, help = 'read end that includes the primer; default is 1')
   parser.add_argument('--umiTag', type = str, default = 'Mi', help = 'tag name for normal UMI')
   parser.add_argument('--primerTag', type = str, default = 'pr', help = 'tag name for Primer')
   parser.add_argument('--mqTag', type = str, default = 'MQ', help = 'tag name for MapQ score of mate')
   parser.add_argument('--tagSeparator', type = str, default = '-', help = 'tag seperator for splitting UMI tag') 
   parser.add_argument('--minAltUMI', type = int, default = 3, help = 'minimum requirement of ALT UMIs; default is 3')
   parser.add_argument('--maxAltAllele', type = int, default = 2, help = 'maximum ALT alleles that meet minAltUMI to be reported; default is 2 (tri-allelic variants)')
   parser.add_argument('--refGenome',type = str,help = 'path to the reference fasta file')
   parser.add_argument('--repBed',type = str,help = 'path to the simpleRepeat bed file')
   parser.add_argument('--srBed',type = str,help = 'path to the full repeat bed file')
   parser.add_argument('--ds', type = int, default = 10000, help = 'down sample if number of UMIs greater than this value (RNA only)')
   parser.add_argument('--bamType', type = str, default = 'raw', help = 'raw (default): raw BAM file with UMIs; consensus: consensused BAM file')
   parser.add_argument('--inputVCF', type = str, default = None, help = 'optional input VCF file') 
   parser.add_argument('--sinBkgErrorDist', type = str, default = '/srv/qgen/data/annotation/bkg.error.v2.RData', help = 'background error rate distribution for normal DNA-seq runs') 
   parser.add_argument('--isDuplex', action = 'store_true', help = 'duplex-seq DNA varinat calling only; default is normal DNAseq')
   parser.add_argument('--duplexTag', type = str, default = None, help = 'tag name for duplex UMI')
   parser.add_argument('--minRpu', type = int, default = 2, help = 'minimum read pairs for UMI to be included; default is 1 for normal DNA-seq (before gradually dropping singletons) and 2 for duplex-seq')
   parser.add_argument('--dupBkgErrorDist', type = str, default ='/srv/qgen/data/annotation/duplex.bkg.error.ga.RData', help = 'background G>A error rate distribution for duplex-seq runs')  
   
#--------------------------------------------------------------------------------------
# main function
#--------------------------------------------------------------------------------------
def main(args):
   # log run start
   timeStart = datetime.datetime.now()
   print("started " + str(timeStart))
   
   # if argument parser global not assigned yet, initialize it
   if parser == None:
      argParseInit()
   
   # get arguments passed in via a lambda object (e.g. from upstream pipeline)
   if type(args) is not argparse.Namespace:
      argsList = []
      for argName, argVal in args.iteritems():
         if argName == "isDuplex":
            if argVal: # --isDuplex triggered upstream
               argsList.append("--{0}".format(argName))
         else:
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
   
   # convert VCF to BED if inputVCF is not 'none'
   bedTarget = args.bedTarget if args.inputVCF is None else utils.vcf2bed(args.inputVCF)
   
   # gather repetitive regions information
   hpRegion = utils.getHpInfo(bedTarget, args.refGenome, args.isRna, args.hpLen)
   repRegion = utils.getTrInfo(bedTarget, args.repBed, args.isRna, args.hpLen)
   srRegion = utils.getOtherRepInfo(bedTarget, args.srBed, args.isRna, args.hpLen)
   
   # read in bed file and create a list of positions, annotated with repetitive region
   locList = utils.getLocList(bedTarget, hpRegion, repRegion, srRegion, args.isDuplex)
   
   # calculate rpu if args.rpu = 0
   if args.rpu == 0.0:
      if args.bamType == 'raw':
         rpu = utils.getMeanRpu(args.bamFile, args.umiTag)
         print("rpu = " + str(round(rpu,1)) + ", computed by smCounter2")
      else:
         rpu = 5.0
         print("rpu = " + str(round(rpu,1)) + ", set by smCounter2 when bamType is consensus and rpu is not given by user")
   else:
      rpu = args.rpu
      print("rpu = " + str(round(rpu,1)) + ", given by user")
      
   # set primer side
   primerSide = 'R1' if args.primerSide == 1 else 'R2'
   
   # set type of input BAM file
   bamType = 'raw' if args.bamType == 'raw' else 'consensus'
   
   # select header for normal and duplex-seq runs
   header = header_dup if args.isDuplex else header_sin
   
   #----- loop over locs
   # prepare to save to disk
   outfile_long = open('intermediate/noThres.' + args.outPrefix + '.VariantList.long.txt', 'w')
   bkgFileName = 'intermediate/bkg.' + args.outPrefix + '.txt'
   outfile_bkg = open(bkgFileName, 'w')
   outfile_long.write('\t'.join(header) + '\n')
   outfile_bkg.write('\t'.join(header_bkg) + '\n')   
   
   print('runtime' + '\t' + 'interval')
   pool = multiprocessing.Pool(args.nCPU)
   func = functools.partial(vc_wrapper, (args.bamFile, args.minBQ, args.minMQ, args.hpLen, args.mismatchThr, args.primerDist, args.consThr, rpu, primerSide, args.refGenome, args.minAltUMI, args.maxAltAllele, args.isRna, args.ds, bamType, args.umiTag, args.primerTag, args.mqTag, args.tagSeparator, args.isDuplex, args.duplexTag, args.minRpu))
   
   # process exons/intervals from bed file in parallel
   empty = True
   for interval_result in pool.map(func, locList):
      for base_result in interval_result:
         vcOutline,bkgOutline = base_result
         outfile_long.write(vcOutline)
         outfile_bkg.write(bkgOutline)
         empty = False

   # clear finished pool
   pool.close()
   pool.join()
   # close output file handles
   outfile_long.close()
   outfile_bkg.close()
   
   # calculate p-value or likelihood ratio
   thres = 'likelihood ratios ' if args.isDuplex else 'p-values '
   print("Calculating " + thres + str(datetime.datetime.now()) + "\n")
   outfile1 = 'intermediate/noThres.' + args.outPrefix + '.VariantList.long.txt'

   outfile2 = 'intermediate/' + args.outPrefix + '.VariantList.long.txt'
   outfile_lod = 'intermediate/' + args.outPrefix + '.umi_depths.lod.bedgraph'
   
   if args.isDuplex:
      pValCmd = ' '.join(['Rscript', pValCode_dup, args.dupBkgErrorDist, './', outfile1, outfile2, str(args.minAltUMI)])
   else:
      pValCmd = ' '.join(['Rscript', pValCode_sin, args.sinBkgErrorDist, './', outfile1, bkgFileName, str(seed), str(nsim), outfile2, outfile_lod, args.outPrefix, str(rpu), str(args.minAltUMI), str(args.inputVCF).lower()])

   if not empty:
      subprocess.check_call(pValCmd, shell=True)
      print("completed p-values " + str(datetime.datetime.now()) + "\n")
      # make VCFs
      vcf.makeVcf('./', outfile2, args.outPrefix, args.refGenome, args.isDuplex)
   else:
      print("empty BED or BAM, no variants detected " + str(datetime.datetime.now()) + "\n")
   
   # remove intermediate files
   #os.remove('hp.roi.bed')
   #os.remove('rep.roi.bed')
   #os.remove('sr.roi.bed')
   os.remove(outfile1)

   # log run completion
   timeEnd = datetime.datetime.now()
   print("completed running " + str(timeEnd) + "\n")
   print("total runtime: "+ str(timeEnd-timeStart) + "\n")  
   
#-----------------------------------------------------------------------------------
# pythonism to run from the command line
#-----------------------------------------------------------------------------------
if __name__ == "__main__":
   # init the argumet parser
   argParseInit()  
   # get command line arguments
   args = parser.parse_args()  
   # call main program   
   main(args)
