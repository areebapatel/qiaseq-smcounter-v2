#!/usr/bin/python
# vim: tabstop=9 expandtab shiftwidth=3 softtabstop=3

import os
import subprocess
import multiprocessing
import datetime
import argparse
import traceback

import utils
import vc

# global constants
codePath = os.path.dirname(os.path.abspath(__file__))
pValCode = os.path.join(codePath,'get_pvalue_v2.R')
locChunkLen = 1000
maxDnaReadDepth = 1000000000
downsamplePileupStackThr = 10 ** 5
parser = None
header_1 = ['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'sUMT', 'sForUMT', 'sRevUMT', 'sVMT', 'sForVMT', 'sRevVMT', 'sVMF', 'sForVMF', 'sRevVMF', 'VDP', 'VAF', 'RefForPrimer', 'RefRevPrimer', 'primerOR', 'pLowQ', 'hqUmiEff', 'allUmiEff', 'refMeanRpb', 'altMeanRpb', 'rpbEffectSize', 'repType', 'hpInfo', 'simpleRepeatInfo', 'tandemRepeatInfo', 'DP', 'FR', 'MT', 'UFR', 'sUMT_A', 'sUMT_T', 'sUMT_G', 'sUMT_C', 'FILTER']
header_2 = ['CHROM', 'POS', 'REF', 'A/G', 'G/A', 'C/T', 'T/C', 'A/C', 'C/A', 'A/T', 'T/A', 'C/G', 'G/C', 'G/T', 'T/G', 'negStrand', 'posStrand', 'AllSMT' ]

#----------------------------------------------------------------------------------------------
# global for argument parsing (hack that works when calling from either command line or pipeline)
# this is done inside a function because multiprocessing module imports the script
#------------------------------------------------------------------------------------------------
def argParseInit():  
   global parser
   parser = argparse.ArgumentParser(description='smCounter2: variant calling using Unique Molecular Identifiers')
   parser.add_argument('--runPath', default=None, help='path to working directory')
   parser.add_argument('--bedTarget', default=None, help='BED file')
   parser.add_argument('--bamFile', default=None, help='BAM file')
   parser.add_argument('--outPrefix', default=None, help='file name prefix')
   parser.add_argument('--nCPU', type=int, default=1, help='number of CPU to use in parallel')
   parser.add_argument('--minBQ', type=int, default=25, help='minimum base quality allowed for analysis')
   parser.add_argument('--minMQ', type=int, default=50, help="minimum mapping quality allowed for analysis. If the bam is tagged with its mate's mapq, then the minimum of the R1 and R2 mapq will be used for comparison, if not each read is compared independently.")
   parser.add_argument('--hpLen', type=int, default=10, help='minimum length for homopolymers')
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
   parser.add_argument('--ds', type=int, default=10000, help='down sample if number of UMIs greater than this value (for RNA only)')
   parser.add_argument('--bamType', type=str, default='raw', help='raw (default): raw BAM file with UMIs; consensus: consensused BAM file')
   parser.add_argument('--inputVCF', type=str, default='none', help='optional input VCF file;')
 
#------------------------------------------------------------------------------------------------
# wrapper function for "vc()" - because Python multiprocessing module does not pass stack trace; from runone/smcounter.py by John Dicarlo
#------------------------------------------------------------------------------------------------
def vc_wrapper(*args):
   try:
      output = vc.vc(*args)
   except Exception as e:
      print("Exception thrown in vc() function at genome location:", args[1], args[2])
      output = ("Exception thrown!\n" + traceback.format_exc(),'no_bg')
      print output[0]
      raise Exception(e)
   return output

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

   # convert VCF to BED if inputVCF is not 'none'
   bedTarget = args.bedTarget if args.inputVCF == 'none' else utils.vcf2bed(args.inputVCF)
   
   # gather repetitive regions information
   hpRegion = utils.getHpInfo(bedTarget, args.refGenome, args.isRna, args.hpLen)
   repRegion = utils.getTrInfo(bedTarget, args.repBed, args.isRna, args.hpLen)
   srRegion = utils.getOtherRepInfo(bedTarget, args.srBed, args.isRna, args.hpLen)

   # read in bed file and create a list of positions, annotated with repetitive region
   locList = utils.getLocList(bedTarget, hpRegion, repRegion, srRegion)

   # calculate rpb if args.rpb = 0
   if args.rpb == 0.0:
      if args.bamType == 'raw':
         rpb = utils.getMeanRpb(args.bamFile) 
         print("rpb = " + str(round(rpb,1)) + ", computed by smCounter2")
      else:
         rpb = 5.0
         print("rpb = " + str(round(rpb,1)) + ", set by smCounter2 when bamType is consensus and rpb is not given by user")
   else:
      rpb = args.rpb
      print("rpb = " + str(round(rpb,1)) + ", given by user")
      
   # set primer side
   primerSide = 'R1' if args.primerSide == 1 else 'R2'

   # set type of input BAM file
   bamType = 'raw' if args.bamType == 'raw' else 'consensus'

   #----- loop over locs
   # prepare to save to disk
   outfile_long = open('intermediate/nopval.' + args.outPrefix + '.VariantList.long.txt', 'w')
   bkgFileName = 'intermediate/bkg.' + args.outPrefix + '.txt'
   outfile_bkg = open(bkgFileName, 'w')

   outfile_long.write('\t'.join(header_1) + '\n')
   outfile_bkg.write('\t'.join(header_2) + '\n')

   # process in chunks (for more granular logging)
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
      results = [pool.apply_async(vc_wrapper, args=(args.bamFile, x[0], x[1], x[2], x[3], x[4], x[5], args.minBQ, args.minMQ, args.hpLen, args.mismatchThr, args.primerDist, args.mtThreshold, rpb, primerSide, args.refGenome, args.minAltUMI, args.maxAltAllele, args.isRna, args.ds, args.bamType)) for x in locChunk]
      
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
   pValCmd = ' '.join(['Rscript', pValCode, args.runPath, outfile1, bkgFileName, str(seed), str(nsim), outfile2, outfile_lod, args.outPrefix, str(rpb), str(args.minAltUMI), args.inputVCF])
   subprocess.check_call(pValCmd, shell=True)
   print("completed p-values at " + str(datetime.datetime.now()) + "\n")

   ## make VCFs
   utils.makeVcf(args.runPath, outfile2, args.outPrefix)

   # remove intermediate files
   os.remove('hp.roi.bed')
   os.remove('rep.roi.bed')
   os.remove('sr.roi.bed')
   os.remove(outfile1)

   # log run completion
   timeEnd = datetime.datetime.now()
   print("completed running at " + str(timeEnd) + "\n")
   print("total time: "+ str(timeEnd-timeStart) + "\n")  
   
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
