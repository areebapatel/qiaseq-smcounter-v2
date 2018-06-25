#!/usr/bin/python
# vim: tabstop=9 expandtab shiftwidth=3 softtabstop=3
# get homopolymer regions
# Chang Xu. 23May2016
# bug fix - plus 1 on start; 23JAN2018

import os
import sys
import pysam
from collections import defaultdict

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
               outline = '\t'.join([chrom, str(max(realL, start)), str(min(realR, end)), 'HP', str(hpLen), '1', str(hpLen), origRef[hpL]]) + '\n'
               outfile.write(outline)
            hpL = i 


   outfile.close()


#----------------------------------------------------------------------------------------------
#pythonism to run from the command line
#----------------------------------------------------------------------------------------------
if __name__ == "__main__":
   bedName = sys.argv[1]
   outName = sys.argv[2]
   minLength = int(sys.argv[3])
   refg = sys.argv[4]
   seqType = sys.argv[5]
   findhp(bedName, outName, minLength,refg,seqType)


