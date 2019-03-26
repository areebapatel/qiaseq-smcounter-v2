#!/usr/bin/Rscript
# vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
# likelihood ratio for duplex-seq variant calling
# Chang Xu, 14MAR2019

rm(list=ls())
suppressMessages(library(tidyverse))
options(stringsAsFactors = F)

args <- commandArgs(TRUE)
bkgErrorDist <- args[1]
wd <- args[2]
outlong <- args[3]
outfile <- args[4]
minAltUMI <- as.numeric(args[5])

# set working directory
setwd(wd)

# import prior distributions of error rates
load(bkgErrorDist)

# function to determine background error rates for duplex UMIs; base change specific 
duplex.p <- function(REF, ALT, sForUMT){
  if((REF == 'C' & ALT == 'A' & sForUMT >= 100) | (REF == 'G' & ALT == 'T' & sForUMT < 100)){
    p <- 9.41e-7
  } else if((REF == 'C' & ALT == 'T' & sForUMT >= 100) | (REF == 'G' & ALT == 'A' & sForUMT < 100)){
    p <- 4.28e-7
  } else if((REF == 'G' & ALT == 'A' & sForUMT >= 100) | (REF == 'C' & ALT == 'T' & sForUMT < 100)){
    p <- 5.09e-6
  } else if((REF == 'T' & ALT == 'A' & sForUMT >= 100) | (REF == 'A' & ALT == 'T' & sForUMT < 100)){
    p <- 3.59e-7
  } else if((REF == 'T' & ALT == 'C' & sForUMT >= 100) | (REF == 'A' & ALT == 'G' & sForUMT < 100)){
    p <- 1.20e-7
  } else{
    p <- 8.20e-8
  }
  return(p)
}

# function to compute likelihood under H0: no variant and H1: variant
llhr <- function(REF, ALT, sForUMT, sRevUMT, sForVMT, sRevVMT, dUMT, dVMT, a, b){
  totalN <- sForUMT + sRevUMT
  totalX <- sForVMT + sRevVMT
  
  if(totalX >= minAltUMI & ALT != 'DEL'){
    # log-likelihood of sVMT under H0, excluding terms cancelled out
    logL.sin <- lbeta(a + totalX, b + totalN - totalX) - lbeta(a, b) 
    # log-likelihood of dVMT under H0
    p.dup <- duplex.p(REF, ALT, sForUMT)
    
    logL.dup <- log(dbinom(x = dVMT, size = dUMT, prob = p.dup))
    # final log-likelihood under H0
    logL.H0 <- logL.sin + logL.dup
    
    # log-likelihood of sVMT under H1; cancelled-out terms not computed
    logL.H1 <- lbeta(totalX + dVMT + 1, totalN + dUMT - totalX - dVMT + 1) - lbeta(dVMT + 1, dUMT - dVMT + 1) - log(dUMT + 1) 
    
    # log likelihood ratio
    logLR <- min(10000, max(-10000, logL.H1 - logL.H0))

    # filter to detect descripency between sVMF and dVMF; based on Poisson test of equal means
    vmf.pval <- poisson.test(x = c(totalX, dVMT), T = c(totalN, dUMT), alternative = 'greater')$p.value
    out <- c(vmf.pval, logL.H1, logLR)

  } else{
    out <- rep(-100000, 3)
  }
  return(out)
}

# read and process smCounter output  
dat <- read.table(outlong, header=T) 
tmp <- mutate(dat, sForUMT = as.numeric(sForUMT), 
              sRevUMT = as.numeric(sRevUMT), 
              sForVMT = as.numeric(sForVMT), 
              sRevVMT = as.numeric(sRevVMT),
              dUMT = as.numeric(dUMT), 
              dVMT = as.numeric(dVMT)) %>%
  select(REF, ALT, sForUMT, sRevUMT, sForVMT, sRevVMT, dUMT, dVMT)

# calculate log-likelihood ratio of each position
lr <- mapply(llhr, tmp$REF, tmp$ALT, tmp$sForUMT, tmp$sRevUMT, tmp$sForVMT, tmp$sRevVMT, tmp$dUMT, tmp$dVMT, a = bkg.ga$a, b = bkg.ga$b)

# save to disk
dat <- mutate(dat, plowDupVMF = lr[1,], logLH1 = lr[2,], logLR = lr[3,], 
              FILTER = ifelse(plowDupVMF < 0.05 & plowDupVMF >= 0, ifelse(FILTER=='PASS', 'LowDupVMF', paste0(FILTER, ';LowDupVMF')), FILTER)) %>%
  select(CHROM, POS, REF, ALT, TYPE, sUMT, sVMT, sVMF, dUMT, dVMT, dVMF, DP, VDP, VAF, RefForPrimer, RefRevPrimer, primerOR, pLowQ, hqUmiEff, allUmiEff, refMeanRpb, altMeanRpb, rpbEffectSize, repType, hpInfo, simpleRepeatInfo, tandemRepeatInfo, FR, MT, UFR, sUMT_A, sUMT_T, sUMT_G, sUMT_C, dUMT_A, dUMT_T, dUMT_G, dUMT_C, plowDupVMF, logLH1, logLR, FILTER) %>%
  write_delim(outfile, delim='\t', col_names=T)

