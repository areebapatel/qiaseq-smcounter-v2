#!/usr/bin/Rscript
# vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
# p-values for duplex-seq variant calling based on intersection-union test
# Aggressive impute (smaller prior sigma) when dVMF > 2 * sVMF
# Chang Xu, 30MAY2019

rm(list=ls())
suppressMessages(library(tidyverse))
suppressMessages(library(rmutil))
suppressMessages(library(poibin))

options(stringsAsFactors = F)

args <- commandArgs(TRUE)
bkgErrorDist <- args[1]
wd <- args[2]
outlong <- args[3]
outfile <- args[4]
min.dUMT <- as.numeric(args[5])

# set working directory
setwd(wd)

# import prior distributions of error rates
load(bkgErrorDist)

# function to get complementary base
complementary <- function(x){
  y <- case_when(x == 'a' ~ 't', 
                 x == 'c' ~ 'g', 
                 x == 'g' ~ 'c', 
                 x == 't' ~ 'c', 
                 x == 'A' ~ 'T', 
                 x == 'C' ~ 'G', 
                 x == 'G' ~ 'C', 
                 x == 'T' ~ 'A', 
                 TRUE ~ 'n')
  return(y)
}

# function to estimate Beta parameters based on mean and standard deviation
beta.ab <- function(mu, sigma) {
  a <- mu * (mu * (1 - mu) / sigma^2 - 1)
  b <- (1 - mu) * (mu * (1 - mu) / sigma^2 - 1)
  return(c(a,b))
}

# function to impute when dVMF > sVMF
impute <- function(dVMF, dUMT, sVMT, sUMT){
  # prior sd
  dVMF.raw <- dVMF / 100
  sigma.prior <- sqrt(dVMF.raw * (1 - dVMF.raw) / dUMT)
  
  # get beta distribution (prior)
  ab <- beta.ab(dVMF.raw, sigma.prior)
  a <- ab[1]
  b <- ab[2]
  
  # posterior mean of Beta(x + a, n - x  + b)
  sVMF.post_mean <- (sVMT + a) / (sUMT + a + b)
  
  # imputed sVMT 
  sVMT.imp <- max(sVMT, round(sUMT * sVMF.post_mean))
  return(sVMT.imp)
}

# function to calculate -log10(pval) for dupelx UMI 
mlog10pval.dup <- function(x, n, p){
  pval <- ifelse(x == 0, 1, pbinom(x - 1, n, p, lower.tail = FALSE))
  return(min(200, -log10(pval)))
}

# function to calculate -log10(pval) for single UMI 
mlog10pval.sin <- function(x, n, a, b){
  # convert to mean-dispersion parameterization of Beta
  s <- a + b
  m <- a / (a + b)
  
  # p-value based on beta-binomial distribution
  pval <- ifelse(x == 0, 1, max(0, 1 - pbetabinom(x-1, n, m, s)))
  return(min(200, -log10(pval)))
}

# function to compute the final p-values from data
mlog10pval.final <- function(REF, ALT, TYPE, sUMT, sVMT, sVMF, dForUMT, dForVMT, dRevUMT, dRevVMT, dVMF, dUMT){
  # do not compute p-values at positions within a deletion
  if(ALT == 'DEL'){
    res <- rep(NA, 4)
    return(res)
  }
  
  # assume indels have G>A error profile
  if(TYPE == 'INDEL'){
    REF <- 'G'
    ALT <- 'A'
  }
  
  ### single UMI - face value of REF and ALT - transition: G>A. transversion: C>T
  dist.sin <- eval(parse(text = paste0('error.profile.duplex$', tolower(REF), tolower(ALT))))
  a.sin <- dist.sin[1]
  b.sin <- dist.sin[2]
  
  # impute sVMT 
  imp.cond <- !is.na(dUMT) & !is.na(dVMF) & !is.na(sVMF) & dVMF > 2 * sVMF & dUMT > 0 & sVMF < 0.1
  sVMT.imp <- ifelse(imp.cond, impute(dVMF, dUMT, sVMT, sUMT), sVMT)
  
  # -log10(pval)
  logpval.sin <- mlog10pval.sin(sVMT.imp, sUMT, a.sin, b.sin)
  
  ### duplex UMI - error distribution depends on primer strand 
  dist.dup.for <- eval(parse(text = paste0('error.profile.duplex$', tolower(REF), tolower(ALT))))
  p.dup.for <- dist.dup.for[3]
  dist.dup.rev <- eval(parse(text = paste0('error.profile.duplex$', tolower(complementary(REF)), tolower(complementary(ALT)))))
  p.dup.rev <- dist.dup.rev[3]
  
  if(dForUMT >= min.dUMT & dRevUMT >= min.dUMT){
    # weight of probabilities
    wt <- c(dForUMT + dForVMT, dRevUMT + dRevVMT)
    ps <- c(p.dup.for, p.dup.rev)
    x <- dForVMT + dRevVMT - 1
    cdf <- ifelse(x >= 0, ppoibin(x, ps, method = 'DFT-CF', wts = wt), 0)
    logpval.dup <- min(200, -log10(1 - cdf))
  } else if(dForUMT >= min.dUMT) {
    logpval.dup <- mlog10pval.dup(dForVMT, dForUMT + dForVMT, p.dup.for)
  } else if(dRevUMT >= min.dUMT) {
    logpval.dup <- mlog10pval.dup(dRevVMT, dRevUMT + dRevVMT, p.dup.rev)
  } else{
    logpval.dup <- NA
  } 
  
  # final p
  logpval <- min(logpval.sin, logpval.dup)
  
  # output
  res <- c(logpval.sin, logpval.dup, logpval, sVMT.imp)
  return(res)
}

# read and process smCounter output  
dat <- read.table(outlong, header = T) 
tmp <- mutate(dat, sUMT = as.numeric(sUMT), 
              sVMT = as.numeric(sVMT), 
              sVMF = as.numeric(sVMF), 
              dForUMT = as.numeric(dForUMT), 
              dRevUMT = as.numeric(dRevUMT), 
              dForVMT = as.numeric(dForVMT), 
              dRevVMT = as.numeric(dRevVMT), 
              dVMF = as.numeric(dVMF), 
              dUMT = as.numeric(dUMT)) %>%
  select(REF, ALT, TYPE, sUMT, sVMT, sVMF, dForUMT, dForVMT, dRevUMT, dRevVMT, dVMF, dUMT)

# calculate p-values at each position
pvals <- mapply(mlog10pval.final, tmp$REF, tmp$ALT, tmp$TYPE, tmp$sUMT, tmp$sVMT, tmp$sVMF, tmp$dForUMT, tmp$dForVMT, tmp$dRevUMT, tmp$dRevVMT, tmp$dVMF, tmp$dUMT)

# save to disk
dat <- mutate(dat, logPvalSin = pvals[1,], 
              logPvalDup = pvals[2,], 
              logPval = pvals[3,], 
              # temporary hack to assign imputed sVMT to dUMT_C to keep consistent header
              dUMT_C = pvals[4,]) %>%    
  select(CHROM, POS, REF, ALT, TYPE, sUMT, sVMT, sVMF, dUMT, dVMT, dVMF, DP, VDP, VAF, sForUMT, sForVMT, sRevUMT, sRevVMT, dForUMT, dForVMT, dRevUMT, dRevVMT, primerOR, pLowQ, hqUmiEff, allUmiEff, refMeanRpb, altMeanRpb, rpbEffectSize, repType, hpInfo, simpleRepeatInfo, tandemRepeatInfo, FR, MT, UFR, sUMT_A, sUMT_T, sUMT_G, sUMT_C, dUMT_A, dUMT_T, dUMT_G, dUMT_C, logPvalSin, logPvalDup, logPval, FILTER) %>%
  write_delim(outfile, delim='\t', col_names=T)


