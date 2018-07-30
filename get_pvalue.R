#!/usr/bin/Rscript
# vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
# integrate prior info of bkg error distribution (based on duplex data) and data specific error distribution
# modified the prior distribution using more stringint consensus criteria. 
# different prior distributions for including/excluding 1 read pair UMIs
# allow Beta distribution adjustment when there is no VCF input
# Chang Xu, 10July2018

rm(list=ls())
library(plyr)
bkgErrorDistSimulation <- '/srv/qgen/data/annotation/bkg.error.v2.RData'

##############################
##       Parameters         ##
##############################
args <- commandArgs(TRUE)
wd <- args[1]
outlong <- args[2]
bkgfile <- args[3]
seed <- as.numeric(args[4])
nsim <- as.numeric(args[5])
outfile_pval <- args[6]
outfile_bedgraph <- args[7]
outprefix <- args[8]
rpb <- as.numeric(args[9])
minAltUMI <- as.numeric(args[10])
inputVCF <- args[11]
min.mtDepth <- 1000

# set working directory
setwd(wd)
set.seed(seed)

##############################################
########          Function            ########
########         Definitions          ########
##############################################
# function to calculate standard deviation
beta.sd <- function(a,b) sqrt(a*b) / ((a+b) * sqrt(a+b+1))
# function to estimate a
calc.a <- function(mu, sigma) mu * (mu*(1-mu) / sigma^2 - 1)
# function to estimate b
calc.b <- function(mu, sigma) (1-mu) * (mu*(1-mu) / sigma^2 - 1)
# function to compute p values
calc.pval <- function(TYPE, REF, ALT, sForUMT, sRevUMT, sForVMT, sRevVMT, p.high.final, p.low.final){
  totalN <- sForUMT + sRevUMT
  totalX <- sForVMT + sRevVMT
  if(totalX >= minAltUMI){    
    if(TYPE=='INDEL' | (REF=='A' & ALT=='G') | (REF=='G' & ALT=='A') | (REF=='C' & ALT=='T') | (REF=='T' & ALT=='C')){
      pr <- p.high.final
    } else{
      pr <- p.low.final
    }
    tmp <- pbinom(q=totalX-1, size=totalN, prob=pr, lower.tail=F)
    pval <- ifelse(totalN==0, 1, mean(tmp, na.rm=T))
  } else{
    pval <- 1.0 
  }
  return(pval)
}
# function to find p-value
pval <- function(n, x, p){
#       @param int    n    :   The UMI depth at a particular site
#       @param float  x    :   Number of variant UMIs at that site
#       @param vector p    :   Vector of values simulated from the
#                              background error distribution of transitions 

  if(x >= 3){
    tmp <- pbinom(q=x-1, size=n, prob=p, lower.tail=F)
    pval <- ifelse(n==0, 1, mean(tmp, na.rm=T))
  } else{
    pval <- NA
  }
  return(pval)
}
# function to find the LOD
calc_lod <- function(n,p.high){
#      @param   int n        : The UMI depth to calculate the lod for
#      @param   float p.high : Vector of values simulated from the
#                              background error distribution of transitions

  # high lod
  low <- 3
  up <- n
  x.high <- max(3, round(0.005 * n))
  while(up - low > 1){  
    p <- pval(n, x.high, p.high)
    if(p >= 1e-6){
      low <- x.high
      x.high <- ceiling(mean(c(x.high, up)))
    } else{
      up <- x.high
      x.high <- floor(mean(c(x.high, low)))
    }
  }
  lod.high <- x.high / n

  if(is.na(lod.high)){
    print(n)
    print(p.high)
    print(lod.high)
  }
  return(lod.high)
}
# function to collapse same value(lod/coverage) columns
# and write a bedgraph file
output_bedgraph <- function(df,outfile,header,val_col="foo"){
#                  @param dataframe df      : The input dataframe to iterate over
#                  @param string    val_col : The column name of the value in the bedgraph
#                  @param string    outfile : The output file path
#                  @param string    header  : The header for the bedgraph file

   file_handle <- file(outfile,"w")
   cat(header,file=file_handle)
   prev_val <- NULL
   test <- c("lod","sumt")
   for (row in 1:nrow(df)) {
      val <- df[row,val_col]
      chr <- df[row, "chr"]
      pos <- df[row, "pos"]
      if (is.null(prev_val)) {
	 prev_val <- val
	 prev_chr <- chr
	 prev_pos <- pos
	 init_pos <- pos
	 next # skip first iteration of loop
      }
      else {
	 if (prev_chr != chr) {
            out <- paste(prev_chr,"\t",as.integer(init_pos-1),"\t",as.integer(prev_pos),"\t",round(prev_val,5),"\n",sep="")
            cat(out,file=file_handle)
	    init_pos <- pos
	 }
	 else if (prev_val != val) {	    
            out <- paste(prev_chr,"\t",as.integer(init_pos-1),"\t",as.integer(prev_pos),"\t",round(prev_val,5),"\n",sep="")
            cat(out,file=file_handle)
	    init_pos <- pos
	 }
	 prev_val <- val
	 prev_chr <- chr
	 prev_pos <- pos
      }
   }
   # finish out last line of the file
   out <- paste(prev_chr,"\t",as.integer(init_pos-1),"\t",as.integer(prev_pos),"\t",round(prev_val,3),"\n",sep="")
   cat(out,file=file_handle)
   close(file_handle)
}
################ END OF FUNCTIONS #####################


#########################################################
########       Begin imperative computations     ########
########                  And                    ########
########         Writing out output files        ########
#########################################################
# define constants
cols <- c('chrom', 'pos', 'ref', 'AG', 'GA', 'CT', 'TC', 'AC', 'AT', 'CA', 'CG', 'GC', 'GT', 'TA', 'TG', 'neg.strand', 'pos.strand', 'all.smt')
out <- NULL

# read in smCounter output 
dat <- read.delim(outlong, header=T, stringsAsFactors=F)

# read in prior information
load(bkgErrorDistSimulation)
if(rpb >= 3.0){
  top4 <- bkg.error$top4.exclude.1rpUMI
} else{
  top4 <- bkg.error$top4.include.1rpUMI
}

a.ga.orig <- top4$shape1[2]
b.ga.orig <- top4$shape2[2]
sigma.high <- beta.sd(a.ga.orig, b.ga.orig)

a.ct.orig <- top4$shape1[3]
b.ct.orig <- top4$shape2[3]

# proportion of zeros
p0.high <- top4$p0[2]
p0.low <- top4$p0[3]
n0.high <- floor(nsim * p0.high)
n0.low <- floor(nsim * p0.low)

# read in data-specific background error file - only when input is not VCF
if(inputVCF == 'none'){
  bkg <- read.delim(bkgfile, header=T, stringsAsFactors=F, sep='\t')
  colnames(bkg) <- cols

  ################### bkg errors from the readset ##################
  # A/G error rate from data 
  tmp <- bkg[(bkg$ref=='A' & bkg$neg.strand > min.mtDepth) | (bkg$ref=='T' & bkg$pos.strand > min.mtDepth), ]
  tmp$all <- ifelse(tmp$ref=='A', tmp$neg.strand, tmp$pos.strand)
  d.ag <- tmp$AG / tmp$all
  tmp <- tmp[d.ag < 0.01,]
  mean.ag <- sum(tmp$AG) / sum(tmp$all)
  n.ag <- nrow(tmp)

  # G/A error rate from data
  tmp <- bkg[(bkg$ref=='G' & bkg$neg.strand > min.mtDepth) | (bkg$ref=='C' & bkg$pos.strand > min.mtDepth), ]
  tmp$all <- ifelse(tmp$ref=='G', tmp$neg.strand, tmp$pos.strand)
  d.ga <- tmp$GA / tmp$all
  tmp <- tmp[d.ga < 0.01,]
  mean.ga <- sum(tmp$GA) / sum(tmp$all)
  n.ga <- nrow(tmp)

  # highest error rate
  mu.high <- max(mean.ag, mean.ga)
  n.high <- min(n.ag, n.ga)

  if(is.na(mu.high) | is.na(n.high) | n.high < 100) {
    p.high <- rbeta(n=nsim, shape1=a.ga.orig, shape2=b.ga.orig)
  } else{
    a.high <- calc.a(mu.high, sigma.high)
    b.high <- calc.b(mu.high, sigma.high)
    p.high <- rbeta(n=nsim, shape1=a.high, shape2=b.high)
  }
} else{  # when input is VCF, use the original beta parameters
  p.high <- rbeta(n=nsim, shape1=a.ga.orig, shape2=b.ga.orig)
}

# low error rates always determined by the prior only
p.low <- c(rbeta(n=nsim-n0.low, shape1=a.ct.orig, shape2=b.ct.orig), rep(0, n0.low))

# compute limit of detection (lod)  for binned sUMT values 
# this is the lowest allele fraction variant which can be called for a given UMI depth at a site
bin_width = 10
all_sUMT_bin_vals <- seq(from = min(dat$sUMT), to = min(10000,max(dat$sUMT)), by = bin_width)
all_sUMT_bins <- seq(from=1,to=length(all_sUMT_bin_vals),by=1)
binned_lod_vals <- sapply(all_sUMT_bin_vals, calc_lod, p.high=p.high)
max_bin <- length(all_sUMT_bin_vals)

get_bin_indices <- function(sumt,max_bin){
  if(sumt > 10000) {
    return (max_bin)
  }
  else {
    return (floor((sumt - min(dat$sUMT) + bin_width)/bin_width))
  }
}
lod_for_sUMT <- binned_lod_vals[sapply(dat$sUMT,get_bin_indices,max_bin=max_bin)]
# write lod bedgraph file
lod_df <- data.frame(chr=dat$CHROM,pos=dat$POS,lod=lod_for_sUMT)
header <- sprintf("track type=bedGraph name='%s.variant-calling-lod'\n",outprefix)
outfile <- sprintf("%s.umi_depths.variant-calling-lod.bedgraph",outprefix)
output_bedgraph(lod_df,outfile,header,"lod")
lod.quantiles <- quantile(lod_df$lod,probs=c(0.01,0.05,0.10,0.50,0.90,0.95,0.99))
write.table(lod.quantiles, paste(outfile,".quantiles.txt",sep=""), sep='|', row.names=T, col.names=F, quote=F)
# write sUMT bedgraph file
sumt_df <- data.frame(chr=dat$CHROM,pos=dat$POS,sumt=dat$sUMT)
header <- sprintf("track type=bedGraph name='%s.umi_depths.variant-calling-input'\n",outprefix)
outfile <- sprintf("%s.umi_depths.variant-calling-input.bedgraph",outprefix)
output_bedgraph(sumt_df,outfile,header,"sumt")

# compute p-values
dat$sForUMT <- as.numeric(dat$sForUMT)
dat$sRevUMT <- as.numeric(dat$sRevUMT)
dat$sForVMT <- as.numeric(dat$sForVMT)
dat$sRevVMT <- as.numeric(dat$sRevVMT)
tmp <- subset(dat, select=c(TYPE, REF, ALT, sForUMT, sForVMT, sRevUMT, sRevVMT))
pval <- mdply(tmp, calc.pval, p.high.final=p.high, p.low.final=p.low)

# set mininum at 1e-200 to avoid log(0)
raw.pval <- pmax(1e-200, pval$V1)
# take -log10
dat$logpval <- round(-log10(raw.pval), 2)

# save to disk
dat <- subset(dat, select=c(CHROM, POS, REF, ALT, TYPE, sUMT, sForUMT, sRevUMT, sVMT, sForVMT, sRevVMT, sVMF, sForVMF, sRevVMF, VDP, VAF, RefForPrimer, RefRevPrimer, primerOR, pLowQ, hqUmiEff, allUmiEff, refMeanRpb, altMeanRpb, rpbEffectSize, repType, hpInfo, simpleRepeatInfo, tandemRepeatInfo, DP, FR, MT, UFR, sUMT_A, sUMT_T, sUMT_G, sUMT_C, logpval, FILTER))
write.table(dat, outfile_pval, sep='\t', row.names=F, col.names=T, quote=F)


