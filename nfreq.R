library(Rsamtools)
library(VariantAnnotation)
library(SomaticSignatures)
library(BSgenome.Hsapiens.UCSC.hg19)

###############################################################################
## Functions
pileupFreq <- function(pileupres) {
  nucleotides <- levels(pileupres$nucleotide)
  res <- split(pileupres, pileupres$seqnames)
  res <- lapply(res, function (x) {split(x, x$pos)})
  res <- lapply(res, function (positionsplit) {
    nuctab <- lapply(positionsplit, function(each) {
      chr = as.character(unique(each$seqnames))
      pos = as.character(unique(each$pos))
      tablecounts <- sapply(nucleotides, function (n) {sum(each$count[each$nucleotide == n])})
      c(chr,pos, tablecounts)
    })
    nuctab <- data.frame(do.call("rbind", nuctab),stringsAsFactors=F)
    rownames(nuctab) <- NULL
    nuctab
  })
  res <- data.frame(do.call("rbind", res),stringsAsFactors=F)
  rownames(res) <- NULL
  colnames(res) <- c("seqnames","start",levels(pileupres$nucleotide))
  res[3:ncol(res)] <- apply(res[3:ncol(res)], 2, as.numeric)
  res
}

filtVcfFfpe <- function(vcf) {
  keep <- !((ref(vcf)=='C' & unlist(rowRanges(vcf)$ALT)=='T' & unlist(info(vcf)$AF)>=0.5) |
              (ref(vcf)=='G' & unlist(rowRanges(vcf)$ALT)=='A' & unlist(info(vcf)$AF)>=0.5))
  vcf <- vcf[keep,]
  return(vcf)
}

checkCoverage <- function(bamfile, chr, pos) {
  bf <- BamFile(bamfile)
  param <- ScanBamParam(which=GRanges(chr, IRanges(start=pos, end=pos)))
  p_param <- PileupParam(max_depth=1000, 
                         ignore_query_Ns=FALSE,
                         min_nucleotide_depth=0,
                         distinguish_strand=FALSE,
                         min_mapq=0,
                         min_base_quality=0,
                         min_minor_allele_depth=0
                         )
  res <- pileup(bf, scanBamParam=param, pileupParam=p_param)
  if (dim(res)[1]==0) {
    res <- rbind(res, data.frame(seqnames=chr,
                                 pos=pos,
                                 nucleotide=NA,
                                 count=0,
                                 which_label=paste(chr,":",pos,"-",pos, sep='')
                                 )
                 )
    return(res)
  } else {
    return(res)
  }
}

pertNeg <- function(vcf, bam) {
  cov<- NA
  for(i in 1:length(vcf)) {
    res <- checkCoverage(bam, 
                         as.vector(seqnames(vcf[i])),
                         ranges(vcf[i])@start
    )
    cov <- rbind(cov,res)
  }
  cov <- cov[-1,]
  pres <-pileupFreq(cov)  
  return(pres)
}

#compute b allele frequency for given position
ballelfreq <- function(bam, chr, pos) {
  ref <- as.vector(as.vector(getSeq(Hsapiens, chr, pos, pos)))
  res <- pileupFreq(checkCoverage(bam, chr, pos))
  covalt <- (res$A + res$C + res$G + res$T) - res[[ref]]
  baf <- covalt / (covalt + res[[ref]])
  return(baf)
}

###############################################################################

# get the bam files
bamfile1 <- "/home/tobias/AWS/storage/t47d_rna/variants/T47D-100ng_S7/merged.conv.sort.rd.split.realigned.recal.bam"
bamfile2 <- "/home/tobias/AWS/storage/t47d_rna/variants/T47D-400ng_S8/merged.conv.sort.rd.split.realigned.recal.bam"
bf1 <- BamFile(bamfile1)
bf2 <- BamFile(bamfile2)

# get the vcf files
vcf1 <- readVcf("/home/tobias/AWS/storage/t47d_rna/variants/T47D-100ng_S7/T47D-100ng_S7_final_variants_sort.vcf.gz", "hg19")
vcf2 <- readVcf("/home/tobias/AWS/storage/t47d_rna/variants/T47D-400ng_S8/T47D-400ng_S8_final_variants_sort.vcf.gz", "hg19")

# only SNVs
vcf1.filt <- vcf1[which(ranges(vcf1)@width==1), ]
vcf2.filt <- vcf2[which(ranges(vcf2)@width==1), ]

# FFPE hard filter, >=50% AF fo C>T / G>A
vcf1.filt <- filtVcfFfpe(vcf1.filt)
vcf2.filt <- filtVcfFfpe(vcf2.filt)

diff1v2 <- setdiff(rowRanges(vcf1.filt), rowRanges(vcf2.filt))
diff2v1 <- setdiff(rowRanges(vcf2.filt), rowRanges(vcf1.filt))
common <- intersect(rowRanges(vcf1.filt), rowRanges(vcf2.filt))

# variants present in replicate 1, but not in replicare 2 ...
x1 <- pertNeg(diff1v2, bamfile2)
depthSample2 <- rowSums(x1[,c(3,4,5,6)])
table(is.na(depthSample2) | depthSample2==1 | depthSample2==2) / length(depthSample2) * 100

# variants present in replicate 2, but not in replicare 1 ...
x2 <- pertNeg(diff2v1, bamfile1)
depthSample1 <- rowSums(x2[,c(3,4,5,6)])
table(is.na(depthSample1) | depthSample1==1 | depthSample1==2) / length(depthSample1) * 100

#################################
# compute b allel freq for vcf1
tab <- data.frame(pos=as.numeric(start(vcf1)), 
                  chr=as.vector(seqnames(vcf1)))

bfrq1 <- apply(tab, 1, function(x) {
  ballelfreq(bamfile1, x['chr'], as.numeric(x['pos']))
})
  
