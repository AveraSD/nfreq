library(Rsamtools)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)
library(myvariant)

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

getCosmic <- function(n=1000){
  cosmic <- list()
  i <- 0
  while (i < n){
    q <- queryVariant(q="_exists_:cosmic AND _exists_:cadd", fields="vcf", skip=i, size=i+1000)
    i <- i + 1000}
  
  cosmic <- c(q, cosmic)
  return(cosmic)
}

# bed to granges
# function from: http://davetang.org/muse/2015/02/04/bed-granges/
bed_to_granges <- function(file){
  df <- read.table(file,
                   header=F,
                   stringsAsFactors=F)
  
  if(length(df) > 6){
    df <- df[,-c(7:length(df))]
  }
  
  if(length(df)<3){
    stop("File has less than 3 columns")
  }
  
  header <- c('chr','start','end','id','score','strand')
  names(df) <- header[1:length(names(df))]
  
  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }
  
  library("GenomicRanges")
  
  if(length(df)==3){
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
  } else if (length(df)==5){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
  } else if (length(df)==6){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
  }
  return(gr)
}