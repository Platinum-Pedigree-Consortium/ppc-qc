#' Function to calculate region ploidy based on the alignments of phased assemblies to the reference genome.
#' It takes as input contig alignment position for haplotype 1 and 2 as in input and reports a set of regions 
#' that has exactly 2n, 1n or multi-copy alignments.
#'
#' @param h1.gr A \code{\link{GRanges-class}} object of genomic regions where haplotype 1 contigs align to the reference genome.
#' @param h1.gr A \code{\link{GRanges-class}} object of genomic regions where haplotype 2 contigs align to the reference genome.
#' @param sample.sex A sex of sample being process, either male or female.
#' @param reference.fai A reference FASTA index to be used to obtain chromosome lengths.
#' @return A \code{list} 
#' @author David Porubsky
#'
getPloidy <- function(h1.gr, h2.gr, sample.sex=NULL, reference.fai=NULL) {
  ## Check user input
  if (!is(h1.gr, 'GRanges')) {
    stop("h1.gr has to be a valid 'GRanges' object !!!")
  }
  if (!is(h2.gr, 'GRanges')) {
    stop("h2.gr has to be a valid 'GRanges' object !!!")
  }
  
  ## Get reference sequence lengths
  if (!is.null(reference.fai)) {
    ref.df <- read.table(reference.fai)
    ref.df <- ref.df[,c(1:2)]
    colnames(ref.df) <- c('seqnames', 'seq.len')
    seq.len <- ref.df$seq.len
    names(seq.len) <- ref.df$seqnames
  }
  
  ## Make sure all seqlengths are defined
  if (any(c(is.na(seqlengths(h1.gr)), is.na(seqlengths(h2.gr))))) {
    if (!is.null(reference.fai)) {
      seqlengths(h1.gr) <- seq.len[seqlevels(h1.gr)]
      seqlengths(h2.gr) <- seq.len[seqlevels(h2.gr)]
    } else {
      stop("Sequence lengths are not defined !!!")
    } 
  }
  
  ## Get coverage per haplotype
  gr.cov.h1 <- as( coverage(h1.gr), 'GRanges')
  gr.cov.h2 <- as( coverage(h2.gr), 'GRanges')
  
  ## Get diploid coverage
  gr.cov.dip <- as( coverage(gr), 'GRanges')
  seqlengths(gr.cov.dip) <- seq.len[seqlevels(gr.cov.dip)]
  ## Get diploid regions (exactly haploid in both haplotypes)
  dip.gr <- as( coverage(c(gr.cov.h1[gr.cov.h1$score == 1], gr.cov.h2[gr.cov.h2$score == 1])), 'GRanges' )
  dip.gr <- dip.gr[dip.gr$score == 2]
  ## Get more then haploid regions
  multi.gr <- as( coverage(c(gr.cov.h1[gr.cov.h1$score > 1], gr.cov.h2[gr.cov.h2$score > 1])), 'GRanges' )
  multi.gr <-  multi.gr[multi.gr$score > 0]
  ## Get haploid regions (exactly haploid in both haplotypes)
  hap.gr <- gr.cov.dip[gr.cov.dip$score == 1]
  if (sample.sex == 'male') {
    hap.aut.gr <- hap.gr[!seqnames(hap.gr) %in% c('chrX', 'chrY')]
    hap.sex.gr <- hap.gr[seqnames(hap.gr) %in% c('chrX', 'chrY')]
    cov.hap <- sum(width(hap.aut.gr))
    cov.sex.hap <- sum(width(hap.sex.gr))
    hap.sex.gr$score <- '1n.male'
  } else if (sample.sex == 'female') {
    hap.aut.gr <- hap.gr
    hap.sex.gr <- GRanges()
    cov.hap <- sum(width(hap.aut.gr))
    cov.sex.hap <- 0
  } else {
    hap.aut.gr <- hap.gr
    hap.sex.gr <- GRanges()
    cov.hap <- sum(width(hap.aut.gr))
    cov.sex.hap <- 0
  }
  
  ## Export diploid regions
  dip.gr$score <- '2n'
  hap.aut.gr$score <- '1n'
  multi.gr$score <- 'multi'
  cov.regions <- c(dip.gr, hap.aut.gr, hap.sex.gr, multi.gr)
  seqlevels(cov.regions) <- chromosomes
  cov.regions <- sort(cov.regions)
  names(mcols(cov.regions)) <- 'ploidy'
  if (!is.null(reference.fai)) {
    seqlengths(cov.regions)  <- seq.len[seqlevels(cov.regions)]
  } else {  
    seqlengths(cov.regions) <- NA
  }  
  
  ## Get number of bases per ploidy
  cov.dip <- sum(width(dip.gr))
  cov.multicov <- sum(width(multi.gr))
  summary.df <- data.frame(sample.id = sample.id,
                           cov.hap = cov.hap,
                           cov.sex.hap = cov.sex.hap,
                           cov.dip = cov.dip,
                           cov.multicov = cov.multicov
  )  
  ## Export
  return( list(cov.regions=cov.regions, summary=summary.df) )
}