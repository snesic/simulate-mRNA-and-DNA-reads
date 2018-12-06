#install if necessary
source("http://bioconductor.org/biocLite.R")

biocLite('ShortRead')
biocLite("Rsamtools")
biocLite('VariantAnnotation')

#load library
library(Rsamtools)
library(VariantAnnotation)
library(ShortRead)
library(rtracklayer)
library(tidyverse)
library(Biostrings)

#------------Functions to use

#---Take n random positions from every exon and create data.frame with all info
create_exon_snps_per_gene <- function(gene, gtf, nrand=2){
  
  g = gtf[gtf$gene_name==gene,]
  e = g[g$type=='exon',]
  exonsRanges = IRanges(start=e$start, end=e$end)
  exonsRanges = exonsRanges[!duplicated(exonsRanges)]
  
  posSNPs = as.vector(apply(as.matrix(exonsRanges), MARGIN = 1, function(x) sample(seq(x[1], x[1]+x[2]-1, 1),nrand, replace = T))) 
  posSNPs = posSNPs[!duplicated(posSNPs)]
  
  snps = data.frame(chr = chr, pos = posSNPs, ref = as.character(extractAt(fa[chr][[1]], at = IRanges(start = posSNPs, end = posSNPs))))
  bases = c('A','T','C','G')
  snps$alt = as.character(lapply(snps$ref, function(x, bases) sample(setdiff(bases, x),1), bases))
  
  snps$phased = sample(c("0|1","1|0"), nrow(snps), replace = T)
  snps$id=paste(snps$chr,':',snps$pos, '_', snps$ref,'/',snps$alt, sep='')
  
  #---Relate SNPs to trascripts
  
  snps$transcript = apply(snps, MARGIN = 1, function(x) paste(e[e$start<=x['pos'] & e$end>=x['pos'], 'transcript_id'], collapse =';'))
  snps = snps %>% separate_rows(transcript, sep=';')
  snps = snps[snps$transcript!='',]
  
  return(snps)
  
}

#---List of transcript's IRanges exon positions 
where_are_exons <- function(exons){
  whereAreExons = lapply(exons, function(x) IRanges(start=x$start, end=x$end))
  whereAreExons = lapply(whereAreExons, function(x) {end(x[end(x)==max(end(x))]) = end(x[end(x)==max(end(x))]) + 100; return(x)}) # Extend the last exon for 100bp 
  whereAreExons = lapply(whereAreExons, function(x) {start(x[start(x)==min(start(x))]) = start(x[start(x)==min(start(x))]) - 100; return(x)}) # Extend the last exon for 100bp 
  return(whereAreExons)
}
  
  #---------Insert dbsnp SNPs -------------------
#---Take n random positions from every exon and create data.frame with all info
create_exon_dbsnps_per_gene <- function(gene, gtf, dbSNP_path, no_snps=5){
  
  g = gtf[gtf$gene_name==gene,] 
  e = g[g$type=='exon',]
  exonsRanges = IRanges(start=e$start, end=e$end)
  exonsRanges = exonsRanges[!duplicated(exonsRanges)]
  
  genes = GRanges(seqnames = as.character(g$seqid), 
                  ranges = IRanges(start = g$start, end = g$end), 
                  strand = g$strand, 
                  gene_name = g$gene_name)
  
  param <- ScanVcfParam(fixed="ALT", 
                        info=NA, 
                        which = GRanges(genes[genes$gene_name==gene]))
  
  dbSNPs = readVcf(dbSNP_path, param = param)
  dbSNPs = dbSNPs[isSNV(dbSNPs)]
  
  dbSNPs = data.frame(chr = as.character(seqnames(rowRanges(dbSNPs))), 
                      id = row.names(dbSNPs), 
                      pos = start(ranges(rowRanges(dbSNPs))), 
                      ref = as.character(rowRanges(dbSNPs)$REF), 
                      alt = unlist(rowRanges(dbSNPs)$ALT),
                      phased = sample(c("0|1","1|0"), nrow(dbSNPs), replace = T))
  
  #---Relate SNPs to trascripts
  
  dbSNPs$transcript = apply(dbSNPs, MARGIN = 1, function(x) paste(e[e$start<=x['pos'] & e$end>=x['pos'], 'transcript_id'], collapse =';'))
  dbSNPs = dbSNPs %>% separate_rows(transcript, sep=';')
  dbSNPs = dbSNPs[dbSNPs$transcript!='',]
  
  return(dbSNPs)
}

#---Define number of trascripts
create_transcript_count <- function(exons, ar1, ar2){
  noPh1 = sample(ar1, length(exons), replace = T)
  noPh2 = sample(ar2, length(exons), replace = T)
  
  transcriptsCount = data.frame(transcript = names(exons), c1 = noPh1, c2 = noPh2)
  transcriptsCount$p.value = binom.test(min(sum(transcriptsCount$c1),sum(transcriptsCount$c1)), 
                                        sum(transcriptsCount$c1) + sum(transcriptsCount$c2))$p.value
  return(transcriptsCount)
}

#---Use exons per transcript (whereAreExons - transcript list) to create mRNA transcripts from the phased reference (ph - DNAString of a given phased contig), 
#---and introduce their count (c array of length no. of transcripts)
create_mrna <- function(whereAreExons, ph, c){
  mRNAph = lapply(whereAreExons, function(x) extractAt(ph, at=x)) # extract all exons 
  mRNAph = DNAStringSet(sapply(mRNAph, function(x) paste(x, collapse=''))) # merge them together 
  mRNAph = data.frame(reads = mRNAph, count = c, transcript = names(mRNAph))
  return(mRNAph)
}

#---Use gene range to create DNA from the phased reference (ph - DNAString of a given phased contig), 
#---and introduce their count (c array of length no. of transcripts)
create_dna <- function(gtf, ph, upDownStream){
  DNAph = extractAt(ph, at=IRanges(start = max(1,gtf[gtf$type=='gene','start'][1] - upDownStream), 
                                   end = min(length(ph), gtf[gtf$type=='gene','end'][1] + upDownStream),
                                   name=gtf[gtf$type=='gene','gene_name'][1]))
  DNAph = data.frame(row.names=names(DNAph), reads = DNAph, count = 1)
  return(DNAph)
}

#-----------Create FASTQ data.frame
reads_to_fastq <- function(read, read_size = 75, shift = 150, coverage = 10){
  
  # Sample positions of paired-end reads
  no_reads = round(width(read$reads)/read_size * coverage) * read$count
  ss = data.frame(size=pmax(1,(width(read$reads)-read_size-shift-1)), no_reads = no_reads)
  
  # Positions of paired-end reads
  where_1p = IRangesList(apply(ss, MARGIN = 1, function(x) IRanges(start = sample(1:x['size'], x['no_reads'], replace = T), width = read_size)))
  where_2p = shift(where_1p,shift)
  where_2p = restrict(where_2p, end = width(read$reads)) # Restrict the range not to go out the transcript
  
  p1 = extractAt(DNAStringSet(read$reads),at=where_1p)
  p2 = extractAt(DNAStringSet(read$reads),at=where_2p)
  names(p1) = read$transcript#rownames(read)
  names(p2) = read$transcript#rownames(read)
  p1 = unlist(p1)
  p2 = unlist(p2)
  p2 = reverseComplement(p2)
  
  reads = data.frame(pe1=paste('@gene:',names(p1), '_pos:', as.character(unlist(start(where_1p))),'/1\n', as.character(p1),'\n+\n', quals, sep=''),
                     pe2=paste('@gene:',names(p2), '_pos:', as.character(unlist(start(where_1p))),'/2\n', as.character(p2),'\n+\n', quals, sep=''),
                     stringsAsFactors=FALSE)
  
  return(reads) # prebaciti u fastq format 
}

#-- Noise to transcripts 

noise_to_transcripts <- function(m, percent=1){
  m = DNAString(m)
  bases = c('A','T','C','G')
  s = round(percent * nchar(p) / 100 )
  if(s==0) return(as.character(m))
  p = IRanges(start=sample(1:nchar(m), s, replace = T), width = 1)
  let = extractAt(m, at = p)
  let = unlist(lapply(let, function(x) (sample(setdiff(bases,x),1))))
  return(as.character(replaceLetterAt(m, at=start(p), letter = let)))
}


