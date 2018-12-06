source('functions.R')

fileFA = 'path/to/fasta.fa'
fileGTF = 'path/to/gtf'
dbSNP_path = 'path/to/dbsnp.vcf'
fa <- readDNAStringSet(fileFA, format='fasta')
gtf <- readGFF(fileGTF)
gtf = gtf[,c('seqid', 'type', 'gene_biotype', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'transcript_id')]

#--Define read quality
quality = "@ABCDEFGHIJKLMNOPQRSTUVWXYZ"
quals = paste(strsplit(quality,'', fixed=TRUE)[[1]][round(15 - 5/(seq(1,read_size, 1)**0.05))], collapse = '')

read_size = 75
coverage = 10
shift = 150
upDownStream = 100

# List of genes in GTF and their positions

ar1 = 1:10
ar2 = 1:10 + 3

for(i in c(1,6,12,18)){
  chr = as.character(i)
  print(chr)

  #----Find two sets of genes: overlapped and non overlapped
  genes = gtf[gtf$type=='gene' & gtf$gene_biotype=='protein_coding' & gtf_all$seqid==chr ,]

  g = GRanges(seqnames = as.character(genes$seqid), 
              ranges = IRanges(start = genes$start, end = genes$end), 
              strand = genes$strand, 
              gene_name = genes$gene_name)

  overlaps = as.data.frame(findOverlaps(g,g,ignore.strand=T))
  genes$count = overlaps %>% group_by(queryHits) %>% summarise(count=n()) %>% dplyr::select(count)

  listNonOverlappedGenes = genes[genes$count==1,'gene_name']
  listOverlappedGenes = unique(genes[genes$count>1,'gene_name'])


  #for (i in listNonOverlappedGenes[sample(1:length(listNonOverlappedGenes), 100, replace = F)]){
  for (i in listOverlappedGenes[1:min(length(listOverlappedGenes),100)]){
    
    #---Extract exon ranges from the GTF and introduce SNPs 
    gene = i 
    print(gene)
    g = gtf[gtf$gene_name==gene & gtf$seqid==chr,] 
    e = g[g$gene_name==gene & g$type=='exon',]
    exons = e %>% split(., .[,'transcript_id'])
    whereAreExons = where_are_exons(exons) #---List of transcript's IRanges exon positions 

    #---Take n random positions from every exon and create data.frame with all info
    snps = create_exon_snps_per_gene(gene, g) # gene_name, gtf, no_of_snps_per_exon
    dbSNPs = create_exon_dbsnps_per_gene(gene, g, dbSNP_path) # gene_name, gtf, path_to_dbSNPdatabase
    snps = rbind(dbSNPs, snps) %>% arrange(id) %>% distinct(pos, .keep_all = TRUE)

    #-----------Make two (phased) chromosomes--------------       
    ht1 = replaceLetterAt(fa[chr][[1]], at=snps[snps$phased=="1|0","pos"], letter=snps[snps$phased=="1|0","alt"])
    ht2 = replaceLetterAt(fa[chr][[1]], at=snps[snps$phased=="0|1","pos"], letter=snps[snps$phased=="0|1","alt"])

    #-----------Create mRNA reads-------------------
    #---Define number of trascripts
    transcriptsCount = create_transcript_count(exons, ar1, ar2) # sample rundom numbers from the arrays
    snps = snps %>% left_join(transcriptsCount, by = 'transcript')
    snps[snps$phased == "0|1", c("c1", "c2")] <- snps[snps$phased == "0|1", c("c2", "c1")]

    #---Create transripts with their counts
    mRNA1 = create_mrna(whereAreExons, ht1, transcriptsCount$c1)
    mRNA2 = create_mrna(whereAreExons, ht2, transcriptsCount$c2)

    DNA1 = create_dna(g, ht1, upDownStream)
    DNA2 = create_dna(g, ht2, upDownStream)

    #-----------Create mRNA and DNA FASTQ files
    rna = rbind(reads_to_fastq(mRNA1, read_size = 75, shift = 150, coverage = 10),reads_to_fastq(mRNA2, read_size = 75, shift = 150, coverage = 10))
    dna = rbind(reads_to_fastq(DNA1), reads_to_fastq(DNA2))

    write.table(snps, file='test2_rand2_overlap_noNoise_dbsnp.txt', row.names=FALSE, col.names=T, quote=F, append = T)
    write.table(dna[,1, drop=FALSE], file='test2_rand2_overlap_noNoise_dbsnp_dna_1.fq', row.names=FALSE, col.names=FALSE, quote=FALSE, append = TRUE)
    write.table(dna[,2, drop=FALSE], file='test2_rand2_overlap_noNoise_dbsnp_dna_2.fq', row.names=FALSE, col.names=FALSE, quote=FALSE, append = TRUE)
    write.table(rna[,1, drop=FALSE], file='test2_rand2_overlap_noNoise_dbsnp_rna_1.fq', row.names=FALSE, col.names=FALSE, quote=FALSE, append = TRUE)
    write.table(rna[,2, drop=FALSE], file='test2_rand2_overlap_noNoise_dbsnp_rna_2.fq', row.names=FALSE, col.names=FALSE, quote=FALSE, append = TRUE)
  
    #----Add noise--to transcripts-- 
    mRNA1 = mRNA1[rep(row.names(mRNA1), mRNA1$count),]
    mRNA1$count = 1
    mRNA2 = mRNA2[rep(row.names(mRNA2), mRNA2$count),]
    mRNA2$count = 1
  
    mRNA2$reads = unlist(lapply(mRNA2$reads, noise_to_transcripts,0.1))
    mRNA2$reads = unlist(lapply(mRNA2$reads, noise_to_transcripts,0.1))
  
    #-----------Create mRNA and DNA FASTQ files
    rna = rbind(reads_to_fastq(mRNA1, read_size = 75, shift = 150, coverage = 10),reads_to_fastq(mRNA2, read_size = 75, shift = 150, coverage = 10))
    dna = rbind(reads_to_fastq(DNA1), reads_to_fastq(DNA2))
  
    write.table(snps, file='test2_rand2_overlap_noise_dbsnp.txt', row.names=FALSE, col.names=T, quote=F, append = T)
    write.table(dna[,1, drop=FALSE], file='test2_rand2_overlap_noise_dbsnp_dna_1.fq', row.names=FALSE, col.names=FALSE, quote=FALSE, append = TRUE)
    write.table(dna[,2, drop=FALSE], file='test2_rand2_overlap_noise_dbsnp_dna_2.fq', row.names=FALSE, col.names=FALSE, quote=FALSE, append = TRUE)
    write.table(rna[,1, drop=FALSE], file='test2_rand2_overlap_noise_dbsnp_rna_1.fq', row.names=FALSE, col.names=FALSE, quote=FALSE, append = TRUE)
    write.table(rna[,2, drop=FALSE], file='test2_rand2_overlap_noise_dbsnp_rna_2.fq', row.names=FALSE, col.names=FALSE, quote=FALSE, append = TRUE)
    }
}


