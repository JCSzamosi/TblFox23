library('DEXSeq')
library('Rsubread')
library(tidyverse)
library(goseq)
setwd('~/Projects/Clients/Turnbull/CRISPR3_Take2/')

# Set up file locations and sample names
dirs = list.dirs('results/05_map2',recursive = FALSE)
dirs
samps = gsub('results/05_map2/', '', dirs, fixed = TRUE)
samps

# Set up the Gene ID lookup table (csv was produced manually from the GTF file)
gene_id_lookup = read.csv('~/Disk2/CommonData/ReferenceGenomes/Mouse/GRCm39/ncbi-genomes-2022-03-04/unique_gene_ids.csv')
dim(gene_id_lookup)
head(gene_id_lookup)
n_distinct(gene_id_lookup$gene_id)
n_distinct(gene_id_lookup$GeneID)
gene_id_v = gene_id_lookup$GeneID
names(gene_id_v) = gene_id_lookup$gene_id
duplicate_geneIDs = (gene_id_lookup 
                     %>% count(GeneID) 
                     %>% filter(n > 1) 
                     %>% left_join(gene_id_lookup))


# Pre-make the data frames to fill all the numbers in
exon_counts = matrix(NA, ncol = length(samps), nrow = 1520722)
colnames(exon_counts) = samps
head(exon_counts)

# Call it all once to set up the structure
samp = samps[1]
dr = dirs[1]
bamf = paste(dr, 'Aligned.sortedByCoord.out.bam', sep = '/')
gtf = '~/Disk2/CommonData/ReferenceGenomes/Mouse/GRCm39/ncbi-genomes-2022-03-04/GCF_000001635.27_GRCm39_genomic.gtf'
genome = '~/Disk2/CommonData/ReferenceGenomes/Mouse/GRCm39/ncbi-genomes-2022-03-04/GCF_000001635.27_GRCm39_genomic.fna'
call_fc = function(bamf, gtf, genome, useMF = FALSE){
    fc = Rsubread::featureCounts(bamf, 
                       annot.ext = gtf, isGTFAnnotationFile = TRUE,
                       GTF.featureType = 'exon', 
                       GTF.attrType = 'gene_id',
                       GTF.attrType.extra = 'exon_number', 
                       allowMultiOverlap = TRUE,
                       strandSpecific = 2,
                       genome = genome,
                       isPairedEnd = TRUE,
                       countReadPairs = FALSE,
                       requireBothEndsMapped = TRUE,
                       countChimericFragments = FALSE,
                       nthreads = 10,
                       verbose = TRUE,
                       useMetaFeatures = useMF)
    return(fc)
}
# fc = call_fc(bamf, gtf, genome)
# fc_genes = call_fc(bamf, gtf, genome, TRUE)
# stats = (fc$stat
#          %>% data.frame()
#          %>% column_to_rownames('Status'))
# colnames(stats) = samp
# exon_counts[,samp] = fc$counts
# 
# fc_lst = list()
# fc_lst[[samp]] = fc
# fc_lst_gene = list()
# fc_lst_gene[[samp]] = fc_genes

# # Call the same thing on all the mapping files
# 
# for (i in 2:length(dirs)){
#     samp = samps[i]
#     dr = dirs[i]
#     bamf = paste(dr, 'Aligned.sortedByCoord.out.bam', sep = '/')
#     fc = call_fc(bamf, gtf, genome)
#     exon_counts[,samp] = fc$counts
#     stats[fc$stat$Status,samp] = fc$stat$Aligned.sortedByCoord.out.bam
#     fc_lst[[samp]] = fc
#     fc_genes = call_fc(bamf, gtf, genome, useMF = TRUE)
#     fc_lst_gene[[samp]] = fc_genes
# }
# 
# # get gene lengths
# gene_lens = select(fc_genes$annotation, GeneID, Length)
# 
# save(list = c('exon_counts', 'stats', 'fc_lst', 'fc_lst_gene', 'gene_lens'), 
#      file = 'RData/feature_counts.RData')
load('RData/feature_counts.RData')
head(exon_counts)

# Get the gene and exon information for all the files
annot = fc_lst[[1]]$annotation
head(annot)

# Make sure exon numbers are unique and in a format DEXSeq can use
annot = (annot
         %>% data.frame()
         %>% group_by(GeneID)
         %>% mutate(exon_unique = as.character(1:length(GeneID)))
         %>% data.frame())
head(annot)
head(rownames(exon_counts))

# Bring in the sample data and clean it up
samdat = read.csv('data/sample_data/samdat.csv', row.names = 1)
samdat
samps = rownames(samdat)
samps = gsub('-','_',samps,fixed = TRUE)
rownames(samdat) = samps

samps = colnames(exon_counts)
samps = gsub('-','_',samps, fixed = TRUE)
colnames(exon_counts) = samps

# Do the genes with duplicate IDs have any counts associated with them?
rowSums(exon_counts[annot$GeneID %in% duplicate_geneIDs$gene_id,])

# They do. In some cases, it's a pretty big number
dup_cts = exon_counts[annot$GeneID %in% duplicate_geneIDs$gene_id,] 
rownames(dup_cts) = annot$GeneID[annot$GeneID %in% duplicate_geneIDs$gene_id]
dup_cts

annot %>% filter(GeneID %in% duplicate_geneIDs$gene_id)

# These genes can't be merged because they are different, but they have the same
# GeneID so they will be included in the DEXSeq analysis, but excluded from the
# goseq analysis.

# Run DEXSeq on undifferentiated samples
head(exon_counts)

samdat_undif = filter(samdat, Differentiation == 'Naive')
exon_ct_undif = exon_counts[,rownames(samdat_undif)]
samdat_undif
head(exon_ct_undif)

# dxd = DEXSeqDataSet(exon_ct_undif, samdat_undif, 
#                     design = ~sample + exon + exon:Treatment,
#                     featureID = annot$exon_unique,
#                     groupID = annot$GeneID)
# dxd = estimateSizeFactors(dxd, locfunc = genefilter::shorth)
# dxd = estimateDispersions(dxd)
# save(dxd, file = 'RData/dxd_undif.RData')
# plotDispEsts(dxd)
# dxd = testForDEU(dxd)
# dxd = estimateExonFoldChanges(dxd, fitExpToVar = 'Treatment')
# save(dxd, file = 'RData/dxd_undif_deu.RData')
load('RData/dxd_undif_deu.RData')
res = DEXSeqResults(dxd)

# Goseq stuff

# Identify genes with significanly different (padj < 0.05) exon usage in at 
# least one exon
head(res)
sig_genes = (res
           %>% data.frame()
           %>% filter(padj < 0.05)
           %>% select(groupID)
           %>% unique())$groupID

# Check if any of the problem genes have DEU
any(sig_genes %in% duplicate_geneIDs$gene_id)

# Make a binary vector of all genes, 1 for significant, 0 for non
genes = unique(annot$GeneID)
deg_v = genes %in% sig_genes
deg_v = as.numeric(deg_v)
names(deg_v) = genes
deg_v = deg_v[!(genes %in% duplicate_geneIDs$gene_id)]
names(deg_v) = gene_id_v[names(deg_v)]
head(deg_v)
any(is.na(names(deg_v)))

# Get the bias and length estimations for goseq
np = nullp(deg_v, 'mm39', 'refGene')

# Run goseq
gs = goseq(np, 'mm39', 'refGene')

head(gs)
dim(gs)
gs_fixed = (gs
            %>% mutate(padj = p.adjust(over_represented_pvalue, 'BH')))
head(gs_fixed)

# Count the significant categories
sum(sum(gs_fixed$padj < 0.05))
