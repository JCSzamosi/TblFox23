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
#      file = 'feature_counts.RData')
load('feature_counts.RData')
head(exon_counts)
annot = fc_lst[[1]]$annotation
head(annot)
annot = (annot
         %>% data.frame()
         %>% group_by(GeneID)
         %>% mutate(exon_unique = as.character(1:length(GeneID)))
         %>% data.frame())
head(annot)
head(rownames(exon_counts))

samdat = read.csv('data/sample_data/samdat.csv', row.names = 1)
samdat
samps = rownames(samdat)
samps = gsub('-','_',samps,fixed = TRUE)
rownames(samdat) = samps

samps = colnames(exon_counts)
samps = gsub('-','_',samps, fixed = TRUE)
colnames(exon_counts) = samps

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
# save(dxd, file = 'dxd_undif.RData')
# plotDispEsts(dxd)
# dxd = testForDEU(dxd)
# dxd = estimateExonFoldChanges(dxd, fitExpToVar = 'Treatment')
# save(dxd, file = 'dxd_undif_deu.RData')
load('dxd_undif_deu.RData')

res = DEXSeqResults(dxd)
# ================================================

# Goseq stuff

head(res)
sig_genes = (res
           %>% data.frame()
           %>% filter(padj < 0.05)
           %>% select(groupID)
           %>% unique())$groupID

# Bring in the GeneID sheet
gene_id_lookup = read.csv('~/Disk2/CommonData/ReferenceGenomes/Mouse/GRCm39/ncbi-genomes-2022-03-04/unique_gene_ids.csv')
dim(gene_id_lookup)
head(gene_id_lookup)
n_distinct(gene_id_lookup$gene_id)
n_distinct(gene_id_lookup$GeneID)
gene_id_lookup %>% count(GeneID) %>% filter(n > 1) %>% left_join(gene_id_lookup)
gene_id_v = gene_id_lookup$GeneID
names(gene_id_v) = gene_id_lookup$gene_id
genes = unique(annot$GeneID)
deg_v = genes %in% sig_genes
deg_v = as.numeric(deg_v)
names(deg_v) = genes
deg_v = deg_v[!(genes %in% c('Erdr1_1','Erdr1','G530011O06Rik_1','G530011O06Rik'))]
names(deg_v) = gene_id_v[names(deg_v)]
head(deg_v)
any(is.na(names(deg_v)))
np = nullp(deg_v, 'mm39', 'refGene')

gs = goseq(np, 'mm10', 'refGene')

head(gs)
dim(gs)
gs_fixed = (gs
            %>% mutate(padj = p.adjust(over_represented_pvalue, 'BH')))
head(gs_fixed)
#================================
gtf = read.delim(file = gtf, header = FALSE, skip = 4)
gtf_final = gtf$V9
head(gtf_final)
gtf_final = gtf_final[-length(gtf_final)]
gtf1 = gtf_final[1]
gtf1
l1 = strsplit(gtf1, ';')

# ================================================
res_non_zero = filter(data.frame(res), exonBaseMean > 0)

res_non_zero_for_goseq = (res_non_zero
                          %>% group_by(groupID)
                          %>% summarize(padj = min(padj))
                          %>% mutate(Sig = padj < 0.05)
                          %>% select(-padj)
                          %>% unique())
write.csv(res_non_zero_for_goseq, file = 'res_non_zero_for_goseq.csv')
dim(res_non_zero)
res_use = filter(res_non_zero, !is.na(padj))
dim(res_use)
head(res_use)
res_sig = filter(res_use, padj < 0.05)
dim(res_sig)
head(res_sig)
n_distinct(res_sig$groupID)
summary(res_use$log2fold_WT_KO)
res_sig_large = res_sig %>% filter(abs(log2fold_WT_KO) >= 0.5, 
                                   exonBaseMean >= 10)
dim(res_sig_large)
n_distinct(res_sig_large$groupID)
unique(res_sig_large$groupID)

res_sig_large_up = filter(res_sig_large, log2fold_WT_KO > 0)
dim(res_sig_large_up)
unique(res_sig_large_up$groupID)

res_sig_large
write.csv(res_sig_large, file = 'res_sig_large.csv')
res_sig_large_for_goseq = select(res_sig_large, groupID, padj)
res_sig_large_for_goseq = (res_sig_large_for_goseq
                           %>% group_by(groupID)
                           %>% summarize(padj = min(padj))
                           %>% mutate(Sig = padj<0.05)
                           %>% select(-padj)
                           %>% unique())
head(res_sig_large_for_goseq)
write.csv(res_sig_large_for_goseq, file = 'res_sig_large_for_goseq.csv')

head(gene_lens)

genes_for_goseq = (gene_lens
                   %>% mutate(DE = GeneID %in% res_sig$groupID))
head(genes_for_goseq)

genes_for_goseq = (genes_for_goseq
                   %>% mutate(GeneID = toupper(GeneID)))

DEgenes = genes_for_goseq$DE
names(DEgenes) = genes_for_goseq$GeneID

np = nullp(DEgenes, 'mm10', 'geneSymbol')

# plotDEXSeq(res, 'Mff', legend = TRUE, FDR = 0.05, fitExpToVar = 'Treatment',
           # xlim = c(0,1), ylim = c(0,1))
