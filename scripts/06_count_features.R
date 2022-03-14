library('Rsubread')
setwd('~/Projects/Clients/Turnbull/CRISPR3_Take2/')
test_count = featureCounts(
    'results/05_map2/1-C1P_S1/Aligned.sortedByCoord.out.bam',
    annot.ext = '~/Disk2/CommonData/ReferenceGenomes/Mouse/GRCm39/ncbi-genomes-2022-03-04/GCF_000001635.27_GRCm39_genomic.gtf',
    isGTFAnnotationFile = TRUE,
    strandSpecific = 2,
    juncCounts = TRUE,
    genome = '~/Disk2/CommonData/ReferenceGenomes/Mouse/GRCm39/ncbi-genomes-2022-03-04/GCF_000001635.27_GRCm39_genomic.fna',
    isPairedEnd = TRUE,
    countReadPairs = TRUE,
    requireBothEndsMapped = TRUE,
    countChimericFragments = FALSE,
    nthreads = 8,
    reportReadsPath = 'results/06_featureCounts/1-C1P_S1/',
    tmpDir = 'results/06_featureCounts/tmp/',
    verbose = TRUE,
    reportReads = 'CORE')

test_count = featureCounts(
    'results/05_map2/2-S1_S2/Aligned.sortedByCoord.out.bam',
    annot.ext = '~/Disk2/CommonData/ReferenceGenomes/Mouse/GRCm39/ncbi-genomes-2022-03-04/GCF_000001635.27_GRCm39_genomic.gtf',
    isGTFAnnotationFile = TRUE,
    strandSpecific = 2,
    juncCounts = TRUE,
    genome = '~/Disk2/CommonData/ReferenceGenomes/Mouse/GRCm39/ncbi-genomes-2022-03-04/GCF_000001635.27_GRCm39_genomic.fna',
    isPairedEnd = TRUE,
    countReadPairs = TRUE,
    requireBothEndsMapped = TRUE,
    countChimericFragments = FALSE,
    nthreads = 8,
    reportReadsPath = 'results/06_featureCounts/1-C1P_S1/',
    tmpDir = 'results/06_featureCounts/tmp/',
    verbose = TRUE,
    countMultiMappingReads = FALSE)

