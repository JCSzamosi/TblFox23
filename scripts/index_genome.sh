#!/bin/bash

# Index the reference genome
out=/home/jcszamosi/Disk2/CommonData/ReferenceGenomes/Mouse/GRCm39/ncbi-genomes-2022-03-04/STAR2.7.9_100
ref=/home/jcszamosi/Disk2/CommonData/ReferenceGenomes/Mouse/GRCm39/ncbi-genomes-2022-03-04/GCF_000001635.27_GRCm39_genomic.fna
gtf=/home/jcszamosi/Disk2/CommonData/ReferenceGenomes/Mouse/GRCm39/ncbi-genomes-2022-03-04/GCF_000001635.27_GRCm39_genomic.gtf

echo STAR --runThreadN 10 \
	--runMode genomeGenerate \
	--genomeDir $out \
	--sjdbGTFfile $gtf \
	--sjdbOverhang 49 \
	--limitGenomeGenerateRAM 24000000000 \
	--genomeSAsparseD 2 \
	--genomeFastaFiles $ref 

STAR --runThreadN 10 \
	--runMode genomeGenerate \
	--sjdbGTFfile $gtf \
	--sjdbOverhang 49 \
	--genomeDir $out \
	--limitGenomeGenerateRAM 24000000000 \
	--genomeSAsparseD 2 \
	--genomeFastaFiles $ref 
