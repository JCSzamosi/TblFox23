#!/bin/bash

# Mapping
ref=/home/jcszamosi/Disk2/CommonData/ReferenceGenomes/Mouse/GRCm39/ncbi-genomes-2022-03-04/STAR2.7.9_100
gtf=/home/jcszamosi/Disk2/CommonData/ReferenceGenomes/Mouse/GRCm39/ncbi-genomes-2022-03-04/GCF_000001635.27_GRCm39_genomic.gtf

indir=data/
outdir=results/04_map/
r1=_R1.fastq.gz
r2=_R2.fastq.gz

for R1 in $indir/*R1*
do

	sample=`basename $R1 $r1`
	R2=$indir/$sample$r2

	echo mkdir -p $outdir/$sample
	mkdir -p $outdir/$sample

# Run the alignments
	echo STAR --runThreadN 8 \
		--genomeDir $ref \
		--readFilesIn $R1 $R2 \
		--readFilesCommand gunzip -c \
		--outFilterMultimapNmax 20 \
		--alignSJoverhangMin 8 \
		--outFilterMismatchNmax 999 \
		--outFilterMismatchNoverReadLmax 0.04 \
		--alignIntronMin 20 \
		--alignIntronMax 1000000 \
		--alignMatesGapMax 1000000 \
		--alignSJDBoverhangMin 3 \
		--outFileNamePrefix $outdir/$sample/ \
		--outSAMtype BAM SortedByCoordinate \
		--peOverlapMMp 0.1 

	STAR --runThreadN 8 \
		--genomeDir $ref \
		--readFilesIn $R1 $R2 \
		--readFilesCommand gunzip -c \
		--outFilterMultimapNmax 20 \
		--alignSJoverhangMin 8 \
		--outFilterMismatchNmax 999 \
		--outFilterMismatchNoverReadLmax 0.04 \
		--alignIntronMin 20 \
		--alignIntronMax 1000000 \
		--alignMatesGapMax 1000000 \
		--alignSJDBoverhangMin 3 \
		--outFileNamePrefix $outdir/$sample/ \
		--outSAMtype BAM SortedByCoordinate \
		--peOverlapMMp 0.1 

done


# use 2-pass mapping (section 8 in the manual)
# use merging and mapping of overlapping paired ends (section 9 in the manual)
