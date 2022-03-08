#!/bin/bash

# Second pass
ref=/home/jcszamosi/Disk2/CommonData/ReferenceGenomes/Mouse/GRCm39/ncbi-genomes-2022-03-04/STAR2.7.9_100
gtf=/home/jcszamosi/Disk2/CommonData/ReferenceGenomes/Mouse/GRCm39/ncbi-genomes-2022-03-04/GCF_000001635.27_GRCm39_genomic.gtf

indir=data/
SJs=`ls results/04_map/*/*tab`
outdir=results/05_map2/
r1=_R1.fastq.gz
r2=_R2.fastq.gz

for R1 in $indir/*R1*
do

	sample=`basename $R1 $r1`
	R2=$indir/$sample$r2

	echo mkdir -p $outdir/$sample/
	mkdir -p $outdir/$sample/

	echo STAR --runThreadN 8 \
		--genomeDir $ref \
		--readFilesIn $R1 $R2 \
		--sjdbGTFfile $gtf \
		--sjdbFileChrStartEnd $SJs \
		--sjdbOverhang 49 \
		--outFilterMultimapNmax 20 \
		--alignSJoverhangMin 8 \
		--outFilterMismatchNmax 999 \
		--outFilterMismatchNoverReadLmax 0.04 \
		--alignIntronMin 20 \
		--alignIntronMax 1000000 \
		--alignMatesGapMax 1000000 \
		--alignSJDBoverhangMin 1 \
		--outFileNamePrefix $outdir/$sample/ \
		--outSAMtype BAM SortedByCoordinate \
		--peOverlapNbasesMin XX \
		--peOverlapMMp 0.1 \
		--readFilesCommand gunzip -c

	STAR --runThreadN 8 \
		--genomeDir $ref \
		--readFilesIn $R1 $R2 \
		--sjdbGTFfile $gtf \
		--sjdbFileChrStartEnd $SJs \
		--sjdbOverhang 49 \
		--outFilterMultimapNmax 20 \
		--alignSJoverhangMin 8 \
		--outFilterMismatchNmax 999 \
		--outFilterMismatchNoverReadLmax 0.04 \
		--alignIntronMin 20 \
		--alignIntronMax 1000000 \
		--alignMatesGapMax 1000000 \
		--alignSJDBoverhangMin 1 \
		--outFileNamePrefix $outdir/$sample/ \
		--outSAMtype BAM SortedByCoordinate \
		--peOverlapNbasesMin XX \
		--peOverlapMMp 0.1 \
		--readFilesCommand gunzip -c
done
