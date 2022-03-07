#!/bin/bash -i

echo source ~/miniconda3/bin/activate cutadaptenv
source ~/miniconda3/bin/activate cutadaptenv

indir=data/
outdir=results/02_cutadapt/

in1=R1.fastq.gz
in2=R2.fastq.gz

out1=R1_trimmed.fastq.gz
out2=R2_trimmed.fastq.gz

adapt="file:/home/jcszamosi/Disk2/CommonData/LibFiles/cutadapt_adapters.fna"

echo mkdir -p $outdir
mkdir -p $outdir

for R1 in $indir/6*$in1
do
	bnm=`basename $R1 $in1`
	R2=$indir/$bnm$in2
	
	echo cutadapt -a $adapt -A $adapt \
	-o $outdir/$bnm$out1 -p $outdir/$bnm$out2 \
	-j 8 \
	$R1 $R2 \
	-q 28 -m 50

	cutadapt -a $adapt -A $adapt \
	-o $outdir/$bnm$out1 -p $outdir/$bnm$out2 \
	-j 8 \
	$R1 $R2 \
	-q 28 -m 50
done

echo conda deactivate
conda deactivate
