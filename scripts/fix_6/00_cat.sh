#! /bin/bash

# Run from top directory

orig=orig_data
l1=L001_R1_001.fastq.gz

for file in $orig/6*$l1
do
	bnm=`basename $file $l1`
	echo $bnm
	suf=R1.fastq.gz
	cat $orig/$bnm*R1* > data/$bnm$suf
	suf=R2.fastq.gz
	cat $orig/$bnm*R2* > data/$bnm$suf
done
