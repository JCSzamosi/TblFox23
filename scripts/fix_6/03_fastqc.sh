#! /bin/bash

dat=results/02_cutadapt/
outdr=results/03_fastqc/

echo mkdir -p $outdr
mkdir -p $outdr
echo fastqc -o $outdr -t 10 $dat/6*
fastqc -o $outdr -t 10 $dat/6*
