#! /bin/bash

dat=data/
outdr=results/01_fastqc/

echo mkdir -p $outdr
mkdir -p $outdr
echo fastqc -o $outdr -t 10 $dat/*
fastqc -o $outdr -t 10 $dat/*
