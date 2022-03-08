Differential Expression and Differential Splicing after Fox2/3 Knockdown
========================================================================

Steps Overview
--------------

0. Download
	* Downloaded the files by running `bs list biosample` to get a list of
	all the biosamples I had access to, then used `bs download biosample -i
	$sampleID` to download them. To download all the samples in a project at
	once I used `bs list projects` and then `bs download project -i $projectID`
1. Sequence Prep
	1. FastQC 
		* [script](./scripts/01_fastqc.sh)
		* [results](./results/01_fastqc)
		* [logs](./logs/01_fastqc.screenlog)
	2. cutadapt
		* [script](./scripts/02_cutadapt.sh)
		* [results](./results/02_cutadapt)
		* [logs](./logs/02_cutadapt.screenlog)
	3. FastQC pass 2 
		* [script](./scripts/03_fastqc.sh)
		* [results](./results/03_fastqc)
		* [logs](./logs/03_fastqc.screenlog)
2. Mapping
	1. Index the genome with STAR 
		* [script](./scripts/index_genome.sh)
		* [logs](./logs/index_genome.screenlog)
	2. Map first pass with STAR 
		* [script](./scripts/04_map.sh)
		* [results](./results/04_map)
		* [logs](./logs/04_map.screenlog).
	3. Second pass map with STAR

Details
-------

### Sequence Prep

FastQC just displays the sequencing quality and other QC metrics for the
sequences. Everything looked very good.

I used cutadapt to quality-trim the sequences and also to trim off sequencing
adapters. Changes I made from your workflow: I supplied adapter sequences to
cutadapt, so that it could trim them off, and I raised the quality threshold
from 20 to 28.

I then ran FastQC a second time, to check how many sequences were lost. All of
the samples had excellent retention rates.

It was at this stage that I discovered that sample 6 was corrupted during
download and had to re-download it and re-run it to here. Once I did that, it
was fine. Scripts to do thaat were in [fix_6](./scripts/fix_6/).

### Mapping

I downloaded the mouse reference genome GRCm39 from NCBI on March 4, 2022 along
with its gtf file. I indexed this using STAR v2.7.9a, generating a splice
junction database from the GTF file. This is different from your workflow in
that I used the GTF file to provide _a priori_ information about the location of
splice junctions.

I mapped using some different parameters. For instance, I filtered out reads
that mapped to more than 20 loci, or who's mapping score was lower than 8. The
full list of my parameter settings are in the [script](./scripts/04_map.sh).
This mapping produced a new set of splice junctions for each sample.

I then did a second pass mapping of each sample, using the tables of splice
junctions generated from _all_ the samples, so that all reads that split across
splice junctions found in other samples could be detected.

