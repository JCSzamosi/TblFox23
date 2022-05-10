Differential Expression and Differential Splicing after Fox2/3 Knockdown
========================================================================

Steps Overview
--------------

0. Download
	1. Downloaded the sample files by running `bs list biosample` to get a list 
	of all the biosamples I had access to, then used `bs download biosample -i
	$sampleID` to download them. To download all the samples in a project at
	once I used `bs list projects` and then `bs download project -i $projectID`
	2. Downloaded the **mm39** reference genome from NCBI
	[here](https://www.ncbi.nlm.nih.gov/genome/52?genome_assembly_id=992563)
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
		* [script](./scripts/05_map2.sh)
		* [results](./results/05_map2)
		* [logs](./logs/05_map2.screenlog)
3. Analysis
	* All analysis is performed in the [R analysis script](./scripts/06_all_analysis.R)
	1. Count features using `featureCount()` from `Rsubreads`
	2. Test for differential exon usage (DEU) and differential gene expression
	(DE) using `DEXSeq`
	3. Pathway analysis of both DEU and DE using `goseq`

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
that mapped to more than 20 loci, or whose mapping score was lower than 8. The
full list of my parameter settings are in the [script](./scripts/04_map.sh).
This mapping produced a new set of splice junctions for each sample.

I then did a second pass mapping of each sample, using the tables of splice
junctions generated from _all_ the samples, so that all reads that split across
splice junctions found in other samples could be detected.

### Analysis

#### Read Counting

I used `Rsubread::featureCounts()` to count reads assigned to each exon and
gene. Full parameters can be found in the [script](./script/06_all_analysis.R).

#### DE/DEU

I used `DEXSeq` to analyze the differential gene expression and differential
exon usage. I analyzed Naive and Differentiated cells separately with the
variable of interest being scramble vs. guide sequence. I ran with default
parameters and the model specified was `~sample + exon + exon:Treatment`. I 
estimated exon fold changes using `estimateExonFoldChange(dxd, fitExpToVar =
'Treatment')`.

#### Goseq

I used Goseq to identify the gene ontologies/functional categories that are
over-represented in terms of differential exon usage. To do this I collapsed the
features to genes, and I marked a gene as differentially used if **at least one
exon** in that gene was differentially used.

I downloaded the UCSC mm39 refGene database from Bioconductor:
`TxDb.Mmusculus.UCSC.mm39.refGene`. This database uses a different naming
convention than the `gene_id` field in the GTF file, so I created a lookup table
to convert the NCBI gene IDs to the refGene IDs, which were also present in the
GTF file under the field `db_xref GeneID`. Of the 2.9M features in the GTF file,
63 did not include GeneID tags. These consisted of tRNA genes and mitochondrial
genes including CytB, COX genes, NADH dehydrogenase genes, etc. These 63 genes
were omitted from `goseq` analysis. Additionally, two pairs of genes (AAA and
BBB) were omitted because each pair shared a single GeneID value, despite being
on different chromosomes. These genes were both located on the sex chromosomes
(one of each on the X and Y chromosomes) and neither contained any DEU.
