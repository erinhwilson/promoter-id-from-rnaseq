# promoter-id-from-rnaseq
A computational framework to identify promoter sequences from RNA-seq datasets

## Overview
This framework contains 3 main steps:
1. Identify a set of top genes from RNA-seq data
1. Extract sequence regions upstream of a set of genes
1. Search upstream regions for promoter motif structures and make predictions

And 3 main outputs:
1. A fasta file of the exact -35::-10 of the best promoter prediction for each locus
1. A tab-delimited file summarizing the top 3 best promoters for each locus, with Margin of Victories scores to indicate the robustness of each selections (small margin of victory indicates that there were several close predictions worth reviewing)
1. A tab-delimitied file summarzing all possible promoter predictions for each locus

<img src="comp_framework_diagram.jpg" alt="Computational Framework" width="600"/>

## Workflow Instructions

### Obtain RNA-seq data matrix
Obtain a data matrix where each row is a genome locus, each column is an RNA-seq sample, and each value is the RNA-seq read count in transcripts per million (TPM). The choice of workflow that transforms raw RNA-seq data (fastq) to such a matrix is flexible. Here we used barrelseq (code available here).

### Select a set of top genes
Once data is properly formatted as a matrix of genes by experimental samples, we can use these data to select a set of highly expressed genes that remain high across conditions.

Inputs:
1. TPM data matrix
1. sample2condition mapping file. 
  * The columns of the TPM data matrix should reflect unique sample names (distinct RNA-seq experiments). Some of these samples may be replicates, or just separate experiments that fall under the same experimental category (e.g., "Low Methane", or "High Copper" etc). The experimental conditions can be called anything, but the sample2condition file is the formal way to specify which sample belongs to which category.

example format for sample2condition.txt
| | |
| ------------ | ------------- |
| 5GB1_FM12_TR2 | lowCH4 |
| 5GB1_FM23_TR3 | MeOH |
| 5GB1_FM40_T0_TR1 | NoCu |
| 5GB1_FM34_T8_TR1 | HighCu |
