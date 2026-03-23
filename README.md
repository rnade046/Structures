# LESuMoN: Local Enrichment of Structural Motifs in biological Networks

## Overview

#### Description
LESuMoN is a tool for RNA structural motif discovery in the context of protein co-localization networks (e.g. Human Cell Map, Go, Knight *et al.* Nature 2021). Leverage RNAfold (citation) and BayesPairing2 (citation) to predict structural motif int 3'UTR sequences associated to protein in network. The algorithm measures the clustering of proteins with a shared 3'UTR structrual motif within the network. The significance of this clustering measure is then estimated from a normal distribution approximated from a Monte Carlo Sampling Distribution. To correct for multiple hypothesis testing, we assess the false discovery rate at various thresholds against a null model of randomized sequences.

## Prerequisites
* Java Version 8+
* Required Java libraries, included in the repository: 
   * `commons-math3-3.6.1.jar`
   * `json-20230227.jar`
* Python 3.10+
* [RNABayesPairing2](https://jwgitlab.cs.mcgill.ca/sarrazin/rnabayespairing2) and it's dependencies
* [ViennaRNA package](https://www.tbi.univie.ac.at/RNA/)

## Input files
Example input files are available in our [Zenodo repository](https://doi.org/10.5281/zenodo.19162062).

#### 1. Sequences (e.g. `hg38_3utr.fasta` and `hg38_cds.fasta`)
Two separate FASTA files are required: one for 3'UTR sequences and one for coding sequences. These files should contain sequences corresponding to proteins in the network. Additional sequences may be included, but they will be ignored if they are not associated with proteins in the input ID file.

#### 2. Protein Ids (e.g. `mapping.tsv`) 
Tab-separated file mapping protein names used in the network to RefSeq identifiers used in the FASTA files. If a protein has multiple identifiers, they can be listed on separate lines. [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) is a useful tool for generating this file.

| Protein Name | RefSeqIds |
|---|---|
| ABCB1 | NM_001348946 |
| ABCB1 | NM_001348944 |
| BDP1 | NM_018429 |

#### 3. Network (e.g. `network.tsv`)
User provided .tsv file of a protein-protein interaction network. 

   Protein1 | Protein2 | Score
   --- | --- | ---
   ABC | ACC | 0.75

#### 4. Properties file (e.g. `LESuMoN.properties` included in repo) 

This file is passed to the program as a command-line argument. It must be updated with the correct working directory, file paths, and analysis parameters.

All file paths below are specified relative to `working_directory` unless otherwise noted.

| Property Name           | Type      | Description |
|------------------------|----------|-------------|
| `project_name`         | `String` | Unique mame identifying your project
| `nullModel`       | `Boolean` | Specify if using randomized sequences (`true`) or biological sequences (`false`) Example: `true` |
| `working_directory`   | `String` | Directory which contains input files and where output files will be stored. Example: `/path/to/project/sequences/` |
| `codingSeqFile`   | `String` | Path to fasta file relative to working directory. Example: `input_files/sequences.fasta` |
| `utrSeqFile`   | `String` | Path to fasta file relative to working directory. Example: `input_files/sequences.fasta` |
| `geneIdsFile`   | `String` | Path to geneIdsFile relative to working directory. Example: `input_files/geneIds.tsv` |
| `networkRepositoryFile`   | `String` | Path to network file relative to working directory. Example: `input_files/network.tsv` |
| `clusteringMeasure` | `Integer` | clustering measure to use in analysis, TPD (`0`), NWTPD(`4`) . Example: `4` |
| `bpScoreThreshold`  | `Float` | Minimum BayesPairing2 score to consider in analysis. Example: `4.0` |
| `rescaleThreshold`  | `Float` | Minimum percentile of all scores to consider for rescaling procedure. Required for NWTPD. Example: `0.85` |
| `lowerBoundToSample` | `Integer` | Smallest sample size to assess for Monte Carlo Sampling procedure. Recommended: `3` |
| `upperBoundToSample` | `Integer` | Largest sample size to assess for Monte Carlo Sampling procedure. Recommended: `2000` |
| `numberOfSamplings` | `Integer` | Number of times to sample network for Monte Carlo Sampling. Recommended: `1000000` |
    
# Pipeline
![graphical representation of pipeline](pipeline-workflow.png)

To execute the full pipeline, run LESuMoN once on the biological sequences and run LESuMoN_rand multiple times on randomized sequences. Fifteen null-model runs were used in our analyses.

## Set up
Follow these steps to set up the LESuMoN project on your local machine or server:

#### 1. Clone the repository

```bash
git clone https://github.com/LavalleeAdamLab/LESuMoN
```
#### 2. Organize file structure
Create a project directory outside the repository to store input files and analysis outputs.

```text
/path/to/project/
│
├── LESuMoN/                      # Cloned repository
│   ├── README.md
│   ├── LESuMoN.sh
│   ├── LESuMoN_rand.sh
│   ├── fdr.sh
│   ├── summarize_motifs.sh
│   ├── commons-cli-1.9.0.jar
│   ├── commons-math3-3.6.1.jar
│   ├── json-20230227.jar
│   ├── post/
│   ├── localEnrich/
│   ├── annotateMotifs/
│   ├── fdr/
│   └── randomizeSeq/
│
├── input_files/                  # Shared input data
│   ├── coding_sequences.fasta
│   ├── utr_sequences.fasta
│   ├── network.tsv
│   └── geneIDs.tsv
│
├── biological_run/               # LESuMoN run on biological sequences
│   └── LESuMoN.properties
│
└── null_model_runs/              # Repeated null-model runs
    ├── run_01/
    │   └── LESuMoN.properties
    ├── run_02/
    │   └── LESuMoN.properties
    ├── run_03/
    │   └── LESuMoN.properties
    └── ...

```

#### 3. Configure the Properties Files

#### 4. Run All Commands from the Project Directory

```bash
# Program must be executed from LESuMoN/ repository
cd path/to/directory/LESuMoN/
```

## Biological Sequence Analysis
Run LESuMoN structural motif enrichment on the biological sequences.

#### Implementation
```bash session
bash LESuMoN.sh /path/to/project/biological_run/LESuMoN.properties

```
Note : given the size of motifs to test and time complexity of this algorithm, it is recommended to execute this part of the program on a server with job parallelization capabilities.
#### Output
* `annotations.tsv`: Maps proteins associated with predicted structural motifs and their BayesPairing2 scores. Required for motif summarization.
* `pvals_fwd.tsv`: List of computed significance values, required for FDR estimation
* `localEnrichment_MotifDetails.tsv`: Summary of clustering measures and significance values for each motif.

## 2. Null Model Analysis
Run LESuMoN on randomized sequences to estimate the null distribution. This step should be repeated multiple times. Fifteen runs were used in our analyses.

#### Execution
```bash
bash LESuMoN_rand.sh /path/to/project/null_model_runs/run_01/LESuMoN.properties
```

#### Output
* Randomized sequences in a FASTA format (e.g.`shuffled_seq.fasta`)
* `pvals_rand.tsv` required for FDR calculation

## 3. FDR estimation 
Estimate the false discovery rate by comparing motif clustering significance in biological sequences against motif clustering significance in randomized sequences.

#### Execution
```bash
bash fdr.sh
```
#### Output
* `fdr.tsv` Table file indicating the FDR estimated at given *p*-value thresholds and number of motifs significantly clustered at each threshold.

## 4. Motif summary
Generate the final summary table for significantly clustered motifs at a user-defined p-value threshold (e.g. 0.001).

#### Execution
```bash
bash summarize_motifs.sh 0.001
```
#### Output
* Summary file output as .tsv format (e.g.`SupplementaryTable_S1.xlsx`).