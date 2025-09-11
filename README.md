# Amphioxus *de novo* mutation rate
[![bioinformatics](https://img.shields.io/badge/bioRxiv-Genomics-orange)](https://www.biorxiv.org/content/10.1101/2025.07.14.664012v3)  
The repo is the source code for the manuscript:
* **Title**: Germline *de* *novo* mutation rate of the highly heterozygous amphioxus genome
* **DOI**: https://doi.org/10.1101/2025.07.14.664012

1. [Requerements](#requirements)
2. [Contents](#contents)
3. [Analysis workflow](#analysis-workflow)
4. [Cites](#cites)
5. [Contact](#contact)  
   
## Requirements
The code was tested on CentOS Linux 7.  
Packages managed with mamba and singularity.  
ALL the code was run on the cluster using slurm.  
ALL softwares used in the scripts are open source and available on the internet.

## Contents
1. Genome assembly
2. Parents genomes alignment
3. Depth distribution similarity analysis
4. Alignment and variant calling
5. Simulation
6. Effective population size
7. PZMs detection
8. Amplicon sequencing
* plot
* scripts

## Analysis workflow

### 1. Genome_assembly

Genome servey and allele-aware diploid genome assembly using Platanus-allee (v2.2.2) (Kajitani, et al. 2019).  

Assembly evaluation using QUAST (v5.3.0) (Gurevich, et al. 2013) and  SeqKit (v2.9.0) (Shen, et al. 2016).  

### 2. Parents genomes alignment
Parents genomes alignment using minimap2 (v2.25) (Li, et al. 2021) and pick the potential gene conversion site using the script.  

### 3. Depth distribution similarity analysis
Depth distribution similarity analysis using python script and the callable genome detection using the scripts.  

### 4. Alignment and variant calling
Alignment was performed using BWA-MEM2 (v2.2.1) (Vasimuddin, et al. 2019) and variant calling was performed using BCFtools (v1.21) (Danecek, et al. 2021) and freebayes (v1.3.10) (Garrison and Marth 2012).  

### 5. Simulation
For evaluating the false negative rate (FNR) in our DNM detection pipeline, we generated simulated data based on callable genome regions of three selected offspring using ART (Huang, et al. 2012).  

### 6. Effective population size
To estimate nucleotide diversity at synonymous sites (πS), we re-analyzed previously published population data, and calculated πS using the script dNdSpiNpiS (v1.0) (https://kimura.univ-montp2.fr/PopPhyl/index.php?section=tools) with the default parameters.  

### 7. PZM detection
Identifying the candidate post-zygotic mutations (PZMs) were performed using the script.

### 8. Amplicon sequencing
To assess the reliability of DNMs across different data sources, we randomly selected a subset for validation by amplifying target regions with specific primers, followed by amplicon sequencing for genotyping. 

### plot
R scripts for generating the figures in the manuscript.  

### scripts
Independent scripts used in the analysis, including scripts for IGV visualization and gene conversion events identification.

## Cites
* Germline *de novo* mutation rate of the highly heterozygous amphioxus genome. Jing Xue, Lei Tao, Junwei Cao, Guang Li, Cai Li bioRxiv 2025.07.14.664012; doi: https://doi.org/10.1101/2025.07.14.664012

## Contact
For reporting issues or requests related to the package, please use GitHub Issues or write to jingxue1992@gmail.com.