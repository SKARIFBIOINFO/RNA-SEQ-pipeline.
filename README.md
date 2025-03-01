# RNA-SEQ-pipeline.


# Overview
This repository provides a complete RNA-Seq analysis pipeline designed to process and analyze sequencing data from Escherichia coli samples under different experimental conditions. The pipeline encompasses data preprocessing, alignment, quantification, and differential gene expression analysis, utilizing widely adopted bioinformatics tools within an Anaconda environment and automated through Bash scripting.


# Background

In this project, we aim to develop a comprehensive RNA-Seq analysis pipeline to investigate how Escherichia coli (E. coli) adapts its gene expression when utilizing different carbon sources, specifically glucose and glycerol. By implementing this pipeline, we seek to identify differentially expressed genes and understand the bacterial response to different environmental factors, thereby enhancing our knowledge of E. coli's adaptive mechanisms.

Data Source:

The dataset utilized for this analysis is sourced from the Gene Expression Omnibus (GEO) under the accession number **GSE94117**. This study, titled "RNA-Seq analysis of E. coli under multiple growth conditions," examines the bacterial transcriptome across various environmental settings.

Sample Details:

The study comprises 152 samples, each representing E. coli cultures subjected to different growth conditions, including:

Growth Phases: Exponential, stationary, and late-stationary phases.

Carbon Sources: Glucose, gluconate, lactate, and glycerol.

Salt Stresses: Variations in magnesium (Mg²⁺) and sodium (Na⁺) concentrations.

For our analysis, we have selected the following subsets:

Glucose as Carbon Source:

Samples: **SRR5206993, SRR5206994, SRR5206995, SRR5206996, SRR5206997**

Description: E. coli cultures grown with glucose as the primary carbon source.

Glycerol as Carbon Source:

Samples:** SRR5207019, SRR5207020, SRR5207021, SRR5207022, SRR5207023**

Description: E. coli cultures grown with glycerol as the primary carbon source.

These specific samples enable a focused comparison of gene expression profiles between E. coli grown on glucose versus glycerol, shedding light on metabolic adaptations to different carbon sources.
**
In this project, I have established an Anaconda environment to facilitate the RNA-Seq analysis pipeline, incorporating essential tools such as FastQC, Trimmomatic, Bowtie2, and RSEM.**


## Pipeline Workflow

1. **Quality Control**
   
   └── *FastQC*: Assess the quality of raw sequencing reads.
   
       └── Output: Quality reports identifying potential issues (e.g., low-quality scores, adapter contamination).

3. **Read Trimming**
   
   └── *Trimmomatic*: Remove low-quality bases and adapter sequences from reads.
   
       └── Input: Raw sequencing reads.
   
       └── Output: Cleaned reads suitable for alignment.

5. **Alignment**
   
   └── *Bowtie2*: Align cleaned reads to the *E. coli* reference genome.


       └── Input: Cleaned reads.

       └── Output: Alignment files (BAM format) indicating read positions on the genome.

7. **Quantification**
   
   └── *RSEM*: Estimate gene and isoform expression levels from alignment data.


       └── Input: Alignment files.

      └── Output: Gene and isoform expression matrices.

9. **Differential Expression Analysis**
    
   └── *DESeq2*: Identify differentially expressed genes between experimental conditions.


      └── Input: Expression matrices.

      └── Output: Lists of genes with significant expression changes, including statistical metrics.

11. **Visualization and Interpretation**
    
   └── Generate plots (e.g., heatmaps, MA plots) and perform pathway analysis to interpret biological significance.
       
       └── Input: Differential expression results.
       
       └── Output: Visual representations and functional insights into gene expression changes.

