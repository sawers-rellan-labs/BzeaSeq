# BzeaSeq: Teosinte HapMap Integration and WideSeq Analysis

This repository contains pipelines for integrating wild maize relatives into the HapMap3 SNP database and performing ancestry segment calling using the WideSeq approach.

## Overview

This project consists of two main pipelines:

1. **HapMap Expansion Pipeline**: Simulates short reads from wild maize relative assemblies, aligns them to the B73 reference, calls SNPs, and integrates them with HapMap3.

2. **WideSeq Analysis Pipeline**: Processes WideSeq data to identify ancestry segments by aligning to B73, calling SNPs, comparing to HapMap3, and calculating bin frequencies and haplotype similarities.

Both pipelines are optimized for high-performance computing environments using LSF job scheduling.

## Directory Structure

- `scripts/`: Contains all pipeline scripts
  - `hapmap_expansion/`: Scripts for HapMap3 expansion pipeline
  - `wideseq/`: Scripts for WideSeq ancestry segment calling
  - `utilities/`: Helper scripts for monitoring and job management
- `data/`: Input data
  - `reference/`: B73 reference genome
  - `wild_relatives/`: Wild relative genome assemblies
  - `wideseq_fastq/`: WideSeq sequencing data
  - `hapmap/`: HapMap3 VCF files
- `results/`: Pipeline outputs
- `logs/`: Log files from pipeline runs
- `envs/`: Conda environment files
- `docs/`: Documentation

## Workflow Diagram

```mermaid
graph TD
    subgraph "HapMap Expansion Pipeline"
        A1[Wild Maize Assemblies] --> A2[Read Simulation]
        A2 --> A3["BWA Mapping (by chromosome)"]
        A3 --> A4["Variant Calling (by chromosome)"]
        A4 --> A5["SNP Filtering (by chromosome)"]
        A5 --> A6["Intersect with HapMap3 (by chromosome)"]
        A6 --> A7[Merge Chromosomes]
        A7 --> A8[Extended HapMap Database]
    end

    subgraph "WideSeq Analysis Pipeline"
        B1[WideSeq Raw FASTQ] --> B2[Trimmomatic]
        B2 --> B3["BWA Mapping (by chromosome)"]
        B3 --> B4["BCFtools mpileup (by chromosome)"]
        B4 --> B5["Filter with Extended HapMap (by chromosome)"]
        
        A8 -.-> B5
        
        B5 --> B6["Bin Frequency Calculation (by chromosome)"]
        B5 --> B7["Extract Non-reference SNPs (by chromosome)"]
        B7 --> B8["HapMap Panel Subsetting (by chromosome)"]
        B8 --> B9["Jaccard Index Calculation (by chromosome)"]
        
        B6 --> B10[Merge Bin Frequencies]
        B9 --> B11[Merge Jaccard Results]
        
        B10 --> B12[Genome-wide Visualization]
        B11 --> B13[Top Matches Visualization]
        
        B12 --> B14[Final Report]
        B13 --> B14
    end
    
    subgraph "Resource Management"
        C1[LSF Job Submission]
        C2[Resource Monitoring]
        C3[Conda Environment]
        
        C1 -.-> A2
        C1 -.-> A3
        C1 -.-> A4
        C1 -.-> A5
        C1 -.-> A6
        C1 -.-> A7
        
        C1 -.-> B2
        C1 -.-> B3
        C1 -.-> B4
        C1 -.-> B5
        C1 -.-> B6
        C1 -.-> B7
        C1 -.-> B8
        C1 -.-> B9
        C1 -.-> B10
        C1 -.-> B11
        C1 -.-> B12
        C1 -.-> B13
        C1 -.-> B14
        
        C2 -.-> C1
        C3 -.-> C1
    end
    
    style A8 fill:#f9f,stroke:#333,stroke-width:2px
    style B14 fill:#f9f,stroke:#333,stroke-width:2px