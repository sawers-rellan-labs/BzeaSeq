# Getting Teosinte Variants from Schnable2023 Dataset

## Table of Contents
- [1. Overview](#1-overview)
- [2. Data Acquisition](#2-data-acquisition)
- [3. Data Exploration and Validation](#3-data-exploration-and-validation)
- [4. Sample Selection and Filtering](#4-sample-selection-and-filtering)
- [5. Variant Processing Pipeline](#5-variant-processing-pipeline)
- [6. Quality Control and Statistics](#6-quality-control-and-statistics)
- [7. Integration with WideSeq Pipeline](#7-integration-with-wideseq-pipeline)

## 1. Overview

This document describes the detailed process for acquiring and processing teosinte variant data from the Schnable2023 dataset for use in the BzeaSeq WideSeq analysis pipeline. The Schnable2023 study provides genome-wide variant data from both maize and teosinte samples already aligned to the B73 reference genome version 5 (Zm-B73-REFERENCE-NAM-5.0), eliminating the need for coordinate liftover that would be required with other datasets like Chen2022.

### Key Advantages of Schnable2023 Dataset
- Pre-aligned to B73v5 reference genome
- Includes comprehensive teosinte representation across subspecies
- High-quality variant calls with consistent methodology
- Chromosome-specific VCF files for efficient processing

## 2. Data Acquisition

### 2.1 Download Infrastructure

The VCF files are distributed across 10 chromosome-specific files available from the SNPVersity 2.0 repository. Download URLs and checksums are provided in the `inst/extdata/checksums.tab` file.

```bash
# Example checksums.tab format:
# filename	description	size	download_time	md5	url
schnable2023_chr1.vcf.gz	High quality filtered variants with imputation identified in Chromsome 1.	1.2 GB	2 minutes / 16 minutes	5b576f00da277a8ee00c09013b47eb92	https://ars-usda.app.box.com/s/rv1vsf5k36x5d5fu28p3bzhxykfpyybs/file/1650670571309
```

### 2.2 Download Script Implementation

```bash
#!/bin/bash
# download_schnable2023_SNPs.sh
#
# Downloads VCF files for all 10 chromosomes from the Schnable2023 dataset.
# Includes MD5 checksum verification after download.

# Configuration
OUTPUT_DIR="./schnable2023"
CHECKSUMS_FILE="checksums.tab"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Process each chromosome
for chr in {1..10}; do
    echo "Downloading chromosome ${chr} VCF file..."
    
    # Extract download URL and MD5 checksum from checksums file
    URL=$(grep "schnable2023_chr${chr}.vcf.gz" ${CHECKSUMS_FILE} | awk -F'\t' '{print $6}')
    EXPECTED_MD5=$(grep "schnable2023_chr${chr}.vcf.gz" ${CHECKSUMS_FILE} | awk -F'\t' '{print $5}')
    
    if [ -z "${URL}" ] || [ -z "${EXPECTED_MD5}" ]; then
        echo "Error: Could not find URL or MD5 for chromosome ${chr}"
        continue
    fi
    
    OUTFILE="${OUTPUT_DIR}/schnable2023_chr${chr}.vcf.gz"
    
    # Download the file
    curl -L -o ${OUTFILE} ${URL}
    
    # Verify MD5 checksum
    ACTUAL_MD5=$(md5sum ${OUTFILE} | awk '{print $1}')
    
    if [ "${ACTUAL_MD5}" = "${EXPECTED_MD5}" ]; then
        echo "Verification successful for chromosome ${chr}"
    else
        echo "Verification failed for chromosome ${chr}"
        echo "Expected: ${EXPECTED_MD5}"
        echo "Actual: ${ACTUAL_MD5}"
    fi
done

echo "Download completed."
```

### 2.3 File Naming Corrections

Some downloaded files may contain naming inconsistencies that need correction:

```bash
# Fix file naming inconsistencies using bash parameter substitution
for file in schanble2023_*.vcf.gz; do 
    mv -v "$file" "${file/schanble/schnable}"
done
```

This ensures consistent `schnable2023_` prefixes throughout the dataset.

## 3. Data Exploration and Validation

### 3.1 Basic Dataset Statistics

First, examine the dataset composition to understand sample count and variant distribution:

```bash
# Count total samples in the dataset
bcftools query -l schnable2023/schnable2023_chr10.vcf.gz | wc -l
```

**Output:**
```
744
```

```bash
# Generate a complete sample list
bcftools query -l schnable2023/schnable2023_chr10.vcf.gz > schnable2023_id.list

# Get basic statistics for one chromosome
bcftools stats schnable2023/schnable2023_chr10.vcf.gz > schnable2023_chr10.stats
more schnable2023_chr10.stats
```

**Output:**
```
SN	0	number of samples:	744
SN	0	number of records:	2054020
SN	0	number of no-ALTs:	0
SN	0	number of SNPs:	2054020
SN	0	number of MNPs:	0
SN	0	number of indels:	0
SN	0	number of others:	0
SN	0	number of multiallelic sites:	125043
SN	0	number of multiallelic SNP sites:	125043
```

### 3.2 Metadata Analysis and Processing

The metadata for the Schnable2023 dataset requires careful curation due to inconsistencies in the original supplementary materials:

```bash
# Taxonomic distribution analysis
grep -v "Z. mays" schnable2023/schnable2023_metadata.tab \
  | cut -f2 | grep -v "Species" \
  | sort | uniq -c
```

**Output:**
```
      1 New taxonomyb
     14 Teosinte (mix)
     20 Zea diploperennis
     14 Zea luxurians
      5 Zea mays subsp. huehuetenangensis
     81 Zea mays subsp. mexicana
     70 Zea mays subsp. parviglumis
     14 Zea nicaraguensis
     19 Zea perennis
```

### 3.3 Reference Sample Validation

Confirm that B73 (the reference genotype) is included in the dataset:

```bash
grep B73 schnable2023_metadata.tab
```

**Output:**
```
B73	Z. mays	Inbred line	USA	285169576	277505548	97.31	97.47	19.72	60259366	21.13	-	-	Zea mays subsp. mays	Zea mays subsp. mays (TEM)
```

## 4. Sample Selection and Filtering

### 4.1 Teosinte Sample Identification

Create comprehensive lists of samples for the WideSeq reference panel:

```bash
# Extract wild relatives (non-maize subspecies, excluding Tripsacum)
grep -v "Z. mays" schnable2023/schnable2023_metadata.tab \
  | grep -v "Tripsacum" | cut -f1 | grep -v "ID" \
  > wild_relatives_id.list

# Count teosinte samples
grep -v "subsp. mays" schnable2023/schnable2023_metadata.tab \
  | grep -v "Tripsacum" | wc -l
```

**Output:**
```
238
```

### 4.2 Reference Panel Construction

```bash
# Add B73 to the reference sample list
grep "B73" schnable2023/schnable2023_metadata.tab \
  | cut -f1 > B73_id.list

# Combine B73 and teosinte samples, excluding Tripsacum
cat B73_id.list wild_relatives_id.list \
  | grep -v tripsacum | sort | uniq > wideseq_ref_id.list

# Verify final sample count
wc -l wideseq_ref_id.list
```

**Output:**
```
239 wideseq_ref_id.list
```

### 4.3 Dataset Compatibility Analysis

Compare with other datasets for cross-validation:

```bash
# Check overlap with Chen2022 dataset (if available)
grep -w -f chen2022_id.list schnable2023_id.list | wc -l

# Verify reference samples are present
grep -w -f wideseq_ref_id.list schnable2023_id.list | wc -l
```

## 5. Variant Processing Pipeline

### 5.1 Sample Extraction and Filtering Script

The main processing script applies sample filtering and minor allele frequency thresholds:

```bash
#!/bin/bash
# get_variants_from_schnable2023.sh
#
# Extracts teosinte and reference samples from Schnable2023 VCF files
# and applies MAF filtering.

INPUT_DIR="schnable2023"
OUTPUT_DIR="wideseq_ref"
SAMPLE_LIST="wideseq_ref_id.list"

# MAF threshold - this will bias against rare alleles like huehuetenangensis
MAF_THRESHOLD="0.05"

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Process each chromosome
for chr in {1..10}; do
    echo "Processing chromosome ${chr}..."
    
    # Input and output file paths
    INPUT_VCF="${INPUT_DIR}/schnable2023_chr${chr}.vcf.gz"
    SAMPLE_FILTER_VCF="${OUTPUT_DIR}/wideseq_taxa_chr${chr}.vcf.gz"
    OUTPUT_VCF="${OUTPUT_DIR}/wideseq_chr${chr}.vcf.gz"
    
    # Extract samples and filter for biallelic SNPs only
    bcftools view -S ${SAMPLE_LIST} -m2 -M2 -v snps ${INPUT_VCF} \
      -Oz -o ${SAMPLE_FILTER_VCF}
    
    # Apply MAF filtering   
    bcftools view --min-af ${MAF_THRESHOLD}:minor -v snps ${SAMPLE_FILTER_VCF} \
      -Oz -o ${OUTPUT_VCF}
    
    # Index the output file for both bcftools and GATK
    bcftools index ${OUTPUT_VCF}
    gatk IndexFeatureFile -I ${OUTPUT_VCF}
    
    echo "Completed chromosome ${chr}"
done

echo "All chromosomes processed."
```

### 5.2 Processing Steps Explanation

The pipeline performs several critical filtering steps:

1. **Sample Selection** (`-S ${SAMPLE_LIST}`): Extracts only teosinte and B73 reference samples
2. **Biallelic Filtering** (`-m2 -M2`): Ensures only sites with exactly 2 alleles
3. **SNP-only Filtering** (`-v snps`): Removes indels and complex variants
4. **MAF Filtering** (`--min-af 0.05:minor`): Requires minor allele frequency ≥ 5%

**Important Note on MAF Filtering:** The 5% MAF threshold may bias against rare teosinte subspecies like *Z. mays* subsp. *huehuetenangensis* which has only 5 accessions in the dataset. This trade-off ensures sufficient statistical power for ancestry analysis while potentially missing subspecies-specific variants.

## 6. Quality Control and Statistics

### 6.1 Variant Count Analysis

Calculate comprehensive statistics for the filtered variant set:

```bash
#!/bin/bash
# calculate_variant_stats.sh
#
# Calculates statistics for the filtered teosinte variants.

INPUT_DIR="./wideseq_ref"
OUTPUT_DIR="./stats"

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Calculate total variants across all chromosomes
total_variants=0

# Process each chromosome
for chr in {1..10}; do
    echo "Calculating statistics for chromosome ${chr}..."
    
    # Input file
    INPUT_VCF="${INPUT_DIR}/wideseq_chr${chr}.vcf.gz"
    
    # Calculate statistics
    bcftools stats ${INPUT_VCF} > ${OUTPUT_DIR}/wideseq_chr${chr}.stats
    
    # Count variants
    variants=$(grep "number of SNPs:" ${OUTPUT_DIR}/wideseq_chr${chr}.stats | awk '{print $6}')
    total_variants=$((total_variants + variants))
    
    echo "Chromosome ${chr}: ${variants} variants"
done

echo "Total variants across all chromosomes: ${total_variants}"

# Calculate average variant density per 100kb
genome_size=2131846805  # Zm-B73-REFERENCE-NAM-5.0 nuclear chromosome length in bp
density_per_100kb=$(echo "scale=2; ${total_variants} * 100000 / ${genome_size}" | bc)

echo "Average variant density: ${density_per_100kb} variants per 100kb"
```

**Expected Output:**
```
Calculating statistics for chromosome 1...
Chromosome 1: 4006441 variants
Calculating statistics for chromosome 2...
Chromosome 2: 3100759 variants
Calculating statistics for chromosome 3...
Chromosome 3: 3128544 variants
Calculating statistics for chromosome 4...
Chromosome 4: 3529519 variants
Calculating statistics for chromosome 5...
Chromosome 5: 2878472 variants
Calculating statistics for chromosome 6...
Chromosome 6: 2143655 variants
Calculating statistics for chromosome 7...
Chromosome 7: 2297952 variants
Calculating statistics for chromosome 8...
Chromosome 8: 2283288 variants
Calculating statistics for chromosome 9...
Chromosome 9: 2181127 variants
Calculating statistics for chromosome 10...
Chromosome 10: 2054020 variants
Total variants across all chromosomes: 27603777
Average variant density: 1294.82 variants per 100kb
```

### 6.2 Quality Metrics

The final teosinte reference variant set provides:
- **27.6 million SNPs** across the nuclear genome
- **~1,295 variants per 100kb** average density
- **Comprehensive teosinte representation** across 8 subspecies/species
- **High-quality variants** with MAF ≥ 5% in the reference panel

## 7. Integration with WideSeq Pipeline

### 7.1 File Organization

The processed variants are organized for seamless integration with the WideSeq analysis pipeline:

```
wideseq_ref/
├── wideseq_chr1.vcf.gz
├── wideseq_chr1.vcf.gz.tbi
├── wideseq_chr1.vcf.gz.idx
├── wideseq_chr2.vcf.gz
├── wideseq_chr2.vcf.gz.tbi
├── wideseq_chr2.vcf.gz.idx
...
└── wideseq_chr10.vcf.gz
    ├── wideseq_chr10.vcf.gz.tbi
    └── wideseq_chr10.vcf.gz.idx
```

### 7.2 Index Files

Both bcftools (.tbi) and GATK (.idx) index files are created to ensure compatibility with downstream tools:

```bash
# bcftools index for general VCF operations
bcftools index ${OUTPUT_VCF}

# GATK index for CollectAllelicCounts and other GATK tools
gatk IndexFeatureFile -I ${OUTPUT_VCF}
```

### 7.3 Validation for WideSeq Pipeline

Verify the processed files are ready for the WideSeq pipeline:

```bash
# Test GATK compatibility
gatk ValidateVariants \
  -R ../ref/Zm-B73-REFERENCE-NAM-5.0.fa \
  -V wideseq_ref/wideseq_chr1.vcf.gz \
  --validation-type-to-exclude ALL

# Test bcftools compatibility
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
  wideseq_ref/wideseq_chr1.vcf.gz | head -5
```

### 7.4 Pipeline Integration Notes

The teosinte reference variant set integrates with the WideSeq pipeline through:

1. **Interval List Creation**: VCF files serve as interval lists for GATK CollectAllelicCounts
2. **Position-specific Analysis**: Only positions present in the teosinte panel are analyzed
3. **Ancestry Classification**: Variants provide the reference framework for ancestry calling
4. **Phylogenetic Analysis**: Teosinte samples enable ancestry assignment through similarity calculations

## References

1. Schnable, J.C., et al. (2023). Exploring the pan-genome of the *Zea* genus through genome-wide association studies. The Plant Journal, 116(1), 225-238. https://doi.org/10.1111/tpj.16123

2. Hufford, M.B., Seetharam, A.S., Woodhouse, M.R., et al. (2021). De novo assembly, annotation, and comparative analysis of 26 diverse maize genomes. Science, 373(6555), 655-662. https://doi.org/10.1126/science.abg5289

3. Andorf, C.M., Ross-Ibarra, J., Seetharam, A.S., et al. (2025). A unified VCF dataset from nearly 1,500 diverse maize accessions and resources to explore the genomic landscape of maize. G3 Genes|Genomes|Genetics, 15(2), jkae281. https://doi.org/10.1093/g3journal/jkae281
