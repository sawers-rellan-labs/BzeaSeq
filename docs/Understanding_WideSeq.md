# Understanding WideSeq Read Count Collection Workflow

This document explains the reasoning behind each step in the WideSeq read count collection process, which is designed to efficiently collect allelic counts at teosinte-specific variant positions.

## Workflow Overview

The WideSeq analysis pipeline uses GATK's CollectAllelicCounts tool to gather reference and alternate allele read counts at specific genomic positions. This approach provides a more efficient alternative to traditional variant calling, as it directly counts reads at predetermined positions without the computational overhead of full genotype calling.

## Step-by-Step Explanation

### 1. Creating Sequence Dictionary

```bash
picard CreateSequenceDictionary --REFERENCE="../ref/Zm-B73-REFERENCE-NAM-5.0.fa"  --OUTPUT="../ref/Zm-B73-REFERENCE-NAM-5.0.dict"
```

**Reasoning**: GATK tools require a sequence dictionary (.dict) file alongside the reference genome. This file contains metadata about the reference genome's contigs (chromosomes) including their names and lengths. Creating this dictionary is a one-time preprocessing step necessary for all subsequent GATK operations.

### 2. Adding Read Groups to BAM Files

```bash
picard AddOrReplaceReadGroups \
  -I /path/to/sample.bam \
  -O bam/sample_sorted_rg.bam \
  --RGPL illumina \
  --RGLB Lane1 \
  --RGPU NONE \
  --RGSM sample_name
```

**Reasoning**: GATK tools require read group information in BAM files. Read groups identify which reads belong to the same sequencing run and sample, which is essential for proper variant calling and read count collection. Without read groups, GATK will refuse to process the BAM files. The specific read group values (RGPL, RGLB, etc.) provide metadata about the sequencing platform and sample identity.

### 3. Indexing the Read Group BAM File

```bash
samtools index bam/sample_sorted_rg.bam
```

**Reasoning**: Indexed BAM files enable random access to specific genomic regions, dramatically improving performance for region-specific operations. GATK requires BAM files to be both sorted and indexed. This indexing step creates a .bai file that allows GATK to efficiently locate reads from specific chromosomal regions without scanning the entire file.

### 4. Indexing the VCF File

```bash
gatk IndexFeatureFile -I ./wideseq_ref/wideseq_chr10.vcf.gz
```

**Reasoning**: Similar to BAM indexing, VCF files need to be indexed for random access to positions. This creates a tabix index (.tbi file) that allows GATK and other tools to quickly locate variants in specific regions without parsing the entire VCF file. The index is essential for the next step that uses the VCF as a region selector.

### 5. Selecting Variants as Interval Targets

```bash
gatk SelectVariants \
  -R ../ref/Zm-B73-REFERENCE-NAM-5.0.fa \
  -V /path/to/sample_filtered.vcf.gz \
  -L ./wideseq_ref/wideseq_chr10.vcf.gz \
  -O ./bzea/sample_wideseq_selected_chr10.vcf
```

**Reasoning**: This step creates an intersection between two VCF files:
1. The sample's previously called variants (`sample_filtered.vcf.gz`)
2. The teosinte reference panel variants (`wideseq_chr10.vcf.gz`)

This intersection is critical as it focuses the allelic count collection on positions that are:
- Known to be variable in teosinte populations (from wideseq_ref)
- Also identified as variant in the sample being analyzed

This approach greatly reduces the number of positions that need to be processed, increasing efficiency while maintaining focus on biologically relevant sites. The resulting VCF serves as the position list (intervals) for the next step.

### 6. Collecting Allelic Counts

```bash
gatk CollectAllelicCounts \
    -I bam/sample_sorted_rg.bam \
    -R ../ref/Zm-B73-REFERENCE-NAM-5.0.fa \
    -L ./bzea/sample_wideseq_selected_chr10.vcf \
    -O ./bzea/sample.allelicCounts.tsv
```

**Reasoning**: This is the core step that examines the BAM file and counts how many reads support the reference allele and how many support the alternate allele at each position defined in the intersection VCF. The advantages of this approach include:

1. **Efficiency**: Only specified positions are examined, dramatically reducing computation time
2. **Accuracy**: Raw read counts provide more precise information about allele frequencies than genotype calls
3. **Low-coverage handling**: Even positions with few reads can provide useful information about allele balance
4. **Consistency**: The same positions are analyzed across all samples for better comparability

The output is a tab-separated file with columns for chromosome, position, reference nucleotide, alternate nucleotide, reference count, alternate count, and total count. This format is ideal for downstream analysis in R with data.table, as it provides the raw data needed to calculate allele frequencies and bin-level statistics.

## Conclusion

This workflow efficiently extracts allelic count information at teosinte reference panel positions, focusing only on sites that are also variable in the sample being analyzed. By collecting raw read counts rather than making genotype calls, the approach provides more granular data for downstream ancestry analysis while significantly reducing computational requirements.

The output allelic counts serve as the foundation for bin-level genotype calling (REF/HET/ALT), Jaccard similarity calculations with teosinte accessions, and ultimately the phylogenetic analysis and introgression visualization in the complete WideSeq pipeline.