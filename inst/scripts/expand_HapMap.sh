#!/bin/bash

#BSUB -J hapmap_expansion
#BSUB -n 4
#BSUB -q normal
#BSUB -R "rusage[mem=16GB]"
#BSUB -o hapmap_expansion_%J.out
#BSUB -e hapmap_expansion_%J.err

# Set variables
FASTA_DIR="/path/to/wild_maize_fasta"
REFERENCE="/path/to/B73_v4_reference.fasta"
OUTPUT_DIR="/path/to/output"
HAPMAP_VCF="/path/to/HapMap3.vcf.gz"
THREADS=4

# Create output directories
mkdir -p ${OUTPUT_DIR}/{reads,mapping,variants,hapmap_filtered,merged}

# List of chromosomes in maize
CHROMOSOMES=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10")

# 1. Simulate reads from wild maize assemblies (once per assembly)
for fasta in ${FASTA_DIR}/*.fasta; do
    base=$(basename $fasta .fasta)
    
    # Submit read simulation job
    bsub -J "sim_${base}" \
         -n 4 \
         -R "rusage[mem=16GB]" \
         -o ${OUTPUT_DIR}/logs/sim_${base}.out \
         -e ${OUTPUT_DIR}/logs/sim_${base}.err \
         "art_illumina -ss HS25 -i ${fasta} -p -l 150 -f 30 -o ${OUTPUT_DIR}/reads/${base}_sim"
done

# 2. Index reference genome
bsub -J "index_ref" \
     -n 4 \
     -R "rusage[mem=16GB]" \
     -o ${OUTPUT_DIR}/logs/index_ref.out \
     -e ${OUTPUT_DIR}/logs/index_ref.err \
     "bwa index ${REFERENCE}"

# 3. Wait for read simulation to complete, then submit mapping jobs per chromosome
for fasta in ${FASTA_DIR}/*.fasta; do
    base=$(basename $fasta .fasta)
    
    # Submit a job array for each chromosome
    bsub -J "map_${base}[1-10]" \
         -w "done(sim_${base})" \
         -n 4 \
         -R "rusage[mem=16GB]" \
         -o ${OUTPUT_DIR}/logs/map_${base}_%I.out \
         -e ${OUTPUT_DIR}/logs/map_${base}_%I.err \
         "
         # Get chromosome number from job array index
         idx=\$((LSB_JOBINDEX - 1))
         chrom=${CHROMOSOMES[\$idx]}
         
         # Create chromosome-specific reference
         samtools faidx ${REFERENCE} \${chrom} > ${OUTPUT_DIR}/mapping/\${chrom}.fasta
         bwa index ${OUTPUT_DIR}/mapping/\${chrom}.fasta
         
         # Map reads for this chromosome only
         bwa mem -t ${THREADS} ${OUTPUT_DIR}/mapping/\${chrom}.fasta \
             ${OUTPUT_DIR}/reads/${base}_sim1.fq ${OUTPUT_DIR}/reads/${base}_sim2.fq | \
             samtools view -b -o ${OUTPUT_DIR}/mapping/${base}.\${chrom}.raw.bam -
         
         # Sort and mark duplicates
         samtools sort -o ${OUTPUT_DIR}/mapping/${base}.\${chrom}.sorted.bam ${OUTPUT_DIR}/mapping/${base}.\${chrom}.raw.bam
         samtools markdup ${OUTPUT_DIR}/mapping/${base}.\${chrom}.sorted.bam ${OUTPUT_DIR}/mapping/${base}.\${chrom}.bam
         samtools index ${OUTPUT_DIR}/mapping/${base}.\${chrom}.bam
         
         # Clean up intermediate files
         rm ${OUTPUT_DIR}/mapping/${base}.\${chrom}.raw.bam ${OUTPUT_DIR}/mapping/${base}.\${chrom}.sorted.bam
         "
done

# 4. Call variants per chromosome using GATK
for fasta in ${FASTA_DIR}/*.fasta; do
    base=$(basename $fasta .fasta)
    
    # Submit a job array for each chromosome
    bsub -J "call_${base}[1-10]" \
         -w "done(map_${base}*)" \
         -n 4 \
         -R "rusage[mem=16GB]" \
         -o ${OUTPUT_DIR}/logs/call_${base}_%I.out \
         -e ${OUTPUT_DIR}/logs/call_${base}_%I.err \
         "
         # Get chromosome number from job array index
         idx=\$((LSB_JOBINDEX - 1))
         chrom=${CHROMOSOMES[\$idx]}
         
         # Call variants with HaplotypeCaller
         gatk HaplotypeCaller \
           -R ${OUTPUT_DIR}/mapping/\${chrom}.fasta \
           -I ${OUTPUT_DIR}/mapping/${base}.\${chrom}.bam \
           -O ${OUTPUT_DIR}/variants/${base}.\${chrom}.vcf.gz \
           -ERC GVCF
         "
done

# 5. Combine GVCFs per chromosome
for fasta in ${FASTA_DIR}/*.fasta; do
    base=$(basename $fasta .fasta)
    
    bsub -J "combine[1-10]" \
         -w "done(call_*)" \
         -n 4 \
         -R "rusage[mem=16GB]" \
         -o ${OUTPUT_DIR}/logs/combine_%I.out \
         -e ${OUTPUT_DIR}/logs/combine_%I.err \
         "
         # Get chromosome number from job array index
         idx=\$((LSB_JOBINDEX - 1))
         chrom=${CHROMOSOMES[\$idx]}
         
         # Get all GVCF files for this chromosome
         gvcfs=\$(ls ${OUTPUT_DIR}/variants/*.\${chrom}.vcf.gz)
         
         # Create input arguments for CombineGVCFs
         input_args=''
         for gvcf in \$gvcfs; do
             input_args=\"\$input_args -V \$gvcf\"
         done
         
         # Combine GVCFs
         gatk CombineGVCFs \
           -R ${OUTPUT_DIR}/mapping/\${chrom}.fasta \
           \$input_args \
           -O ${OUTPUT_DIR}/variants/combined.\${chrom}.g.vcf.gz
         
         # Joint genotyping
         gatk GenotypeGVCFs \
           -R ${OUTPUT_DIR}/mapping/\${chrom}.fasta \
           -V ${OUTPUT_DIR}/variants/combined.\${chrom}.g.vcf.gz \
           -O ${OUTPUT_DIR}/variants/raw_variants.\${chrom}.vcf.gz
         
         # Extract and filter SNPs
         gatk SelectVariants \
           -R ${OUTPUT_DIR}/mapping/\${chrom}.fasta \
           -V ${OUTPUT_DIR}/variants/raw_variants.\${chrom}.vcf.gz \
           -select-type SNP \
           -O ${OUTPUT_DIR}/variants/raw_snps.\${chrom}.vcf.gz
         
         gatk VariantFiltration \
           -R ${OUTPUT_DIR}/mapping/\${chrom}.fasta \
           -V ${OUTPUT_DIR}/variants/raw_snps.\${chrom}.vcf.gz \
           -filter \"QD < 2.0\" --filter-name \"QD2\" \
           -filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \
           -filter \"SOR > 3.0\" --filter-name \"SOR3\" \
           -filter \"FS > 60.0\" --filter-name \"FS60\" \
           -filter \"MQ < 40.0\" --filter-name \"MQ40\" \
           -filter \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" \
           -filter \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\" \
           -O ${OUTPUT_DIR}/variants/filtered_snps.\${chrom}.vcf.gz
         "
done

# 6. Intersect with HapMap3 per chromosome
bsub -J "intersect[1-10]" \
     -w "done(combine*)" \
     -n 4 \
     -R "rusage[mem=16GB]" \
     -o ${OUTPUT_DIR}/logs/intersect_%I.out \
     -e ${OUTPUT_DIR}/logs/intersect_%I.err \
     "
     # Get chromosome number from job array index
     idx=\$((LSB_JOBINDEX - 1))
     chrom=${CHROMOSOMES[\$idx]}
     
     # Extract chromosome from HapMap3
     bcftools view -r \${chrom} ${HAPMAP_VCF} -Oz -o ${OUTPUT_DIR}/hapmap_filtered/hapmap3.\${chrom}.vcf.gz
     bcftools index ${OUTPUT_DIR}/hapmap_filtered/hapmap3.\${chrom}.vcf.gz
     
     # Intersect filtered SNPs with HapMap3
     bcftools isec -n =2 -w 1 ${OUTPUT_DIR}/variants/filtered_snps.\${chrom}.vcf.gz \
         ${OUTPUT_DIR}/hapmap_filtered/hapmap3.\${chrom}.vcf.gz \
         -Oz -o ${OUTPUT_DIR}/hapmap_filtered/overlapping_snps.\${chrom}.vcf.gz
     
     bcftools index ${OUTPUT_DIR}/hapmap_filtered/overlapping_snps.\${chrom}.vcf.gz
     "

# 7. Merge chromosomes to create final expanded HapMap
bsub -J "merge_final" \
     -w "done(intersect*)" \
     -n 4 \
     -R "rusage[mem=16GB]" \
     -o ${OUTPUT_DIR}/logs/merge_final.out \
     -e ${OUTPUT_DIR}/logs/merge_final.err \
     "
     # List all chromosome VCFs
     chr_vcfs=\$(ls ${OUTPUT_DIR}/hapmap_filtered/overlapping_snps.chr*.vcf.gz)
     
     # Concatenate VCFs
     bcftools concat \$chr_vcfs -Oz -o ${OUTPUT_DIR}/merged/expanded_hapmap3.vcf.gz
     bcftools index ${OUTPUT_DIR}/merged/expanded_hapmap3.vcf.gz
     "