#!/bin/bash

# Define Parameters
Reference="genome/sbicolor_btx623.dna.toplevel.fa"
Read_R1="raw_data/PI144134_1.fastq.gz"
Read_R2="raw_data/PI144134_2.fastq.gz"
Prefix="PI144134"
threads="8"

# Start analysis
mkdir ${Prefix}

# Define R option
bwa_r_opts=$(echo '"@RG\tID:'${Prefix}'\tLB:'${Prefix}'\tPL:illumina\tPM:nextseq\tSM:'${Prefix}'"')

# Alignment
bwa mem -K 100000000 -Y -v 3 -t ${threads} -R ${bwa_r_opts} ${Reference} ${Read_R1} ${Read_R2} > ${Prefix}_aligned_reads.sam

# Mark duplicates
gatk MarkDuplicatesSpark -I ${Prefix}_aligned_reads.sam -M ${Prefix}_dedup_metrics.txt -O ${Prefix}_sorted_dedup_reads.bam

# Index bam
samtools index ${Prefix}_sorted_dedup_reads.bam

# Picard collect matrices
picard CollectAlignmentSummaryMetrics R=${Reference} I=${Prefix}_sorted_dedup_reads.bam O=${Prefix}_alignment_metrics.txt

# Picard histogram
picard CollectInsertSizeMetrics INPUT=${Prefix}_sorted_dedup_reads.bam OUTPUT=${Prefix}_insert_metrics.txt HISTOGRAM_FILE=${Prefix}_insert_size_histogram.pdf

# Call haplotype
gatk HaplotypeCaller --input ${Prefix}_sorted_dedup_reads.bam --output ${Prefix}_raw_variants.vcf --reference ${Reference}

# Select raw SNPs
gatk SelectVariants --variant ${Prefix}_raw_variants.vcf --reference ${Reference} -select-type SNP --output ${Prefix}_raw_snps.vcf

#select raw IDELs
gatk SelectVariants --variant ${Prefix}_raw_variants.vcf --reference ${Reference} -select-type INDEL --output ${Prefix}_raw_indels.vcf

# Filter SNPs
gatk VariantFiltration --reference ${Reference} --variant ${Prefix}_raw_snps.vcf --output ${Prefix}_filtered_snps.vcf --filter-name "QD_filter" -filter "QD < 2.0" --filter-name "FS_filter" -filter "FS > 60.0" --filter-name "MQ_filter" -filter "MQ < 40.0" --filter-name "SOR_filter" -filter "SOR > 4.0" --filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" --filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

# Filter IDELs
gatk VariantFiltration --reference ${Reference} --variant ${Prefix}_raw_indels.vcf --output ${Prefix}_filtered_indels.vcf --filter-name "QD_filter" -filter "QD < 2.0" --filter-name "FS_filter" -filter "FS > 200.0" --filter-name "SOR_filter" -filter "SOR > 10.0"

# Exclude filtered SNPs
gatk SelectVariants --exclude-filtered --variant ${Prefix}_filtered_snps.vcf --output ${Prefix}_bqsr_snps.vcf

# Exclude filtered IDELs
gatk SelectVariants --exclude-filtered --variant ${Prefix}_filtered_indels.vcf --output ${Prefix}_bqsr_indels.vcf

# base recaliberation
gatk BaseRecalibrator --reference ${Reference} --input ${Prefix}_sorted_dedup_reads.bam --known-sites ${Prefix}_bqsr_snps.vcf --known-sites ${Prefix}_bqsr_indels.vcf --output ${Prefix}_recal_data.table

# Apply BSQR
gatk ApplyBQSR --reference ${Reference} --input ${Prefix}_sorted_dedup_reads.bam --bqsr-recal-file ${Prefix}_recal_data.table --output ${Prefix}_recal_reads.bam

# index BSQR Bam
samtools index ${Prefix}_recal_reads.bam

# base recaliberation post recaliberation
gatk BaseRecalibrator --reference ${Reference} --input ${Prefix}_recal_reads.bam --known-sites ${Prefix}_bqsr_snps.vcf --known-sites ${Prefix}_bqsr_indels.vcf --output ${Prefix}_post_recal_data.table

# Analyze Covariates
gatk AnalyzeCovariates -before ${Prefix}_recal_data.table -after ${Prefix}_post_recal_data.table -plots ${Prefix}_recalibration_plots.pdf

# Call Variants
gatk HaplotypeCaller --reference ${Reference} --input ${Prefix}_recal_reads.bam --output ${Prefix}_raw_variants_recal.vcf

# Extract SNPs
gatk SelectVariants --reference ${Reference} --variant ${Prefix}_raw_variants_recal.vcf -select-type SNP --output ${Prefix}_raw_snps_recal.vcf

# Extract INDELS
gatk SelectVariants --reference ${Reference} --variant ${Prefix}_raw_variants_recal.vcf -select-type INDEL --output ${Prefix}_raw_indels_recal.vcf

# Filter SNPs
gatk VariantFiltration --reference ${Reference} --variant ${Prefix}_raw_snps_recal.vcf --output ${Prefix}/${Prefix}_filtered_snps_final.vcf --filter-name "QD_filter" -filter "QD < 2.0" --filter-name "FS_filter" -filter "FS > 60.0" --filter-name "MQ_filter" -filter "MQ < 40.0" --filter-name "SOR_filter" -filter "SOR > 4.0" --filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" --filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

# Filter INDELs
gatk VariantFiltration --reference ${Reference} --variant ${Prefix}_raw_indels_recal.vcf --output ${Prefix}/${Prefix}_filtered_indels_final.vcf --filter-name "QD_filter" -filter "QD < 2.0" --filter-name "FS_filter" -filter "FS > 200.0" --filter-name "SOR_filter" -filter "SOR > 10.0"

# Compress vcf
bgzip ${Prefix}/${Prefix}_filtered_snps_final.vcf

# Index vcf
bcftools index ${Prefix}/${Prefix}_filtered_snps_final.vcf.gz 

# Annotate SNPs and Predict Effects
snpEff Sorghum_bicolor ${Prefix}/${Prefix}_filtered_snps_final.vcf.gz > ${Prefix}/${Prefix}_filtered_snps_Anno.vcf
snpEff eff Sorghum_bicolor ${Prefix}/${Prefix}_filtered_snps_final.vcf.gz -htmlStats ${Prefix}/${Prefix}_summary.html > ${Prefix}/${Prefix}_filtered_snps_Anno_eff.vcf

# Pipeline status
file1=${Prefix}/${Prefix}'_filtered_snps_Anno.vcf'
file2=${Prefix}/${Prefix}'_filtered_snps_Anno_eff.vcf'

if [ -f "$file1" ] && [ -f "$file2" ]; then
echo '"'${Prefix}'_filtered_snps_Anno.vcf" and "'${Prefix}'_filtered_snps_Anno_eff.vcf" are exists'
echo 'All done Good...'
else
echo "SnpEff step failed please check inputs"
fi
