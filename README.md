# Variant call pipeline using GATK4

Here we represent a workflow for variant call using GATK4. This includes alignment, base qual-score recaliberation, variant call, filtering SNPs-Indes and predict their effects.

# Getting started
1. Install needed softwares by creating conda environment
   ```
   mamba env create -f environment.yml
   ```
3. Obtain reference genome and resquencing reads in fastq or fastq.gz format
4. Add path for reference genome in file named "prepare_and_index_genome.sh" and run following command
   ```
   bash prepare_and_index_genome.sh
   ``` 
   this will index reference genome for BWA and Picard
6. Add details for reference genome, raw fastqs and sample prefix to file named "Run_variant_call.sh" 
7. Now this workflow is ready for run using bash
   ```
   bash Run_variant_call.sh
   ```
9. We have exexcuted this workflow on sorghum so we used "Sorghum_bicolor" as snpeff database. This parameter will vary users need.

# Thanks for using our variant call pipeline
