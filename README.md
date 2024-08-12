# Variant call pipeline using GATK4

Here we represent a workflow for variant call using GATK4. This includes alignment, base qual-score recaliberation, variant call, filtering SNPs-Indes and predict their effects.

# Getting started
1. Download the pipeline as follows
   ```
   git clone https://github.com/nikhilshinde0909/Variant_call_pipeline.git
   ```
   
2. Install needed softwares by creating conda environment and activate it
   ```
   mamba env create -f environment.yml
   ```
   
3. Activate conda environment
   ```
   conda activate variant_env
   ```
   
4. Obtain reference genome and resquencing reads in fastq or fastq.gz format
   
5. Add path for reference genome in file named "prepare_and_index_genome.sh" and run following command
   ```
   bash prepare_and_index_genome.sh
   ``` 
   this will index reference genome for BWA and Picard
   
6. Add details for reference genome, raw fastqs and sample prefix to file named "Run_variant_call.sh"
   
7. Now this workflow is ready for run using bash
   ```
   bash Run_variant_call.sh
   ```
     
8. We have exexcuted this workflow on sorghum so we used "Sorghum_bicolor" as snpeff database. This parameter will vary as per users interest.

# Thanks for using our variant call pipeline
