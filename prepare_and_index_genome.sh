#!/bin/bash

# Define the genome file path
genome="genome/sbicolor_btx623.dna.toplevel.fa"

# BWA
bwa index ${genome}

# Picard
picard CreateSequenceDictionary -R ${genome}

# Samtools
samtools faidx ${genome}

