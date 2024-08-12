#!/bin/bash

#SBATCH --account=bgmp                    
#SBATCH --partition=bgmp                  
#SBATCH --cpus-per-task=4                 
#SBATCH --mem=16GB                          

conda activate bgmp_py312
/usr/bin/time -v ./demux.py \
       -f1 "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz" \
       -f2 "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz" \
       -f3 "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz" \
       -f4 "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz" \