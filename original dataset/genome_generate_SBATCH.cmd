#!/bin/bash
#SBATCH --job-name="genomeGenerate"
#SBATCH --workdir=.
#SBATCH --output=ggenerate_%j.out
#SBATCH --error=ggenerate_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=12:00:00
#SBATCH --qos=bsc_ls

module load star/2.7.0f
STAR --runMode genomeGenerate --runThreadN 47 --genomeDir genome/ --genomeFastaFiles grch38/ncbi-genomes-2022-05-09/grch38 --sjdbGTFfile grch38_gtf/ncbi-genomes-2022-05-09/grch38_gtf --sjdbOverhang 150
