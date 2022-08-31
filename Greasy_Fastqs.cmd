#!/bin/bash
#SBATCH --job-name="FastQProcessing"
#SBATCH --workdir=/gpfs/projects/bsc08/shared_projects/Vascular_Disease/cicBIOgune/AC-41_TotalRNAseq/FASTQs
#SBATCH --output=FastQProcessing_%j.out
#SBATCH --error=FastQProcessing_%j.err
#SBATCH --ntasks=37
#SBATCH --cpus-per-task=48
#SBATCH --time=48:00:00
#SBATCH --qos=bsc_ls

module load star/2.7.0f
/apps/GREASY/latest/INTEL/IMPI/bin/greasy ../pair_ended/Run.txt
