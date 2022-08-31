#!/bin/bash
#SBATCH --job-name="FastQAdapterRm"
#SBATCH --workdir=/gpfs/projects/bsc08/shared_projects/Vascular_Disease/cicBIOgune/AC-41_TotalRNAseq/FASTQs
#SBATCH --output=AdapterRm%j.out
#SBATCH --error=AdapterRm%j.err
#SBATCH --ntasks=37
#SBATCH --cpus-per-task=48
#SBATCH --time=48:00:00
#SBATCH --qos=bsc_ls

/apps/GREASY/latest/INTEL/IMPI/bin/greasy ../pair_ended/trim.txt
