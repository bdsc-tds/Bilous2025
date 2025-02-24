#!/bin/bash -l
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000
#SBATCH --time=24:00:00
#SBATCH --output=workflow/log/run.out

rm -rf .snakemake/slurm_logs
rm -rf .snakemake/tmp.*
mamba activate sm
snakemake --sdm=conda --executor slurm --jobs 500 --keep-going \
          --rerun-incomplete --rerun-triggers mtime  \
          --workflow-profile profiles/slurm/ "$@"