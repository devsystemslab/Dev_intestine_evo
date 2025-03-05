#!/bin/bash

#SBATCH -n 5
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --job-name=${1}_${2}_chrombpnet
#SBATCH --output=slurm_logs/%x_grid_%A_%a.out  # Output file
#SBATCH --gpus=1
#SBATCH --gres=gpumem:20G

model_path="$1"
regions="$2"
genome_fasta="$3"
genome_sizes="$4"
output_prefix="$5"

# load software modules
#module load eth_proxy gcc/8.2.0 cuda/11.2.2 cudnn/8.8.1.3
source $HOME/.bashrc

# activate a mamba environment
echo 'Activating the environment'
mamba activate chrombpnet2
export LD_LIBRARY_PATH=/cluster/project/treutlein/jjans/software/miniforge3/envs/cuda11_env/lib:$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

chrombpnet contribs_bw \
        -m $model_path \
        -r $regions \
        -g $genome_fasta \
        -c $genome_sizes \
        -op $output_prefix
