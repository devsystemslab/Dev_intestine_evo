#!/bin/bash

#SBATCH -n 5
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --job-name={1}_{2}_chrombpnet
#SBATCH --output=slurm_logs/{1}_{2}_chrombpnet_grid_%A_%a.out
#SBATCH --gpus=1
#SBATCH --gres=gpumem:20G

bam_file="$1"
peaks="$2"
negative_peaks="$3"
fold_json="$4"
output_model="$5"
genome_fasta="$6"
genome_sizes="$7"
bias_model="$8"

# load software modules
#module load eth_proxy gcc/8.2.0 cuda/11.2.2 cudnn/8.8.1.3
source $HOME/.bashrc

# activate a mamba environment
echo 'Activating the environment'
mamba activate chrombpnet2
export LD_LIBRARY_PATH=/cluster/project/treutlein/jjans/software/miniforge3/envs/cuda11_env/lib:$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

chrombpnet pipeline \
        -ibam $bam_file \
        -d "ATAC" \
        -g $genome_fasta \
        -c $genome_sizes \
        -p $peaks \
        -n $negative_peaks \
        -fl $fold_json \
        -b $bias_model \
        -o $output_model