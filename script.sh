#!/bin/bash

#SBATCH --job-name="nematode-microbiome"
#SBATCH --partition=bik_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#SBATCH --time=05-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e nematode-microbiome.err-%N
#SBATCH -o nematode-microbiome.out-%N

#Path Variables
INPUT=/home/ad14556/snakemake-pipelines/metagenomicsPipeline

#Activate QIIME2
module load Miniforge3/24.1.2-0
source activate /home/ad14556/conda-envs/envs/snakemake

#export TMPDIR=/scratch/ad14556/tmp
#tmpdir=/scratch/ad14556/tmp/

snakemake --ignore-incomplete --jobs 15 --until metagenome_assembly --use-conda --cluster-config config/cluster.yaml --keep-going --latency-wait 5 \
    --cluster "sbatch --parsable --partition={cluster.partition} --job-name={cluster.name} --mem={cluster.mem}G --time={cluster.time} --nodes={cluster.nodes} --cpus-per-task={cluster.cpu} --ntasks={cluster.tasks} --mail-user={cluster.mail_user} --mail-type={cluster.mail_type}"
