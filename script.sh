#!/bin/bash

#SBATCH --job-name="nematode-microbiome"
#SBATCH --partition=bik_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#SBATCH --time=07-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e nematode-microbiome.err-%N
#SBATCH -o nematode-microbiome.out-%N

#Path Variables
INPUT=/home/ad14556/snakemake-pipelines/metagenomicsPipeline

#Activate QIIME2
module load Miniconda3
source activate /home/ad14556/conda-envs/envs/snakemake

snakemake --jobs 6 --use-conda --cluster-config config/cluster.yaml --keep-going \
    --cluster "sbatch --parsable --qos=unlim --partition={cluster.partition} --job-name={cluster.name} --mem={cluster.mem}gb --time={cluster.time} --nodes={cluster.nodes} --cpus-per-task={cluster.cpu} --ntasks={cluster.tasks} --mail-user={cluster.mail_user} --mail-type={cluster.mail_type} --error={cluster.name}.err --output={cluster.name}.out"
