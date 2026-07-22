#!/bin/sh

#SBATCH --job-name="blast"
#SBATCH --partition=bik_p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=12G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ad14556@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -e blast.err-%N
#SBATCH -o blast.out-%N


module load BLAST+

INPUT=/scratch/ad14556/nematode-microbiome-final/01-analysis/09-barrnap

for FILE in ${INPUT}/*; do
  SAMPLE=$(basename ${FILE} .fasta)
  blastn -query ${INPUT}/${SAMPLE}.fasta -db nt -remote -outfmt "6 qseqid sseqid staxids pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -out ${INPUT}/${SAMPLE}-blast.txt -max_target_seqs 1
done

