rule megahit:
    input:
        multiqc=scratch_dir + "01-analysis/05-trimmed-multiqc/multiqc.html",
        R1=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz", sample=SAMPLES),
        R2=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz", sample=SAMPLES),
        unpaired_R1=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_unpaired.fastq.gz", sample=SAMPLES),
        unpaired_R2=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_unpaired.fastq.gz", sample=SAMPLES)
    output:
        dir=directory(scratch_dir + "01-analysis/10-assembled-megahit/")
    run:
        R1_paired_list = ",".join(map(str, input.R1))
        R2_paired_list = ",".join(map(str, input.R2))
        R1_unpaired_list = ",".join(map(str, input.unpaired_R1))
        R2_unpaired_list = ",".join(map(str, input.unpaired_R2))
        shell("module load Miniconda3")
        shell("source activate /home/ad14556/conda-envs/envs/megahit")
        shell("megahit -1 {R1_paired_list} -2 {R2_paired_list} -r {R1_unpaired_list},{R2_unpaired_list} -o {output.dir}")