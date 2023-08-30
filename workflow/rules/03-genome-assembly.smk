rule megahit:
    input:
        multiqc=scratch_dir + "01-analysis/05-trimmed-multiqc/multiqc.html",
        R1=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz", sample=SAMPLES),
        R2=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz", sample=SAMPLES),
        unpaired_R1=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_unpaired.fastq.gz", sample=SAMPLES),
        unpaired_R2=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_unpaired.fastq.gz", sample=SAMPLES)
    output:
        dir=scratch_dir + "01-analysis/10-assembled-megahit/"
#        contigs=scratch_dir + "01-analysis/10-assembled-megahit/final.contigs.fa"
#    conda:
#        "../envs/megahit.yaml"
    run:
        R1_paired_list = ",".join(map(str, input.R1))
        R2_paired_list = ",".join(map(str, input.R2))
        R1_unpaired_list = ",".join(map(str, input.unpaired_R1))
        R2_unpaired_list = ",".join(map(str, input.unpaired_R2))
        shell("module load MEGAHIT")
        shell("megahit -1 {R1_paired_list} -2 {R2_paired_list} -r {R1_unpaired_list},{R2_unpaired_list} -o {output.dir} --presets meta-sensitive")


#rule megahit:
#    input:
#        multiqc=scratch_dir + "01-analysis/05-trimmed-multiqc/multiqc.html",
#        R1=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz",
#        R2=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz",
#        unpaired_R1=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_unpaired.fastq.gz",
#        unpaired_R2=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_unpaired.fastq.gz"
#    output:
#        dir=directory(scratch_dir + "01-analysis/10-assembled-megahit/{sample}")
#    conda:
#        "../envs/megahit.yaml"
#    shell:
#        '''
#        megahit -1 {input.R1} -2 {input.R2} -r {input.unpaired_R1},{input.unpaired_R1} -t 24 -o {output.dir}
#        '''

#rule metaspades:
#    input:
#        multiqc=scratch_dir + "01-analysis/05-trimmed-multiqc/multiqc.html",
#        R1=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz",
#        R2=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz",
#        unpaired_R1=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_unpaired.fastq.gz",
#        unpaired_R2=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_unpaired.fastq.gz"
#    output:
#        merged_unpaired=scratch_dir + "01-analysis/03-trimmomatic/{sample}_merged_unpaired.fastq.gz",
#        dir=directory(scratch_dir + "01-analysis/09-assembled-metaspades/{sample}")
#    params:
#    conda:
#        "../envs/spades.yaml"
#    shell:
#        '''
#        cat {input.unpaired_R1} {input.unpaired_R2} > {output.merged_unpaired}
#        spades.py --meta -1 {input.R1} -2 {input.R2} -s {output.merged_unpaired} --threads 16 -m 1000 -o {output.dir}
#        '''
