rule megahit:
    input:
        multiqc=scratch_dir + "01-analysis/05-trimmed-multiqc/multiqc.html",
        R1=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz", sample=SAMPLES),
        R2=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz", sample=SAMPLES),
        unpaired_R1=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_unpaired.fastq.gz", sample=SAMPLES),
        unpaired_R2=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_unpaired.fastq.gz", sample=SAMPLES)
    output:
        dir=directory(scratch_dir + "01-analysis/10-assembled-megahit/")
    params:
        R1_paired_list = lambda wildcards, input: ','.join(input.R1),
        R2_paired_list = lambda wildcards, input: ','.join(input.R2),
#        R1_unpaired_list = lambda wildcards, input: ','.join(input.unpaired_R1),
#        R2_unpaired_list = lambda wildcards, input: ','.join(input.unpaired_R2)
    conda:
        "../envs/megahit.yaml"
    shell:
        '''
        megahit -1 {params.R1_paired_list} -2 {params.R2_paired_list} -o {output.dir} -t 24 --presets meta-sensitive
        '''

