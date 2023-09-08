#R1_paired_list = ",".join(map(str,input.R1))
#R2_paired_list = ",".join(map(str,input.R2))
#R1_unpaired_list = ",".join(map(str,input.unpaired_R1))
#R2_unpaired_list = ",".join(map(str,input.unpaired_R2))

rule megahit:
    input:
        multiqc=scratch_dir + "01-analysis/05-trimmed-multiqc/multiqc.html",
        R1=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz", sample=SAMPLES),
        R2=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz", sample=SAMPLES),
        unpaired_R1=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_unpaired.fastq.gz", sample=SAMPLES),
        unpaired_R2=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_unpaired.fastq.gz", sample=SAMPLES)
    output:
        dir=directory(scratch_dir + "01-analysis/10-assembled-megahit/")
    conda:
        "../envs/megahit.yaml"
    shell:
        '''
        megahit -1 {R1} -2 {R2} -r {unpaired_R1},{unpaired_R2} -o {output.dir}
        '''