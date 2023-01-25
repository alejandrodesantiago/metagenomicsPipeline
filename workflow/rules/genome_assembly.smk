rule metaspades:
    input:
        multiqc=scratch_dir + "01-analysis/05-trimmed-multiqc/multiqc.html",
        R1=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz",
        R2=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz",
        unpaired_R1=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_unpaired.fastq.gz",
        unpaired_R2=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_unpaired.fastq.gz"
    output:
        merged_unpaired=scratch_dir + "01-analysis/03-trimmomatic/{sample}_merged_unpaired.fastq.gz",
        dir=directory(scratch_dir + "01-analysis/06-assembled-metaspades/{sample}")
    params:
    conda:
        "../envs/spades.yaml"
    shell:
        '''
        cat {input.unpaired_R1} {input.unpaired_R2} > {output.merged_unpaired}
        spades.py --meta -1 {input.R1} -2 {input.R2} -s {output.merged_unpaired} --threads 16 -m 500 -o {$
        '''
