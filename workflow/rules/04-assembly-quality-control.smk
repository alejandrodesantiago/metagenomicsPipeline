rule metaquast:
    input:
        megahit=scratch_dir + "01-analysis/10-assembled-megahit/final.contigs.fa",
        R1=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz",
        R2=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz",
        unpaired_R1=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_unpaired.fastq.gz",
        unpaired_R2=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_unpaired.fastq.gz"
    output:
        dir=directory(scratch_dir + "01-analysis/11-metaquast/{sample}_assembly_quality")
    params: "final.contigs.fa"
    conda:
        "../envs/quast.yaml"
    shell:
        '''
        metaquast {input.megahit}/{params} -1 {input.R1} -2 {input.R2} -o {output} --max-ref-number 0
        '''

# multiqc for quast
rule assembly_multiqc:
    input:
        expand(scratch_dir + "01-analysis/11-metaquast/{sample}_assembly_quality", sample=SAMPLES)
    output:
        scratch_dir + "01-analysis/12-assembly-multiqc/multiqc.html"
    params: ""
    log:
        scratch_dir + "03-log/12-assembly-multiqc.log"
    wrapper:
        "v1.5.0/bio/multiqc"
