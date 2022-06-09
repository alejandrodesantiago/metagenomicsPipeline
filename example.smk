# Author: Alejandro De Santiago
# Affiliation: University of Georgia
# How to Run: snakemake --cores 1 --use-conda -s example.smk
configfile: "config.yaml"

input_dir       = config['input_path']
output_dir      = config['output_path']
scratch_dir     = config['scratch_path']
project_name    = config['project_name']
adapters        = config['adapters']

FILES,= glob_wildcards(input_dir + "{file}.fastq.gz")
SAMPLES, = glob_wildcards(input_dir + "{sample}.R1.fastq.gz")

rule all:
    input:
        initial_multiqc=scratch_dir + "01-analysis/02-initial-multiqc/multiqc.html", # needed to run initial multiqc
        trimmed_multiqc=scratch_dir + "01-analysis/05-trimmed-multiqc/multiqc.html", # needed to run multiqc on trimmed dataset
#        spades=expand(scratch_dir + "01-analysis/06-assembled-metaspades/{sample}", sample=SAMPLES) # run metaspades
        metaquast=expand(scratch_dir + "01-analysis/06-metaquast/{sample}_assembly_quality", sample=SAMPLES)

# quality control visualization (fastqc and multiqc)
rule initial_fastqc:
    input:
        input_dir + "{file}.fastq.gz"
    output:
        html=scratch_dir + "01-analysis/01-intial-fastqc/{file}.html",
        zip=scratch_dir + "01-analysis/01-initial-fastqc/{file}_fastqc.zip"
    params: ""
    log:
        scratch_dir + "03-log/01-initial-fastqc/{file}.log"
    wrapper:
        "v1.5.0/bio/fastqc"

rule initial_multiqc:
    input:
        expand(scratch_dir + "01-analysis/01-initial-fastqc/{file}_fastqc.zip", file=FILES)
    output:
        scratch_dir + "01-analysis/02-initial-multiqc/multiqc.html"
    params: ""
    log:
        scratch_dir + "03-log/02-initial-multiqc.log"
    wrapper:
        "v1.5.0/bio/multiqc"

# trimming data
rule trimmomatic:
    input:
        r1=input_dir + "{sample}.R1.fastq.gz",
        r2=input_dir + "{sample}.R2.fastq.gz"
    output:
        r1=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz",
        r2=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz",
        r1_unpaired=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_unpaired.fastq.gz",
        r2_unpaired=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_unpaired.fastq.gz"
    params:
        trimmer=["ILLUMINACLIP:{}:2:30:7".format(adapters), "LEADING:2", "TRAILING:2", "SLIDINGWINDOW:4:20", "MINLEN:100"],
        extra="",
        compression_level="-9"
    log:
        scratch_dir + "03-log/03-trimmomatic/{sample}.log"
    threads:
        32
    wrapper:
        "v1.5.0/bio/trimmomatic/pe"

# quality control visualization of trimmed data
rule trimmed_fastqc:
    input:
        scratch_dir + "01-analysis/03-trimmomatic/{sample}_{dir}_{pair}.fastq.gz"
    output:
        html=scratch_dir + "01-analysis/04-trimmed-fastqc/{sample}_{dir}_{pair}_fastqc.html",
        zip=scratch_dir + "01-analysis/04-trimmed-fastqc/{sample}_{dir}_{pair}_fastqc.zip"
    params: ""
    log:
        scratch_dir + "03-log/04-trimmed-fastqc/{sample}_{dir}_{pair}.log"
    threads: 2
    wrapper:
        "v1.5.0/bio/fastqc"

rule trimmed_multiqc:
    input:
        expand(scratch_dir + "01-analysis/04-trimmed-fastqc/{sample}_{dir}_{pair}_fastqc.zip", sample=SAMPLES, dir=["R1", "R2"], pair=["paired", "unpaired"])
    output:
        scratch_dir + "01-analysis/05-trimmed-multiqc/multiqc.html"
    params: ""
    log:
        scratch_dir + "03-log/05-initial-multiqc.log"
    wrapper:
        "v1.5.0/bio/multiqc"

# metagenomic assembly
rule metaspades:
    input:
        R1=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz",
        R2=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz"
    output:
        dir=directory(scratch_dir + "01-analysis/05-assembled-metaspades/{sample}")
    params: ""
#    log: ""
    conda:
        "envs/spades.yaml"
    shell:
        '''
        spades.py --meta --pe1-1 {input.R1} --pe1-2 {input.R2} --threads 4 -o {output.dir}
        '''
# assembly quality
rule metaquast:
    input:
        contig=expand(scratch_dir + "01-analysis/06-assembled-metaspades/{sample}", sample=SAMPLES),
        R1=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz", sample=SAMPLES),
        R2=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz", sample=SAMPLES)
    output:
        dir=directory(scratch_dir + "01-analysis/06-metaquast/{sample}_assembly_quality")
    params: ""
    conda:
        "envs/quast.yaml"
    shell:
        '''
        metaquast {input.contig} -1 ${input.R1}.fastq -2 {input.R2} -o {output}
        '''