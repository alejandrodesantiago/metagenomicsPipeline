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
        trimmomatic=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_r1_paired.fastq.gz", sample=SAMPLES) # run trimmomatic

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
        multiqc=expand(scratch_dir + "01-analysis/02-initial-multiqc/multiqc.html"),
        r1=input_dir + "{sample}.R1.fastq.gz",
        r2=input_dir + "{sample}.R2.fastq.gz"
    output:
        r1_paired=scratch_dir + "01-analysis/03-trimmomatic/{sample}_r1_paired.fastq.gz",
        r2_paired=scratch_dir + "01-analysis/03-trimmomatic/{sample}_r2_paired.fastq.gz",
        r1_unpaired=scratch_dir + "01-analysis/03-trimmomatic/{sample}_r1_unpaired.fastq.gz",
        r2_unpaired=scratch_dir + "01-analysis/03-trimmomatic/{sample}_r2_unpaired.fastq.gz"
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