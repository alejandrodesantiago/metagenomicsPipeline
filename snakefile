# Author: Alejandro De Santiago
# Affiliation: University of Georgia
# How to Run: snakemake --cores 1 --use-conda -s snakemakeFileName.smk

configfile: "config/config.yaml" # edit config file with your project information

input_dir       = config['input_path']
output_dir      = config['output_path']
scratch_dir     = config['scratch_path']
project_name    = config['project_name']
adapters        = config['adapters']

FILES,= glob_wildcards(input_dir + "{file}.fastq.gz")
SAMPLES, = glob_wildcards(input_dir + "{sample}.R1.fastq.gz")

rule all:
    input:
        initial_multiqc=scratch_dir + "01-analysis/02-initial-multiqc/multiqc.html", # run initial multiqc
        trimmed_multiqc=scratch_dir + "01-analysis/05-trimmed-multiqc/multiqc.html", # run multiqc on trimmed dataset
        metaspades=expand(scratch_dir + "01-analysis/06-assembled-metaspades/{sample}", sample=SAMPLES), # run genome assemblies
        metaphlan=expand(scratch + "01-analysis/08-metaphlan/{sample}_taxonomy_profile.txt", sample=SAMPLES) # taxonomy profiling of short reads

##### load rules #####
include: "workflow/rules/quality_control.smk"		# step 1-5
include: "workflow/rules/taxonomic_profiling.smk"	# step 6-7
include: "workflow/rules/genome_assembly.smk"		# step 8
include: "workflow/rules/assembly_qc.smk"		# step 9
include: "workflow/rules/binning.smk"			# step 10-17


# metagenomic assembly
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
    params: ""
    conda:
        "workflow/envs/spades.yaml"
    shell:
        '''
        cat {input.unpaired_R1} {input.unpaired_R2} > {output.merged_unpaired}
        spades.py --meta -1 {input.R1} -2 {input.R2} -s {output.merged_unpaired} --threads 16 -m 500 -o {output.dir}
        '''
