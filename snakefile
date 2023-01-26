# Author: Alejandro De Santiago
# Affiliation: University of Georgia
# How to Run: snakemake --cores 1 --use-conda -s snakemakeFileName.smk

configfile: "config/config.yaml" # edit config file with your project information

input_dir       = config['input_path']
output_dir      = config['output_path']
scratch_dir     = config['scratch_path']
project_name    = config['project_name']
adapters        = config['adapters']

(FILES,) = glob_wildcards(input_dir + "{file}.fastq.gz")
(SAMPLES,) = glob_wildcards(input_dir + "{sample}.R1.fastq.gz")

rule all:
    input:
        initial_multiqc=scratch_dir + "01-analysis/02-initial-multiqc/multiqc.html", # run initial multiqc
        trimmed_multiqc=scratch_dir + "01-analysis/05-trimmed-multiqc/multiqc.html", # run multiqc on trimmed dataset
#        metaspades=expand(scratch_dir + "01-analysis/06-assembled-metaspades/{sample}", sample=SAMPLES), # run genome assemblies
#        flash=expand(scratch_dir + "01-analysis/06-flash/{sample}_merged.fastq.gz", sample=SAMPLES),
        metaphlan=expand(scratch_dir + "01-analysis/07-metaphlan/{sample}_taxonomy_profile.txt", sample=SAMPLES) # taxonomy profiling of short reads

##### load rules #####
include: "workflow/rules/quality_control.smk"		# step 1-5
include: "workflow/rules/taxonomic_profiling.smk"	# step 6-7
#include: "workflow/rules/genome_assembly.smk"		# step 8
#include: "workflow/rules/assembly_qc.smk"		# step 9
#include: "workflow/rules/binning.smk"			# step 10-17
