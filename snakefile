# Author: Alejandro De Santiago
# Affiliation: University of Georgia
# How to Run: snakemake --cores 1 --use-conda -s snakemakeFileName.smk

configfile: "config/config.yaml" # edit config file with your project information
clusterfile: "config/cluster.yaml" # edit cluster file

input_dir       = config['input_path']
output_dir      = config['output_path']
scratch_dir     = config['scratch_path']
project_name    = config['project_name']
adapters        = config['adapters']

(FILES,) = glob_wildcards(input_dir + "{file}.fastq.gz")
(SAMPLES,) = glob_wildcards(input_dir + "{sample}.R1.fastq.gz")

rule all:
    input:
        # quality control with multiqc and trimmomatic
        multiqc_rawdata=scratch_dir + "01-analysis/02-initial-multiqc/multiqc.html",
        multiqc_trimmed=scratch_dir + "01-analysis/05-trimmed-multiqc/multiqc.html",
        # assembly quality
        metaquast=scratch_dir + "01-analysis/11-assembly-multiqc/multiqc.html",
        # binning eukaryotes
        metabat=expand(scratch_dir + "01-analysis/14-eukmags/04-dastool/{sample}.metabat.scaffolds2bin.tsv", sample=SAMPLES),
        concoct=expand(scratch_dir + "01-analysis/14-eukmags/04-dastool/{sample}.concoct.scaffolds2bin.tsv", sample=SAMPLES),
        dastool=expand(scratch_dir + "01-analysis/14-eukmags/04-dastool/{sample}",sample=SAMPLES)

##### load rules #####
include: "workflow/rules/01-quality-control.smk"			# step 1 - Quality Control Using Trimmomatic and MultiQC
include: "workflow/rules/02-taxonomic-profiling.smk"		# step 2 - Taxonomic profiling using Kraken and Metaphlan
include: "workflow/rules/03-genome-assembly.smk"			# step 3 - Assembly using metaSPAdes and MEGAHIT
include: "workflow/rules/04-assembly-quality-control.smk"	# step 4 - Assembly Quality using MetaQuast and MultiQC
include: "workflow/rules/05-binning-eukrep.smk"				# step 5 - Bin eukaryote and prokaryote contigs with Eukrep
include: "workflow/rules/06-binning-euk.smk"				# step 6 - Bin eukaryote reads
#include: "workflow/rules/07-binning-pro.smk"
#include: "workflow/rules/08-annotate-euk.smk"
#include: "workflow/rules/09-annotate-pro.smk"
#include: "workflow/rules/10-binning-quality-euk.smk"
#include: "workflow/rules/11-binning-quality-pro.smk"
#include: "workflow/rules/12-extract-18S.smk" 
#include: "workflow/rules/13-visualization-anvio.smk"
