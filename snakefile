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
(SAMPLES,) = glob_wildcards(input_dir + "{sample}_R1.fastq.gz")

rule all:
    input:
        ## QC WITH MULTIQC AND TRIMMOMATIC ##
        multiqc_rawdata=scratch_dir + "01-analysis/02-initial-multiqc/multiqc.html",
        multiqc_trimmed=scratch_dir + "01-analysis/05-trimmed-multiqc/multiqc.html",
        ## METAPHLAN4 FOR RAPID TAXONOMIC ID
        metaphlan=expand(scratch_dir + "01-analysis/07-metaphlan/{sample}_taxonomy_profile.txt", sample=SAMPLES),
        ## GENOME ASSEMBLY, MITOCHONDRIA, RIBOSOMAL GENES, AND QUALITY MULTIQC AND METAQUAST ##
        mitoz=expand(scratch_dir + "01-analysis/09-mitoz/{sample}", sample=SAMPLES),
	barrnap=expand(scratch_dir + "01-analysis/10-barrnap/{sample}.fasta", sample=SAMPLES),
        busco=expand(scratch_dir + "01-analysis/11-busco/{sample}", sample=SAMPLES),
        ## BINNING EUKARYOTES WITH DASTOOL, METABAT, CONCOCT, AND MAXBIN2 ##
        gtdbTK=expand(scratch_dir + "01-analysis/14-gtdb/{sample}", sample=SAMPLES),
        checkm=expand(scratch_dir + "01-analysis/15-checkm/{sample}", sample=SAMPLES)
#         depth=scratch_dir + "01-analysis/14-eukmags/02-metabat2/depth.txt"
#        euk_metabat=expand(scratch_dir + "01-analysis/14-eukmags/04-dastool/{sample}.metabat.scaffolds2bin.tsv", sample=SAMPLES),
#        euk_concoct=expand(scratch_dir + "01-analysis/14-eukmags/04-dastool/{sample}.concoct.scaffolds2bin.tsv", sample=SAMPLES),
#        euk_dastool=expand(scratch_dir + "01-analysis/14-eukmags/04-dastool/{sample}",sample=SAMPLES),
#        ## BINNING PROKARYOTES WITH DASTOOL, METABAT, CONCOCT, AND MAXBIN2 ##
#        pro_metabat=expand(scratch_dir + "01-analysis/15-bacmags/04-dastool/{sample}.metabat.scaffolds2bin.tsv", sample=SAMPLES),
#        pro_concoct=expand(scratch_dir + "01-analysis/15-bacmags/04-dastool/{sample}.concoct.scaffolds2bin.tsv", sample=SAMPLES),
#        pro_dastool=expand(scratch_dir + "01-analysis/15-bacmags/04-dastool/{sample}",sample=SAMPLES),

##### load rules #####
include: "workflow/rules/01-quality-control.smk"		# step 1 - Quality Control Using Trimmomatic and MultiQC
include: "workflow/rules/02-taxonomic-profiling.smk"		# step 2 - Taxonomic profiling using Kraken and Metaphlan
include: "workflow/rules/03-genome-assembly.smk"		# step 3 - Assembly using MEGAHIT
include: "workflow/rules/04-binning-pro.smk"			# step 4 - Bin Prokaryote Contigs
include: "workflow/rules/05-binning-quality-pro.smk"            # step 5 - Bacterial MAG Quality
#include: "workflow/rules/06-annotate-pro.smk"			# step 9 - Annotate Prokaryote Bins Using DRAM
