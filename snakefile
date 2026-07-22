# Author: Alejandro De Santiago
# Affiliation: University of Georgia
# How to Run: snakemake --cores 1 --use-conda -s snakemakeFileName.smk

configfile: "config/config.yaml" # edit config file with your project information
clusterfile: "config/cluster.yaml" # edit cluster file

# paths
input_dir       = config['input_path']
output_dir      = config['output_path']
scratch_dir     = config['scratch_path']

# databases
adapters        = config['adapters']
checkm_db       = config['checkm_db']
gtdbtk_db       = config['gtdbtk_db']
kraken_db       = config['kraken_db']
busco_db        = config['busco_db']
busco_db_name   = config['busco_db_name']
completebin_db  = config['completebin_db']
magpurify_db    = config['magpurify_db']

# bin qual
completeness    = config['completeness']
contamination   = config['contamination']

# sample and file names
(FILES,) = glob_wildcards(input_dir + "{file}.fastq.gz")
(SAMPLES,) = glob_wildcards(input_dir + "{sample}_R1.fastq.gz")

rule all:
    input:
        ## QC WITH MULTIQC AND TRIMMOMATIC
        reads_quality_trim=output_dir + "01-reads-quality-trimmed-multiqc-trimmed.html",
        ## Assemblying Metagenomes with MEGAHIT 
        metagenome_assembly=expand(output_dir + "02-metagenome-assembly/{sample}-metagenome.fa", sample=SAMPLES),
        ## HOST SKIMMING WITH BUSCO, BARRNAP, AND MITOZ
        host_genome_skim=expand(output_dir + "03-host-genome-skim/{sample}/{sample}-busco-summary.txt", sample=SAMPLES),
        ## BINNING WITH COMEBIN,METABAT2,AND DASTOOL
        microbiome_bins=expand(scratch_dir + "01-analysis/13-bacmags/{sample}-microbiome-bins-done.txt", sample=SAMPLES),
        microbiome_bins_quality=output_dir + "05-microbiome-bins-quality.tsv",
        ## QUANTIFICATION OF MICROBIAL BINS
        microbiome_bins_quant=output_dir + "06-microbiome-bins-quant-rpkm-no-log.csv"

##### load rules #####
include: "workflow/rules/01-quality-control.smk"		# step 1 - Quality Control Using Trimmomatic and MultiQC
include: "workflow/rules/02-genome-assembly.smk"		# step 2 - Assembly using MEGAHIT
include: "workflow/rules/03-host-genome.smk"			# step 3 - Host Genome Skimming
include: "workflow/rules/04-binning-pro.smk"			# step 4 - Binning Microbiome
include: "workflow/rules/05-binning-quality-pro.smk"		# step 5 - Filtering High Quality Bins
include: "workflow/rules/06-binning-read-map.smk"		# step 6 - Quantify Microbiome
