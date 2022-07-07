# Author: Alejandro De Santiago
# Affiliation: University of Georgia
# How to Run: snakemake --cores 1 --use-conda -s snakemakeFileName.smk

configfile: "config.yaml" # edit config file with your project information

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
        assembly_multiqc=scratch_dir + "01-analysis/07-assembly-multiqc/multiqc.html", # need to run multiqc for assembly quality
        metaspades=expand(scratch_dir + "01-analysis/05-assembled-metaspades/{sample}", sample=SAMPLES),
        dastool_euk=expand(scratch_dir + "01-analysis/09-binned-euk/{sample}/dastool/{sample}_dastool", sample=SAMPLES)

#        metaquast=expand(scratch_dir + "01-analysis/06-metaquast/{sample}_assembly_quality", sample=SAMPLES), # need to run metaquast for assembly quality
#        eukrep=expand(scratch_dir + "01-analysis/08-EukRep/{sample}/{sample}_euk.fasta", sample=SAMPLES)


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
#        file=scratch_dir + "01-analysis/06-assembled-metaspades/{sample}/contigs.fasta",
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
# may need to provide a reference
rule metaquast:
    input:
        R1=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz", sample=SAMPLES),
        R2=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz", sample=SAMPLES)
    output:
        dir=directory(scratch_dir + "01-analysis/06-metaquast/{sample}_assembly_quality")
    params:
        input=expand(scratch_dir + "01-analysis/05-assembled-metaspades/{sample}/contigs.fasta", sample=SAMPLES)
    conda:
        "envs/quast.yaml"
    shell:
        '''
        metaquast {params.input} -1 {input.R1} -2 {input.R2} -o {output}
        '''

# multiqc for quast
rule assembly_multiqc:
    input:
        expand(scratch_dir + "01-analysis/06-metaquast/{sample}_assembly_quality", sample=SAMPLES)
    output:
        scratch_dir + "01-analysis/07-assembly-multiqc/multiqc.html"
    params: ""
    log:
        scratch_dir + "03-log/07-assembly-multiqc.log"
    wrapper:
        "v1.5.0/bio/multiqc"

# binning eukaryotes
rule eukrep:
 #   input:
 #       contig=expand(scratch_dir + "01-analysis/06-assembled-metaspades/{sample}/contigs.fasta", sample=SAMPLES)
    output:
       euk=scratch_dir + "01-analysis/08-EukRep/{sample}/{sample}_euk.fasta",
       pro=scratch_dir + "01-analysis/08-EukRep/{sample}/{sample}_pro.fasta"
    params:
        min_contig = 1000, # due to fragmented genomes
        file=expand(scratch_dir + "01-analysis/05-assembled-metaspades/{sample}/contigs.fasta", sample=SAMPLES)
    conda:
        "envs/eukrep.yaml"
    shell:
        '''
        EukRep -i {params.file} -o {output.euk} --prokarya {output.pro} --min {params.min_contig}
        '''


# binning eukaryotes
rule mapReads:
    input:
        euk=expand(scratch_dir + "01-analysis/08-EukRep/{sample}/{sample}_euk.fasta", sample=SAMPLES),
        r1=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz", sample=SAMPLES),
        r2=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz", sample=SAMPLES),
        r1_unpaired=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_unpaired.fastq.gz", sample=SAMPLES),
        r2_unpaired=expand(scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_unpaired.fastq.gz", sample=SAMPLES)
    output:
        euk = scratch_dir + "01-analysis/09-binned-euk/{sample}/mappedReads/alignment_euk_merged_final.bam",
        euk_paired = scratch_dir + "01-analysis/09-binned-euk/{sample}/mappedReads/alignment_euk_paired.bam",
        euk_unpaired = scratch_dir + "01-analysis/09-binned-euk/{sample}/mappedReads/alignment_euk_unpaired.bam",
        single_reads = scratch_dir + "01-analysis/03-trimmomatic/{sample}_single_reads.fastq.gz"
    params:
        threads=1
#   log: ""
    conda:
        "envs/bwa-samtools.yaml"
    shell:
        '''
        cat {input.r1_unpaired} {input.r2_unpaired} > {output.single_reads}
        bwa index {input.euk} 
        bwa mem -t {params.threads} {input.euk} {input.r1} {input.r2} | samtools sort -o {output.euk_paired}
        bwa mem -t {params.threads} {input.euk} {output.single_reads} | samtools sort -o {output.euk_unpaired}
        samtools merge {output.euk} {output.euk_paired} {output.euk_unpaired}
        samtools index {output.euk}
        '''

rule euk_metabat:
    input:
        fasta=scratch_dir + "01-analysis/08-EukRep/{sample}/{sample}_euk.fasta",
        bam=scratch_dir + "01-analysis/09-binned-euk/{sample}/mappedReads/alignment_euk_merged_final.bam"
    output:
        depth=scratch_dir + "01-analysis/09-binned-euk/{sample}/metabat2/depth.txt"
    params:
        bin=scratch_dir + "01-analysis/09-binned-euk/{sample}/metabat2/bin/{sample}_bin"
#    log: ""
    conda:
        "envs/metabat2.yaml"
    shell:
        '''
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}
        metabat2 -i {input.fasta} -a {output.depth} -o {params.bin}
        '''

rule concoct:
    input:
        fasta = scratch_dir + "01-analysis/08-EukRep/{sample}/{sample}_euk.fasta",
        bam = scratch_dir + "01-analysis/09-binned-euk/{sample}/mappedReads/alignment_euk_merged_final.bam"
    output:
        bed = scratch_dir + "01-analysis/09-binned-euk/{sample}/concoct/concoct.bed",
        fasta = scratch_dir + "01-analysis/09-binned-euk/{sample}/concoct/concoct.fa",
        depth = scratch_dir + "01-analysis/09-binned-euk/{sample}/concoct/depth.tsv",
        dir = directory(scratch_dir + "01-analysis/09-binned-euk/{sample}/concoct/{sample}"),
        bin = directory(scratch_dir + "01-analysis/09-binned-euk/{sample}/concoct/bin/{sample}")
#    params: ""
#    log: ""
    conda:
        "envs/concoct.yaml"
    shell:
        '''
        cut_up_fasta.py {input.fasta} -c 10000 -o 0 --merge_last -b {output.bed} > {output.fasta}
        concoct_coverage_table.py {output.bed} {input.bam} > {output.depth}
        concoct --composition_file {output.fasta} --coverage_file {output.depth} -b {output.dir}
        merge_cutup_clustering.py {output.dir}_clustering_gt1000.csv > {output.dir}_clustering_merged.csv
        extract_fasta_bins.py {input.fasta} {output.dir}_clustering_merged.csv --output_path {output.bin}
        '''

rule euk_dastool:
    input:
        contigs = expand(scratch_dir + "01-analysis/08-EukRep/{sample}/{sample}_euk.fasta", sample=SAMPLES),
        concoct = expand(scratch_dir + "01-analysis/09-binned-euk/{sample}/concoct/bin", sample = SAMPLES),
        metabat = expand(scratch_dir + "01-analysis/09-binned-euk/{sample}/metabat2/depth.txt", sample = SAMPLES)
    output:
        metabat=scratch_dir + "01-analysis/09-binned-euk/{sample}/dastool/metabat.scaffolds2bin.tsv",
        concoct=scratch_dir + "01-analysis/09-binned-euk/{sample}/dastool/concoct.scaffolds2bin.tsv",
        dastool=scratch_dir + "01-analysis/09-binned-euk/{sample}/dastool/{sample}_dastool"
    params:
        metabat = expand(scratch_dir + "01-analysis/09-binned-euk/{sample}/metabat2/bin",sample=SAMPLES)
#    log: ""
    conda:
        "envs/dastool.yaml"
    shell:
        '''
        Fasta_to_Scaffolds2Bin.sh -i {params.metabat} -e fa > {output.metabat}
        Fasta_to_Scaffolds2Bin.sh -i {input.concoct} -e fa > {output.concoct}
        DAS_Tool -i {output.metabat},{output.concoct} -l metabat,concoct -c {input.contigs} -o {params} --write_bins 1
        '''