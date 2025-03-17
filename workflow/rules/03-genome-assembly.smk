rule megahit:
    input:
        multiqc=scratch_dir + "01-analysis/05-trimmed-multiqc/multiqc.html",
        R1=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz", 
        R2=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz"
    output:
        dir=directory(scratch_dir + "01-analysis/08-assembled-megahit/{sample}")
    conda:
        "../envs/megahit.yaml"
    shell:
        '''
        megahit -1 {input.R1} -2 {input.R2} -o {output.dir} -t 24 --presets meta-sensitive
        '''

rule busco:
    input:
        scratch_dir + "01-analysis/08-assembled-megahit/{sample}"
    output:
        path=directory(scratch_dir + "01-analysis/11-busco/{sample}")
    log:
        scratch_dir + "03-log/11-busco/{sample}_busco_euk.log"
    params:
        name="busco_host",
        input="final.contigs.fa",
        mode="genome",
        extra="--auto-lineage-euk nematoda_od10"
    threads: 8
    shell:
        '''
        module load BUSCO
        busco -f -i {input}/{params.input} --mode {params.mode} -l /home/ad14556/snakemake-pipelines/metagenomicsPipeline/workflow/databases/nematoda_odb10 -o {params.name} --out_path {output.path}
        '''

rule mitoz:
    input:
        dir=scratch_dir + "01-analysis/08-assembled-megahit/{sample}"
    output:
        dir=directory(scratch_dir + "01-analysis/09-mitoz/{sample}")
    params:
        contigfile="final.contigs.fa",
        prefix="{sample}"
    shell:
        '''
        module load MitoZ
        mkdir -p {output.dir}
        cd {output.dir}
        mitoz findmitoscaf --fastafile {input.dir}/{params.contigfile} --workdir {output.dir} --outprefix {output.dir}/{params.prefix} \
            --clade Nematoda --requiring_taxa Nematoda --min_abundance 0 --slow_search
        '''

rule euk_barrnap:
    input:
        dir=scratch_dir + "01-analysis/08-assembled-megahit/{sample}"
    output:
        fasta=scratch_dir + "01-analysis/10-barrnap/{sample}.fasta",
        gff=scratch_dir + "01-analysis/10-barrnap/{sample}.gff"
    params:
        kingdom="euk", # options are bac, arc, euk, mito
        fasta="final.contigs.fa"
#    conda:
#        "../envs/barrnap.yaml"
    shell:
        '''
        module load barrnap
        barrnap -k euk -o {output.fasta} < {input.dir}/{params.fasta} > {output.gff}
        '''
