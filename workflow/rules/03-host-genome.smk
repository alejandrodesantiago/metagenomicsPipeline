rule busco:
    input:
        contigs=scratch_dir + "01-analysis/07-assembled-metaspades/{sample}/scaffolds_1000.fa"        
    output:
        txt=scratch_dir + "01-analysis/07-busco/{sample}/busco_host/short_summary.specific." + busco_db_name + ".busco_host.txt"
    log:
        err=scratch_dir + "03-log/07-busco/{sample}-busco-euk.err",
        out=scratch_dir + "03-log/07-busco/{sample}-busco-euk.out"
    params:
        dir=scratch_dir + "01-analysis/07-busco/{sample}",
        name="busco_host",
        mode="genome",
        extra=busco_db 
    threads: config['busco_threads']
    shell:
        '''
        ml BUSCO
        busco -f -i {input.contigs} --cpu {threads} --metaeuk --mode {params.mode} -l {params.extra} -o {params.name} --out_path {params.dir} > {log.out} 2> {log.err}
        '''

rule mitoz:
    input:
        contigs=scratch_dir + "01-analysis/07-assembled-metaspades/{sample}/scaffolds_1000.fa"
    output:
        cir=scratch_dir + "01-analysis/08-mitoz/{sample}/{sample}.start2end_for-circular-mt-only",
        hit=scratch_dir + "01-analysis/08-mitoz/{sample}/{sample}.result/{sample}.hmmtblout.besthit.sim.filtered-by-taxa.fa",
        mit=scratch_dir + "01-analysis/08-mitoz/{sample}/{sample}.result/{sample}.mitogenome.fa",
        lap=scratch_dir + "01-analysis/08-mitoz/{sample}/{sample}.result/{sample}.overlap_information"
    log:
        err=scratch_dir + "03-log/08-mitoz/{sample}-mitoz-find.err",
        out=scratch_dir + "03-log/08-mitoz/{sample}-mitoz-find.out"
    params:
        dir=scratch_dir + "01-analysis/08-mitoz/{sample}",
        prefix="{sample}",
        taxa=config['host_taxa'],
        abundance=config['mitoz_abundance']
    threads: config['mitoz_threads']
    conda: 
        "../envs/mitoz.yaml"
    shell:
        '''
        mkdir -p {params.dir}
        cd {params.dir}
        mitoz findmitoscaf --fastafile {input.contigs} --workdir {params.dir} --outprefix {params.dir}/{params.prefix} \
            --clade {params.taxa} --requiring_taxa {params.taxa} --min_abundance {params.abundance} --slow_search --thread_number {threads} > {log.out} 2> {log.err}
        touch {output.mit}
        touch {output.cir}
        touch {output.hit}
        touch {output.lap}
        '''

rule euk_barrnap:
    input:
        contigs=scratch_dir + "01-analysis/07-assembled-metaspades/{sample}/scaffolds_1000.fa"
    output:
        fasta=scratch_dir + "01-analysis/09-barrnap/{sample}.fasta",
        gff=scratch_dir + "01-analysis/09-barrnap/{sample}.gff"
    params:
        kingdom="euk", # options are bac, arc, euk, mito
    conda:
        "../envs/barrnap.yaml"
    shell:
        '''
        barrnap -k euk -o {output.fasta} < {input.contigs} > {output.gff}
        '''

rule all_genome_skims_processed:
    input:
        barrnap=expand(scratch_dir + "01-analysis/09-barrnap/{sample}.fasta", sample=SAMPLES),
        mitoz=expand(scratch_dir + "01-analysis/08-mitoz/{sample}/{sample}.result/{sample}.mitogenome.fa", sample=SAMPLES),
        busco=expand(scratch_dir + "01-analysis/07-busco/{sample}/busco_host/short_summary.specific." + busco_db_name + ".busco_host.txt", sample=SAMPLES)
    output:
        done=scratch_dir + "01-analysis/10-all-samples-finished-processing.txt"
    shell:
        "echo 'All genome skims have been processed.' > {output.done}"

rule host_genome_skim:
    input:
        processed=scratch_dir + "01-analysis/10-all-samples-finished-processing.txt",
        busco=scratch_dir + "01-analysis/07-busco/{sample}/busco_host/short_summary.specific." + busco_db_name + ".busco_host.txt",
        mitoz_cir=scratch_dir + "01-analysis/08-mitoz/{sample}/{sample}.start2end_for-circular-mt-only",
        mitoz_hit=scratch_dir + "01-analysis/08-mitoz/{sample}/{sample}.result/{sample}.hmmtblout.besthit.sim.filtered-by-taxa.fa",
        mitoz_mit=scratch_dir + "01-analysis/08-mitoz/{sample}/{sample}.result/{sample}.mitogenome.fa",
        mitoz_lap=scratch_dir + "01-analysis/08-mitoz/{sample}/{sample}.result/{sample}.overlap_information",
        bar_fasta=scratch_dir + "01-analysis/09-barrnap/{sample}.fasta",
        bar_gff=scratch_dir + "01-analysis/09-barrnap/{sample}.gff"
    output:
        complete=directory(output_dir + "03-host-genome-skim/{sample}/{sample}-busco/single_copy_busco_sequences/"),
        fragment=directory(output_dir + "03-host-genome-skim/{sample}/{sample}-busco/fragmented_busco_sequences/"),
        multicopy=directory(output_dir + "03-host-genome-skim/{sample}/{sample}-busco/multi_copy_busco_sequences/"),
        busco=output_dir + "03-host-genome-skim/{sample}/{sample}-busco-summary.txt",
        mitoz_cir=output_dir + "03-host-genome-skim/{sample}/{sample}-mitoz-circular-mito.fa",
        mitoz_hit=output_dir + "03-host-genome-skim/{sample}/{sample}-mitoz-besthits-mito.txt",
        mitoz_mit=output_dir + "03-host-genome-skim/{sample}/{sample}-mitoz-mitogenome.fa",
        mitoz_lap=output_dir + "03-host-genome-skim/{sample}/{sample}-mitoz-overlap-info.txt",
        bar_fasta=output_dir + "03-host-genome-skim/{sample}/{sample}-barrnap-rRNA.fasta",
        bar_gff=output_dir + "03-host-genome-skim/{sample}/{sample}-barrnap-rRNA.gff"
    params:
        complete=scratch_dir + "01-analysis/07-busco/{sample}/busco_host/run_" + busco_db_name + "/busco_sequences/single_copy_busco_sequences/",
        fragment=scratch_dir + "01-analysis/07-busco/{sample}/busco_host/run_" + busco_db_name + "/busco_sequences/fragmented_busco_sequences/",
        multicopy=scratch_dir + "01-analysis/07-busco/{sample}/busco_host/run_" + busco_db_name + "/busco_sequences/multi_copy_busco_sequences/",
    shell:
        '''
        cp -r {params.complete} {output.complete}
        cp -r {params.fragment} {output.fragment}
        cp -r {params.multicopy} {output.multicopy}
        cp {input.busco} {output.busco}
        cp {input.mitoz_cir} {output.mitoz_cir}
        cp {input.mitoz_hit} {output.mitoz_hit}
        cp {input.mitoz_mit} {output.mitoz_mit}
        cp {input.mitoz_lap} {output.mitoz_lap}
        cp {input.bar_fasta} {output.bar_fasta}
        cp {input.bar_gff} {output.bar_gff}
        '''

