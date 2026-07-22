rule bwa_index:
    input:
        fasta=scratch_dir + "01-analysis/07-assembled-metaspades/{sample}/scaffolds_1000.fa"
    output:
        done=scratch_dir + "01-analysis/13-bacmags/01-alignment/{sample}-index-finished.txt"
    conda:
        "../envs/bwa-samtools.yaml"
    shell:
        '''
        bwa index {input.fasta}
        touch {output.done}
        '''

rule pro_mapReads:
    input:
        fasta=scratch_dir + "01-analysis/07-assembled-metaspades/{sample}/scaffolds_1000.fa",
        r1=scratch_dir + "01-analysis/04-deduplicated/{sample}_R1_paired.fastq.gz",
        r2=scratch_dir + "01-analysis/04-deduplicated/{sample}_R2_paired.fastq.gz",
        index=scratch_dir + "01-analysis/13-bacmags/01-alignment/{sample}-index-finished.txt"
    output:
        map=scratch_dir + "01-analysis/13-bacmags/01-alignment/{sample}/alignment.bam"
    params:
        dir=scratch_dir + "01-analysis/13-bacmags/01-alignment/{sample}",
        threads=24,
    conda:
        "../envs/bwa-samtools.yaml"
    shell:
        '''
        mkdir -p {params.dir}
        bwa mem -t {params.threads} {input.fasta} {input.r1} {input.r2} | samtools sort -o {output.map} --threads {params.threads}
        samtools index -@ {params.threads} {output.map}
        '''

rule pro_completebin:
    input:
        fasta=scratch_dir + "01-analysis/07-assembled-metaspades/{sample}/scaffolds_1000.fa",
        map=scratch_dir + "01-analysis/13-bacmags/01-alignment/{sample}/alignment.bam",
        db=completebin_db
    output:
        dir=directory(scratch_dir + "01-analysis/13-bacmags/03-completebin/{sample}-final/"),
        temp=directory(scratch_dir + "01-analysis/13-bacmags/03-completebin/{sample}-temp/"),
        done=scratch_dir + "01-analysis/13-bacmags/03-completebin/{sample}-completebin-done.txt"
    log:
        out=scratch_dir + "03-log/completebin/{sample}.out",
        err=scratch_dir + "03-log/completebin/{sample}.err"
    conda:
        "../envs/completebin.yaml"
    threads: 48
    shell:
        '''
        completebin -c {input.fasta} -b {input.map} -o {output.dir} -temp {output.temp} -db {input.db} --min_contig_length 1000 --num_workers {threads} --sec_clu_algo mix || true
        touch {output.done}
        '''

rule pro_comebin:
    input:
        fasta=scratch_dir + "01-analysis/07-assembled-metaspades/{sample}/scaffolds_1000.fa",
        map=scratch_dir + "01-analysis/13-bacmags/01-alignment/{sample}/alignment.bam"
    output:
        dir=directory(scratch_dir + "01-analysis/13-bacmags/02-comebin/{sample}/"),
        done=scratch_dir + "01-analysis/13-bacmags/02-comebin/{sample}-comebin-done.txt"
    params:
        map=scratch_dir + "01-analysis/13-bacmags/01-alignment/{sample}/"
    conda:
        "../envs/comebin.yaml"
    threads: 24
    shell:
        '''
        run_comebin.sh -a {input.fasta} -o {output.dir} -p {params.map} -t {threads} || true
        touch {output.done}
        '''

rule pro_metabat:
    input:
        fasta=scratch_dir + "01-analysis/07-assembled-metaspades/{sample}/scaffolds_1000.fa",
        map=scratch_dir + "01-analysis/13-bacmags/01-alignment/{sample}/alignment.bam"
    output:
        depth=scratch_dir + "01-analysis/13-bacmags/05-metabat2/{sample}/{sample}_depth.txt",
        dir=directory(scratch_dir + "01-analysis/13-bacmags/05-metabat2/{sample}/bins")
    params:
        bins="{sample}"
    conda:
        "../envs/metabat2.yaml"
    shell:
        '''
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.map} 
        metabat2 -i {input.fasta} -a {output.depth} -o {output.dir}/{params.bins} -t 0 --minContig 1500
        '''

rule contigtocomplete:
    input:
        completebin=scratch_dir + "01-analysis/13-bacmags/03-completebin/{sample}-final/",
        done=scratch_dir + "01-analysis/13-bacmags/03-completebin/{sample}-completebin-done.txt"
    output:
        completebin=scratch_dir + "01-analysis/13-bacmags/08-dastool/{sample}.completebin.scaffolds2bin.tsv"
    params:
        dir=scratch_dir + "01-analysis/13-bacmags/08-dastool/",
        dastool="/home/ad14556/snakemake-pipelines/metagenomicsPipeline/workflow/scripts/DAS_Tool/src"
    shell:
        '''
        mkdir -p {params.dir}
        if [[ -d {input.completebin} ]]; then
          {params.dastool}/./Fasta_to_Contig2Bin.sh -i {input.completebin} -e fasta > {output.completebin}
        else
          touch {output.completebin}
        fi
        '''

rule contigtometabat: 
    input:
        metabat=scratch_dir + "01-analysis/13-bacmags/05-metabat2/{sample}/bins/"
    output:
        metabat=scratch_dir + "01-analysis/13-bacmags/08-dastool/{sample}.metabat.scaffolds2bin.tsv"
    params:
        dir=scratch_dir + "01-analysis/13-bacmags/08-dastool/",
        dastool="/home/ad14556/snakemake-pipelines/metagenomicsPipeline/workflow/scripts/DAS_Tool/src"
    shell:
        '''
        mkdir -p {params.dir}
        if [[ -d {input.metabat} ]]; then
          {params.dastool}/./Fasta_to_Contig2Bin.sh -i {input.metabat} -e fa > {output.metabat}
        else
          touch {output.metabat}
        fi
        '''

rule contigtocomebin:
    input:
        comebin=scratch_dir + "01-analysis/13-bacmags/02-comebin/{sample}-comebin-done.txt"
    output:
        comebin=scratch_dir + "01-analysis/13-bacmags/08-dastool/{sample}.comebin.scaffolds2bin.tsv"
    params:
        dir=scratch_dir + "01-analysis/13-bacmags/08-dastool/",
        comebin=scratch_dir + "01-analysis/13-bacmags/02-comebin/{sample}/comebin_res/comebin_res_bins/",
        dastool="/home/ad14556/snakemake-pipelines/metagenomicsPipeline/workflow/scripts/DAS_Tool/src"
    shell:
        '''
        mkdir -p {params.dir}
        if [[ -d {params.comebin} ]]; then
          {params.dastool}/./Fasta_to_Contig2Bin.sh -i {params.comebin} -e fa > {output.comebin}
        else
          touch {output.comebin}
        fi
        '''

rule pro_dastool:
    input:
        contigs=scratch_dir + "01-analysis/07-assembled-metaspades/{sample}/scaffolds_1000.fa",
        completebin=scratch_dir + "01-analysis/13-bacmags/08-dastool/{sample}.completebin.scaffolds2bin.tsv",
        comebin=scratch_dir + "01-analysis/13-bacmags/08-dastool/{sample}.comebin.scaffolds2bin.tsv",
        metabat=scratch_dir + "01-analysis/13-bacmags/08-dastool/{sample}.metabat.scaffolds2bin.tsv"
    output:
        dastool=directory(scratch_dir + "01-analysis/13-bacmags/08-dastool/{sample}"),
        done=scratch_dir + "01-analysis/13-bacmags/08-dastool/{sample}-dastool-done.txt"
    params:
        basename=scratch_dir + "01-analysis/13-bacmags/08-dastool/{sample}/{sample}",
        dastool="/home/ad14556/snakemake-pipelines/metagenomicsPipeline/workflow/scripts/DAS_Tool/src"
    log:
        out=scratch_dir + "03-log/dastool/{sample}.out",
        err=scratch_dir + "03-log/dastool/{sample}.err"
    conda:
        "../envs/dastool.yaml"
    shell:
        '''
        mkdir -p {output.dastool}
        Rscript {params.dastool}/DAS_Tool.R -i {input.completebin},{input.comebin},{input.metabat} -l completebin,comebin,metabat -c {input.contigs} \
            -o {params.basename} --write_bins --write_bin_evals -t 12 --score_threshold=0 > {log.out} 2> {log.err}; touch {output.done}
        '''

rule prepare_bins:
    input:
        dastool=scratch_dir + "01-analysis/13-bacmags/08-dastool/{sample}-dastool-done.txt"
    output:
        bins=scratch_dir + "01-analysis/13-bacmags/{sample}-microbiome-bins-done.txt"
    params:
        dir=scratch_dir + "01-analysis/13-bacmags/09-bins",
        final_bins=scratch_dir + "01-analysis/13-bacmags/09-bins/{sample}_MAG_",
        table=scratch_dir + "01-analysis/13-bacmags/08-dastool/{sample}/{sample}_DASTool_summary.tsv",
        dastool_bins=scratch_dir + "01-analysis/13-bacmags/08-dastool/{sample}/{sample}_DASTool_bins",
        comp=completeness,
        cont=contamination
    shell:
        """
        mkdir -p {params.dir}

        python workflow/scripts/process_data_2.py \
            {params.table} \
            {params.final_bins} \
            --completeness_threshold {params.comp} \
            --contamination_threshold {params.cont} \
            --source_dir {params.dastool_bins}

        touch {output.bins}
        """

rule magpurify:
    input:
        bin=scratch_dir + "01-analysis/13-bacmags/09-bins/{bin}.fa",
        db=magpurify_db
    output:
        mag=scratch_dir + "01-analysis/13-bacmags/10-magpurify/{bin}.fa"
    params:
        temp=lambda wc: scratch_dir + f"01-analysis/13-bacmags/10-magpurify-temp/{wc.bin}"
    conda:
        "../envs/magpurify.yaml"
    shell:
        """
        mkdir -p "{params.temp}"
        mkdir -p "$(dirname "{output.mag}")"

        magpurify phylo-markers "{input.bin}" "{params.temp}" --db "{input.db}"
        magpurify clade-markers "{input.bin}" "{params.temp}" --db "{input.db}"
        magpurify known-contam "{input.bin}" "{params.temp}" --db "{input.db}"
        magpurify gc-content "{input.bin}" "{params.temp}" --db "{input.db}"
        magpurify clean-bin "{input.bin}" "{params.temp}" "{output.mag}"
        """

def get_magpurify_bins(wc):
    bins = [
        os.path.basename(x).replace(".fa", "")
        for x in glob.glob(
            scratch_dir + f"01-analysis/13-bacmags/09-bins/{wc.sample}_MAG_*.fa"
        )
    ]

    return expand(
        scratch_dir + "01-analysis/13-bacmags/10-magpurify/{bin}.fa",
        bin=bins
    )


rule microbiome_bins:
    input:
        mags=get_magpurify_bins
    output:
        done=scratch_dir + "01-analysis/13-bacmags/10-magpurify/{sample}-magpurify-done.txt"
    params:
        dir=output_dir + "04-microbiome-bins"
    shell:
        r"""
        mkdir -p "{params.dir}"

        cp {input.mags} "{params.dir}/"

        touch "{output.done}"
        """

#rule prepare_bins:
#    input:
#        dastool=scratch_dir + "01-analysis/13-bacmags/08-dastool/{sample}-dastool-done.txt",
#    output:
#       bins=scratch_dir + "01-analysis/13-bacmags/{sample}-microbiome-bins-done.txt",
#    params:
#        dir=scratch_dir + "01-analysis/13-bacmags/09-bins",
#        final_bins=scratch_dir + "01-analysis/13-bacmags/09-bins/{sample}_MAG_", 
#        table=scratch_dir + "01-analysis/13-bacmags/08-dastool/{sample}/{sample}_DASTool_summary.tsv",
#        dastool_bins=scratch_dir + "01-analysis/13-bacmags/08-dastool/{sample}/{sample}_DASTool_bins",
#        comp=completeness,
#        cont=contamination
#    shell:
#        '''
#	mkdir -p {params.dir}
#       python workflow/scripts/process_data_2.py {params.table} {params.final_bins} \
#            --completeness_threshold {params.comp} --contamination_threshold {params.cont} \
#            --source_dir {params.dastool_bins}
#        touch {output.bins}
#        '''
#
#rule magpurify:
#    input:
#        bin=scratch_dir + "01-analysis/13-bacmags/09-bins/{sample}.fa",
#        db=magpurify_db
#    output:
#        mag=scratch_dir + "01-analysis/13-bacmags/09-mag-purify/{sample}.fa"
#    params:
#        temp=lambda wc: scratch_dir + f"01-analysis/13-bacmags/09-mag-purify-temp/{wc.sample}"
#    conda: 
#        "../envs/magpurify.yaml"
#    shell:
#        r"""
#        mkdir -p "{params.temp}"
#
#        magpurify phylo-markers "{input.bin}" "{params.temp}" --db "{input.db}"
#        magpurify clade-markers "{input.bin}" "{params.temp}" --db "{input.db}"
#        magpurify known-contam "{input.bin}" "{params.temp}" --db "{input.db}"
#        magpurify gc-content "{input.bin}" "{params.temp}" --db "{input.db}"
#        magpurify clean-bin "{input.bin}" "{params.temp}" "{output.mag}"
#        """
#
#rule magpurify:
#    input:
#        bin=scratch_dir + "01-analysis/13-bacmags/{sample}-microbiome-bins-done.txt",
#        db=magpurify_db
#    output:
#        mag=directory(scratch_dir + "01-analysis/13-bacmags/09-mag-purify"),
#        temp=directory(scratch_dir + "01-analysis/13-bacmags/09-mag-purify-temp"),
#        done=scratch_dir + "01-analysis/13-bacmags/09-magpurify/{sample}-magpurify-done.txt"
#    params:
#        bins=scratch_dir + "01-analysis/13-bacmags/09-bins"
#    log:
#        out=scratch_dir + "03-log/magpurify/{sample}.out",
#        err=scratch_dir + "03-log/magpurify/{sample}.err"
#    conda:
#        "../envs/magpurify.yaml"
#    shell:
#        '''
#        for file in {params.bins}/*; do
#          sample=$(basename ${file} .fa)
#          mkdir {output.temp}/${sample}
#          magpurify phylo-markers ${file} {output.temp}/${sample} --db {input.db}
#          magpurify clade-markers ${file} {output.temp}/${sample} --db {input.db}	  
#          magpurify known-contam ${file} {output.temp}/${sample} --db {input.db}
#	  magpurify gc-content ${file} {output.temp}/${sample} --db {input.db}
#          magpurify clean-bin ${file} {output.temp}/${sample} {output.mag}/${sample}.fa
#        done
#        touch {output.done}
#        '''
#
#rule microbiome_bins:
#    input:
#        magpurify_done=scratch_dir + "01-analysis/13-bacmags/08-magpurify/{sample}-magpurify-done.txt",
#        magpurify=directory(scratch_dir + "01-analysis/13-bacmags/09-mag-purify")
#    output:
#        done=scratch_dir + "01-analysis/13-bacmags/09-magpurify/{sample}-magpurify-done.txt"
#    params:
#        dir=output_dir + "04-microbiome-bins"
#    shell:
#        '''
#	touch {input.magpurify}/* && cp {input.magpurify}/* {output.dir}
#        '''
