rule pro_mapReads:
    input:
        fasta=scratch_dir + "01-analysis/08-assembled-megahit/{sample}",
        r1_paired=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz",
        r2_paired=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz"
    output:
        map=directory(scratch_dir + "01-analysis/13-bacmags/01-alignment/{sample}")
    params:
        map="alignment.bam",
        threads=48,
        file="final.contigs.fa"
    conda:
        "../envs/bwa-samtools.yaml"
    shell:
        '''
        mkdir {output.map}
        bwa index {input.fasta}/{params.file}
        bwa mem -t {params.threads} {input.fasta}/{params.file} {input.r1_paired} {input.r2_paired} | samtools sort -o {output.map}/{params.map} --threads {params.threads}
        samtools index -@ {params.threads} {output.map}/{params.map}
        '''

rule pro_comebin:
    input:
        fasta=scratch_dir + "01-analysis/08-assembled-megahit/{sample}",
        map=scratch_dir + "01-analysis/13-bacmags/01-alignment/{sample}"
    output:
        dir=directory(scratch_dir + "01-analysis/13-bacmags/02-comebin/{sample}/")
    params:
        bins="{sample}",
        file="final.contigs.fa"
    conda:
        "../envs/comebin.yaml"
    shell:
        '''
        run_comebin.sh -a {input.fasta}/{params.file} -o {output.dir} -p {input.map} -t 40
        '''

rule pro_metabat:
    input:
        fasta=scratch_dir + "01-analysis/08-assembled-megahit/{sample}",
        map=scratch_dir + "01-analysis/13-bacmags/01-alignment/{sample}"
    output:
        depth=scratch_dir + "01-analysis/13-bacmags/03-metabat2/{sample}/{sample}_depth.txt",
        dir=directory(scratch_dir + "01-analysis/13-bacmags/03-metabat2/{sample}/bins")
    params:
        bins="{sample}",
        file="final.contigs.fa",
	map="alignment.bam"
    conda:
        "../envs/metabat2.yaml"
    shell:
        '''
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.map}/{params.map}
        metabat2 -i {input.fasta}/{params.file} -a {output.depth} -o {output.dir}/{params.bins}
        '''

rule pro_semibin:
    input:
        fasta=scratch_dir + "01-analysis/08-assembled-megahit/{sample}",
        map=scratch_dir + "01-analysis/13-bacmags/01-alignment/{sample}"
    output:
        dir=directory(scratch_dir + "01-analysis/13-bacmags/04-semibin/{sample}/")
    params:
        bins="{sample}",
        file="final.contigs.fa",
        map="alignment.bam"
    conda:
        "../envs/semibin2.yaml"
    shell:
        '''
        SemiBin2 single_easy_bin --input-fasta {input.fasta}/{params.file} --input-bam {input.map}/{params.map} --environment ocean --output {output.dir} --threads 40 --compression none
        '''

rule pro_dastool:
    input:
        contigs=scratch_dir + "01-analysis/08-assembled-megahit/{sample}/final.contigs.fa",
        comebin=scratch_dir + "01-analysis/13-bacmags/02-comebin/{sample}/",
        metabat=scratch_dir + "01-analysis/13-bacmags/03-metabat2/{sample}/bins/",
        semibin=scratch_dir + "01-analysis/13-bacmags/04-semibin/{sample}/"
    output:
        metabat=scratch_dir + "01-analysis/13-bacmags/05-dastool/{sample}.metabat.scaffolds2bin.tsv",
        comebin=scratch_dir + "01-analysis/13-bacmags/05-dastool/{sample}.comebin.scaffolds2bin.tsv",
        semibin=scratch_dir + "01-analysis/13-bacmags/05-dastool/{sample}.semibin.scaffolds2bin.tsv",
        dastool=directory(scratch_dir + "01-analysis/13-bacmags/05-dastool/{sample}")
    params:
        basename="{sample}",
        comebin="comebin_res/comebin_res_bins/",
        semibin="output_bins/",
        dastool="/home/ad14556/snakemake-pipelines/metagenomicsPipeline/workflow/scripts/DAS_Tool/src"
    conda:
        "../envs/dastool.yaml"
    shell:
        '''
        mkdir -p {output.dastool}
        {params.dastool}/./Fasta_to_Contig2Bin.sh -i {input.metabat} -e fa > {output.metabat}
        {params.dastool}/./Fasta_to_Contig2Bin.sh -i {input.comebin}/{params.comebin} -e fa > {output.comebin}
        {params.dastool}/./Fasta_to_Contig2Bin.sh -i {input.semibin}/{params.semibin} -e fa > {output.semibin}
        Rscript {params.dastool}/DAS_Tool.R -i {output.metabat},{output.comebin},{output.semibin} -l metabat,comebin,semibin -c {input.contigs} -o {output.dastool}/{params.basename} --write_bins --write_bin_evals -t 12 --score_threshold=0
        '''


