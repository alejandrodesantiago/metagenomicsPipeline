rule pro_mapReads:
    input:
        euk=scratch_dir + "01-analysis/13-eukrep/prokaryotes.fasta",
        r1_paired=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz",
        r2_paired=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz",
        r1_unpaired=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_unpaired.fastq.gz",
        r2_unpaired=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_unpaired.fastq.gz"
    output:
        final=scratch_dir + "01-analysis/15-bacmags/01-alignment/{sample}_final.bam",
        mapped_paired=scratch_dir + "01-analysis/15-bacmags/01-alignment/{sample}_paired_reads.bam",
        mapped_unpaired=scratch_dir + "01-analysis/15-bacmags/01-alignment/{sample}_unpaired_reads.bam",
        unpaired_reads=scratch_dir + "01-analysis/03-trimmomatic/{sample}_single_reads.fastq.gz"
    params:
        threads=2
    conda:
        "../envs/bwa-samtools.yaml"
    shell:
        '''
        cat {input.r1_unpaired} {input.r2_unpaired} > {output.unpaired_reads}
        bwa index {input.euk}
        bwa mem -t {params.threads} {input.euk} {input.r1_paired} {input.r2_paired} | samtools sort -o {output.mapped_paired}
        bwa mem -t {params.threads} {input.euk} {output.unpaired_reads} | samtools sort -o {output.mapped_unpaired}
        samtools merge {output.final} {output.mapped_paired} {output.mapped_unpaired}
        samtools index {output.final}
        '''

rule pro_metabat:
    input:
        fasta=scratch_dir + "01-analysis/13-eukrep/prokaryotes.fasta",
        bam=scratch_dir + "01-analysis/15-bacmags/01-alignment/{sample}_final.bam"
    output:
        depth=scratch_dir + "01-analysis/15-bacmags/02-metabat2/depth.txt"
    params:
        bin=scratch_dir + "01-analysis/15-bacmags/02-metabat2/bins"
    conda:
        "../envs/metabat2.yaml"
    shell:
        '''
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}
        metabat2 -i {input.fasta} -a {output.depth} -o {params.bin}
        '''

rule pro_concoct:
    input:
        fasta=scratch_dir + "01-analysis/13-eukrep/prokaryotes/{sample}_pro.fasta",
        bam=scratch_dir + "01-analysis/15-bacmags/01-alignment/{sample}_final.bam"
    output:
        bed=scratch_dir + "01-analysis/15-bacmags/03-concoct/{sample}_concoct.bed",
        fasta=scratch_dir + "01-analysis/15-bacmags/03-concoct/{sample}_concoct.fasta",
        depth=scratch_dir + "01-analysis/15-bacmags/03-concoct/{sample}_depth.tsv",
        bin=directory(scratch_dir + "01-analysis/15-bacmags/03-concoct/{sample}/bins/")
    params:
        dir = directory(scratch_dir + "01-analysis/15-bacmags/03-concoct/{sample}"),
    conda:
        "../envs/concoct.yaml"
    shell:
        '''
        cut_up_fasta.py {input.fasta} -c 10000 -o 0 --merge_last -b {output.bed} > {output.fasta}
        concoct_coverage_table.py {output.bed} {input.bam} > {output.depth}
        concoct --composition_file {output.fasta} --coverage_file {output.depth} -b {params.dir}
        merge_cutup_clustering.py {params.dir}/clustering_gt1000.csv > {params.dir}/clustering_merged.csv
        mkdir -p {output.bin}
        extract_fasta_bins.py {input.fasta} {params.dir}/clustering_merged.csv --output_path {output.bin}
        '''

rule pro_dastool:
    input:
        contigs=scratch_dir + "01-analysis/13-eukrep/prokaryotes/{sample}_pro.fasta",
        concoct=scratch_dir + "01-analysis/15-bacmags/03-concoct/{sample}/bins/",
        metabat=scratch_dir + "01-analysis/15-bacmags/02-metabat2/{sample}_depth.txt"
    output:
        metabat=scratch_dir + "01-analysis/15-bacmags/04-dastool/{sample}.metabat.scaffolds2bin.tsv",
        concoct=scratch_dir + "01-analysis/15-bacmags/04-dastool/{sample}.concoct.scaffolds2bin.tsv",
        dastool=scratch_dir + "01-analysis/15-bacmags/04-dastool/{sample}"
    params:
        metabat=scratch_dir + "01-analysis/15-bacmags/02-metabat2/{sample}_bin",
        basename="{sample}",
        dastool="/home/ad14556/snakemake-pipelines/metagenomicsPipeline/workflow/scripts/DAS_Tool/src"
    conda:
        "../envs/dastool.yaml"
    shell:
        '''
        mkdir -p {output.dastool}
        mkdir -p {params.metabat}  
        mv {params.metabat}*.fa {params.metabat} || echo "(Error ok if no files found)"
        {params.dastool}/./Fasta_to_Contig2Bin.sh -i {params.metabat} -e fa > {output.metabat}
        {params.dastool}/./Fasta_to_Contig2Bin.sh -i {input.concoct} -e fa > {output.concoct}
        Rscript {params.dastool}/DAS_Tool.R -i {output.metabat},{output.concoct} -l metabat,concoct -c {input.contigs} -o {output.dastool}/{params.basename} --write_bins --write_bin_evals -t 4
        '''


