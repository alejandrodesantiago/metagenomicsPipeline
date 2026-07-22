rule metaspades:
    input:
        R1=scratch_dir + "01-analysis/04-deduplicated/{sample}_R1_paired.fastq.gz",
        R2=scratch_dir + "01-analysis/04-deduplicated/{sample}_R2_paired.fastq.gz"
    output:
        contigs=scratch_dir + "01-analysis/07-assembled-metaspades/{sample}/scaffolds.fasta",
        file=scratch_dir + "01-analysis/07-assembled-metaspades/{sample}-metagenome-assembly-done.txt"
    params:
        dir=scratch_dir + "01-analysis/07-assembled-metaspades/{sample}"
    threads: config['metaspades_threads']
    resources:
        mem_mb=config['metaspades_mem']
    log:
        err=scratch_dir + "03-log/07-metaspades/{sample}-metaspades.err",
        out=scratch_dir + "03-log/07-metaspades/{sample}-metaspades.out"
    conda:
        "../envs/spades.yaml"
    shell:
        '''
        export OMP_NUM_THREADS={threads}
        spades.py -1 {input.R1} -2 {input.R2} --meta --only-assembler -t {threads} -m {resources.mem_mb} -o {params.dir} > {log.out} 2> {log.err}
        touch {output.file}
        '''

rule filter_short:
    input:
        file=scratch_dir + "01-analysis/07-assembled-metaspades/{sample}-metagenome-assembly-done.txt",
        fasta=scratch_dir + "01-analysis/07-assembled-metaspades/{sample}/scaffolds.fasta"
    output:
        fasta=scratch_dir + "01-analysis/07-assembled-metaspades/{sample}/scaffolds_1000.fa"
    conda:
        "../envs/comebin.yaml"
    shell:
        '''
        python3 workflow/scripts/filter_too_short_comebin.py {input.fasta} 1000
        '''

rule metagenome_assembly:
    input:
        contigs=scratch_dir + "01-analysis/07-assembled-metaspades/{sample}/scaffolds_1000.fa",
        file=scratch_dir + "01-analysis/07-assembled-metaspades/{sample}-metagenome-assembly-done.txt"
    output:
        file=output_dir + "02-metagenome-assembly/{sample}-metagenome.fa"
    shell:
        '''
        cp {input.contigs} {output.file}
        '''
