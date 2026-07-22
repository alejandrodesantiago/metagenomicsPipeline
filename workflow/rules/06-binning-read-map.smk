rule extract_bac:
    input:
        db=kraken_db,
        R1=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz",
        R2=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz"
    output:
        R1=scratch_dir + "01-analysis/19-kraken/{sample}_R_1.fastq",
        R2=scratch_dir + "01-analysis/19-kraken/{sample}_R_2.fastq",
        done=scratch_dir + "01-analysis/19-kraken-{sample}-done.txt"
    log:
        err=scratch_dir + "03-logs/19-kraken/19-{sample}-kraken.err",
        out=scratch_dir + "03-logs/19-kraken/19-{sample}-kraken.err"
    params:
        base=scratch_dir + "01-analysis/19-kraken/{sample}"
    conda:
        "../envs/kraken2.yaml"
    shell:
        '''
        kraken2 --db /db/kraken2/20240906/pluspf --paired --classified-out {params.base}_R#.fastq \
            {input.R1} {input.R2} --threads 24 > {log.out} 2> {log.err}
        touch {output.done}
        '''

rule all_samples_processed:
    input:
        done=expand(scratch_dir + "01-analysis/19-kraken-{sample}-done.txt", sample=SAMPLES)
    output:
        done=scratch_dir + "01-analysis/23-all-samples-finished-processing.txt"
    shell:
        "echo 'All samples have been processed. Running RRAP...' > {output}"

rule create_map:
    input:
        rule=scratch_dir + "01-analysis/23-all-samples-finished-processing.txt",
    params:
        dir=scratch_dir + "01-analysis/19-kraken/"
    output:
        file=scratch_dir + "01-analysis/24-rrap-map.txt"
    shell:
        '''
        echo {params.dir} > {output.file}
        '''

rule run_rrap_mags:
    input:
        map=scratch_dir + "01-analysis/24-rrap-map.txt",
#        bac=output_dir + "04-microbiome-bins/"
    output:
        log=scratch_dir + "01-analysis/25-rrap-mags-results/rpkm/MeioBiome/MeioBiome_rpkm_log10.csv",
        raw=scratch_dir + "01-analysis/25-rrap-mags-results/rpkm/MeioBiome/MeioBiome_rpkm_noLog.csv"
    params:
        bins=output_dir + "04-microbiome-bins/",
        dir=scratch_dir + "01-analysis/25-rrap-mags-results/",
        name="MeioBiome",
        suffix="_R_1.fastq"
    conda:
        "../envs/rrap.yaml"
    shell:
        '''
        rrap -i {input.map} -rg {params.bins} -o {params.dir} -n {params.name} -suffix {params.suffix}
        '''

rule microbiome_bins_quant:
    input:
        rrap_mags_nolog=scratch_dir + "01-analysis/25-rrap-mags-results/rpkm/MeioBiome/MeioBiome_rpkm_noLog.csv",
        rrap_mags_log=scratch_dir + "01-analysis/25-rrap-mags-results/rpkm/MeioBiome/MeioBiome_rpkm_log10.csv"
    output:
        rrap_mags_nolog=output_dir + "06-microbiome-bins-quant-rpkm-no-log.csv",
        rrap_mags_log=output_dir + "06-microbiome-bins-quant-rpkm-log.csv"
    shell:
        '''
        cp {input.rrap_mags_nolog} {output.rrap_mags_nolog}
        cp {input.rrap_mags_log} {output.rrap_mags_log}
        '''
