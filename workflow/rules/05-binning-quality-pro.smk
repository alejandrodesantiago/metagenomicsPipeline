rule pro_gtdbtk:
    input:
#        bins=output_dir + "04-microbiome-bins/",
        database=gtdbtk_db
    output:
        gtdb=directory(scratch_dir + "01-analysis/17-gtdb/"),
        done=scratch_dir + "01-analysis/final-bins-gtdbtk-done.txt"
    params:
        bins=output_dir + "04-microbiome-bins/"
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        '''
        conda env config vars set GTDBTK_DATA_PATH={input.database}
        gtdbtk classify_wf --genome_dir {params.bins} --out_dir {output.gtdb} -x fa --cpus 24 --pplacer_cpus 24 --skip_ani_screen
        touch {output.done}
        '''

rule checkm_2:
    input:
#        bins=output_dir + "04-microbiome-bins/",
        database=checkm_db
    output:
        checkm=directory(scratch_dir + "01-analysis/18-checkm"),
        done=scratch_dir + "01-analysis/final-bins-checkm-done.txt"
    params:
        bins=output_dir + "04-microbiome-bins/" 
    log:
        err=scratch_dir + "03-log/18-checkm.err",
        out=scratch_dir + "03-log/18-checkm.out"
#    conda:
#        "../envs/checkm2.yaml"
    shell:
        '''
        module load CheckM2/1.1.0-foss-2024a
        checkm2 predict --force -x .fa --threads 12 --database_path {input.database} \
          --input {params.bins} --output-directory {output.checkm} > {log.out} 2> {log.err}
        touch {output.done}
        '''

rule microbiome_bins_quality:
    input:
        gtdbtk_done=scratch_dir + "01-analysis/final-bins-gtdbtk-done.txt",
        checkm_done=scratch_dir + "01-analysis/final-bins-checkm-done.txt",
        checkm=scratch_dir + "01-analysis/18-checkm",
        gtdbtk=scratch_dir + "01-analysis/17-gtdb",
    output:
        checkm=output_dir + "05-microbiome-bins-quality.tsv",
        bac=output_dir + "05-microbiome-bins-bac-lineages.tsv",
        arc=output_dir + "05-microbiome-bins-arc-lineages.tsv"
    params:
        checkm="quality_report.tsv",
        arc="ar53.summary.tsv",
        bac="bac120.summary.tsv"
    shell:
        '''
        touch {input.checkm}/{params.checkm} && cp {input.checkm}/{params.checkm} {output.checkm}
        touch {input.gtdbtk}/*{params.bac} && cp {input.gtdbtk}/*{params.bac} {output.bac}
        touch {input.gtdbtk}/*{params.arc} && cp {input.gtdbtk}/*{params.arc} {output.arc}
        '''

