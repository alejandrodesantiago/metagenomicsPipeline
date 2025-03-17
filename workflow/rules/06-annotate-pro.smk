rule dram_db:
    input:
    output:
        database=directory(scratch_dir + "02-databases/dram")
    conda:
        "../envs/dram.yaml"
    shell:
        '''
        DRAM-setup.py prepare_databases --output_dir {output.database}
        '''

rule dram_annotate:
    input:
        mags=directory(scratch_dir + "01-analysis/14-bacmags/04-dastool/{sample}")
        database=scratch_dir + "02-databases/dram"
    output:
        annotations=scratch_dir + "01-analysis/16-bacmags-dram/annotations",
        summary=scratch_dir + "01-analysis/16-bacmags-dram/mag-summaries"
    database=directory()
    conda:
           "../envs/dram.yaml
    shell:
        '''
        DRAM.py annotate -i {input.mags} -o {output.annotations}
        DRAM.py distill -i {output.annotations} -o {output.summary} --trna_path {output.annotations}/trnas.tsv --rrna_path {output.annotations}/rrnas.tsv
        '''
