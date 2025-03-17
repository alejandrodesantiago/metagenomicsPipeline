rule pro_gtdbtk:
    input:
        dastool=scratch_dir + "01-analysis/13-bacmags/05-dastool/{sample}"
    output:
        gtdb=directory(scratch_dir + "01-analysis/14-gtdb/{sample}")
    params:
        bins=scratch_dir + "01-analysis/13-bacmags/05-dastool/{sample}/{sample}_DASTool_bins"
    shell:
        '''
        module load GTDB-Tk
        gtdbtk classify_wf --genome_dir {params.bins} --out_dir {output.gtdb} --skip_ani_screen -x fa --cpus 24 --pplacer_cpus 24
        '''

rule checkm:
    input:
        dastool=scratch_dir + "01-analysis/13-bacmags/05-dastool/{sample}",
        database_temp="/work/hmblab/thiosymbion-comparative-genomics/database/CheckM2_database/uniref100.KO.1.dmnd"
    output:
        checkm=directory(scratch_dir + "01-analysis/15-checkm/{sample}")
    params:
        bins=scratch_dir + "01-analysis/13-bacmags/05-dastool/{sample}/{sample}_DASTool_bins"
    shell:
        '''
        module load CheckM2
        checkm2 predict --force -x .fa --threads 12 --database_path {input.database_temp} --input {params.bins} --output-directory {output.checkm}
        '''
