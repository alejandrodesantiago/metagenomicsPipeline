# annotate eukaryote genomes using metaeuk 
rule metaeukDB:
    input:
    output:
        database=scratch_dir + "02-databases/metaeuk/MERC_MMETSP_UNICLUST50.tar.gz"
    params:
    conda:
        "../envs/metaeuk.yaml"
    shell:
        '''
        mkdir - p "/scratch/ad14556/nemMetagenomes/02-databases/metaeuk"
        wget -O https://wwwuser.gwdg.de/~compbiol/metaeuk/2019_11/MERC_MMETSP_Uniclust50_profiles.tar.gz -O {output.database}
        '''
