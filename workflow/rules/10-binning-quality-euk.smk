rule busco:
    input:
        mags=scratch_dir + "01-analysis/14-eukmags/04-dastool/{sample}.bin.{number}.fa"
    output:
        busco=directory(scratch_dir + "01-analysis/14-eukmags/16-busco/{sample}")
    wildcard_constraints:
        number="[0-9]+"
    params:
        lineage="nematoda_odb10", # options are bac, arc, euk, mito
        sample="sample_bin"
    conda:
        "../envs/busco.yaml"
    shell:
        '''
        busco -i {input.mags} -l {params.lineage} -o {params.sample} --out_path {output.busco} -m genome -c 4
        '''

#rule eukcc:
#    input:
#        mags=scratch_dir + "01-analysis/14-eukmags/04-dastool/{sample}.bin.{number}.fa"
#    output:
#        eukcc=directory(scratch_dir + "01-analysis/14-eukmags/16-eukmag_quality/eucc/{sample}.bin.{number}")
#    wildcard_constraints:
#        number="[0-9]+"
#    params:
#        threads="4"
#    conda:
#        "../envs/eukcc.yaml"
#    shell:
#        '''
#        eukcc single --out {output.eukcc} --threads {params.threads} {input.mag}
#        '''