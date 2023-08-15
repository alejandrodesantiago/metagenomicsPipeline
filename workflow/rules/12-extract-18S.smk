rule euk_barrnap:
    input:
        mags=scratch_dir + "01-analysis/14-eukmags/04-dastool/{sample}.bin.{number}.fa"
    output:
        fasta=scratch_dir + "01-analysis/14-eukmags/18-barrnap/{sample}_18S.bin.{number}.fasta",
        gff=scratch_dir + "01-analysis/14-eukmags/18-barrnap/{sample}_18S.bin.{number}.gff"
    wildcard_constraints:
        number="[0-9]+"
    params:
        kingdom="euk" # options are bac, arc, euk, mito
    conda:
        "../envs/barrnap.yaml"
    shell:
        '''
        barrnap -k euk -o {output.fasta} < {input.mags} > {output.gff}
        '''