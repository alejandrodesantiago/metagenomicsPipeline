rule euk_barrnap:
    input:
        mags=scratch_dir + "01-analysis/14-eukmags/04-dastool/{sample}_bin.1.fa}"
    output:
        fasta=scratch_dir + "01-analysis/14-eukmags/18-barrnap/{sample}_18S.bin.1.fasta",
        gff=scratch_dir + "01-analysis/14-eukmags/18-barrnap/{sample}_18S.bin.1.gff"
    params:
        kingdom="euk" # options are bac, arc, euk, mito
    conda:
        "../envs/barrnap.yaml"
    shell:
        '''
        barrnap -k euk -o {output.fasta} < {input.mags} > {output.gff}
        '''

    barrnap-o rrna.fa < contigs.fa > rrna.gff