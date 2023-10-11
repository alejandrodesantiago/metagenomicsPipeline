rule mitoz:
    input:
        dir=scratch_dir + "01-analysis/10-assembled-megahit"
    output:
        dir=directory(scratch_dir + "01-analysis/mitoz")
    params:
        contigfile="final.contigs.fa",
        prefix="nem_mito"
 #   conda:
 #       "../envs/mitoz.yaml"
    shell:
        '''
        module load numpy
        mitoz findmitoscaf --fastafile {input.dir}/{params.contigfile} --workdir {output.dir} --outprefix {params.prefix} --min_abundance 0 --clade Nematoda --requiring_taxa Nematoda
        '''