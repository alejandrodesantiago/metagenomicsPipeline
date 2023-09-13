rule eukrep:
    input:
        contig=scratch_dir + "01-analysis/10-assembled-megahit/"
    output:
       euk=scratch_dir + "01-analysis/13-eukrep/eukaryotes.fasta",
       pro=scratch_dir + "01-analysis/13-eukrep/prokaryotes.fasta"
    params:
       euk=scratch_dir + "01-analysis/13-eukrep",
       pro=scratch_dir + "01-analysis/13-eukrep",
       min_contig = 5000, # set according to N50/L50
       file="final.contigs.fa"
    conda:
        "../envs/eukrep.yaml"
    shell:
        '''
        mkdir -p {params.euk}
        mkdir -p {params.pro}
        EukRep -i {input.contig}/{params.file} -o {output.euk} --prokarya {output.pro} --min {params.min_contig}
        '''
