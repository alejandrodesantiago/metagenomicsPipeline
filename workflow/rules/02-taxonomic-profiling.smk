rule metaphlanDB:
    output:
        database=directory(scratch_dir + "02-databases/metaphlan/")
    conda:
        "../envs/metaphlan.yaml"
    shell:
        '''
        metaphlan --install --bowtie2db {output.database}
        '''

rule metaphlan:
    input:
        R1=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz",
        R2=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz",
        database=scratch_dir + "02-databases/metaphlan/"
    output:
        profile=scratch_dir + "01-analysis/07-metaphlan/{sample}_taxonomy_profile.txt",
        browtie=scratch_dir + "01-analysis/07-metaphlan/{sample}_bowtie2.bz2"
    params:
        out=scratch_dir + "01-analysis/07-metaphlan/"
    conda:
        "../envs/metaphlan.yaml"
    shell:
        '''
        mkdir -p {params.out}
        metaphlan {input.R1},{input.R2} --bowtie2out {output.browtie} --bowtie2db {input.database} --ignore_eukaryotes --nproc 12 --input_type fastq -o {output.profile} 
        '''

rule riboDetector:
    input:
    output:
    params:
    shell:
        '''
        '''
