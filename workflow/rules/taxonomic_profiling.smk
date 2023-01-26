rule flash:
    input:
        R1=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz",
        R2=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz"
    output:
        extendedFrags=scratch_dir + "01-analysis/06-flash/{sample}/{sample}.extendedFrags.fastq.gz",
        notCombined_1=scratch_dir + "01-analysis/06-flash/{sample}/{sample}.notCombined_1.fastq.gz",
        notCombined_2=scratch_dir + "01-analysis/06-flash/{sample}/{sample}.notCombined_2.fastq.gz",
        histogram=scratch_dir + "01-analysis/06-flash/{sample}/{sample}.hist",
        visualHistogram=scratch_dir + "01-analysis/06-flash/{sample}/{sample}.histogram",
        dir=directory(scratch_dir + "01-analysis/06-flash/{sample}")
    params:
        prefix="{sample}"
    conda:
        "../envs/flash.yaml"
    shell:
        '''
        flash2 --compress --max-overlap 125 --output-directory {output.dir} --output-prefix {params.prefix} {input.R1} {input.R2}
        '''

rule taxonomyDB:
    input:
    output:
        database=scratch_dir + "02-databases/metaphlan/"
    params:
    conda:
        "../envs/metaphlan.yaml"
    shell:
        '''
        metaphlan --install --bowtie2db {output.database}
        '''

rule metaphlan:
    input:
        extendedFrags=scratch_dir + "01-analysis/06-flash/{sample}/{sample}.extendedFrags.fastq.gz",
        database=scratch_dir + "02-databases/metaphlan/"
    output:
        profile=scratch_dir + "01-analysis/07-metaphlan/{sample}_taxonomy_profile.txt"
    params:
    conda:
        "../envs/metaphlan.yaml"
    shell:
        '''
        metaphlan {input.extendedFrags} --bowtie2db {input.database} --ignore-eukaryotes --nproc 5 --input_type fastq -o {output.profile} 
        '''
