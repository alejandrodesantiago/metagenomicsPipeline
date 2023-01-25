rule flash:
    input:
        R1=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz",
        R2=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz"
    output:
        merged=scratch_dir + "01-analysis/06-flash/{sample}_merged.fastq.gz"
    params:
    conda:
        ".../envs/flash.yaml"
    shell:
        '''
        flash {input.R1} {input.R2} {output.merged}
        '''

rule metaphlan:
    input:
        reads=expand(scratch_dir + "01-analysis/06-flash/{sample}_merged.fastq.gz", sample=SAMPLES)
    output:
        profile=scratch + "01-analysis/07-metaphlan/{sample}_taxonomy_profile.txt"
    params:
    conda:
        ".../envs/metaphlan.yaml"
    shell:
        '''
        metaphlan {input.profile} --nproc 5 --input_type fastq -o {output.profile}
        '''
