rule initial_fastqc:
    input:
        input_dir + "{file}.fastq.gz"
    output:
        html=scratch_dir + "01-analysis/01-initial-fastqc/{file}.html",
        zip=scratch_dir + "01-analysis/01-initial-fastqc/{file}_fastqc.zip"
    log:
        scratch_dir + "03-log/01-initial-fastqc/{file}.log"
    wrapper:
        "v1.5.0/bio/fastqc"

rule initial_multiqc:
    input:
        expand(scratch_dir + "01-analysis/01-initial-fastqc/{file}_fastqc.zip", file=FILES)
    output:
        scratch_dir + "01-analysis/02-initial-multiqc/multiqc.html"
    log:
        scratch_dir + "03-log/02-initial-multiqc.log"
    wrapper:
        "v1.5.0/bio/multiqc"

rule trimmomatic:
    input:
        r1=input_dir + "{sample}_R1.fastq.gz",
        r2=input_dir + "{sample}_R2.fastq.gz"
    output:
        r1=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz",
        r2=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz",
        r1_unpaired=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_unpaired.fastq.gz",
        r2_unpaired=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_unpaired.fastq.gz"
    params:
        trimmer=["ILLUMINACLIP:{}:2:30:7".format(adapters), "LEADING:2", "TRAILING:2", "SLIDINGWINDOW:4:20", "MINLEN:55"],
        extra="",
        compression_level="-9"
    log:
        scratch_dir + "03-log/03-trimmomatic/{sample}.log"
    threads:
        32
    wrapper:
        "v1.5.0/bio/trimmomatic/pe"

rule fastp:
    input:
        r1=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz",
        r2=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz"
    output:
        r1=scratch_dir + "01-analysis/04-deduplicated/{sample}_R1_paired.fastq.gz",
        r2=scratch_dir + "01-analysis/04-deduplicated/{sample}_R2_paired.fastq.gz"
    log:
        scratch_dir + "03-log/04-deduplicated/{sample}_deduplicated.log"
    shell:
        '''
        module load fastp
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --dedup --dup_calc_accuracy 3
        '''

rule trimmed_fastqc:
    input:
        scratch_dir + "01-analysis/04-deduplicated/{sample}_{dir}_paired.fastq.gz"
    output:
        html=scratch_dir + "01-analysis/05-trimmed-fastqc/{sample}_{dir}_paired_fastqc.html",
        zip=scratch_dir + "01-analysis/05-trimmed-fastqc/{sample}_{dir}_paired_fastqc.zip"
    log:
        scratch_dir + "03-log/05-trimmed-fastqc/{sample}_{dir}_paired.log"
    threads:
        2
    wrapper:
        "v1.5.0/bio/fastqc"

rule trimmed_multiqc:
    input:
        expand(scratch_dir + "01-analysis/05-trimmed-fastqc/{sample}_{dir}_paired_fastqc.zip", sample=SAMPLES, dir=["R1", "R2"])
    output:
        scratch_dir + "01-analysis/06-trimmed-multiqc/multiqc.html"
    log:
        scratch_dir + "03-log/06-initial-multiqc.log"
    wrapper:
        "v1.21.2/bio/multiqc"

#rule process_trimming:
#    input:
#        expand(scratch_dir + "01-analysis/04-trimmed-fastqc/{sample}_{dir}_paired_fastqc.zip", sample=SAMPLES, dir=["R1", "R2"])
#    output:
#        done=scratch_dir + "01-analysis/06-all-samples-quality-trimmed.txt"
#    shell:
#        "echo 'All samples have been quality controlled' > {output.done}"

rule reads_quality_trim:
    input:
#        done=scratch_dir + "01-analysis/06-all-samples-quality-trimmed.txt",
        qc_initial=scratch_dir + "01-analysis/02-initial-multiqc/multiqc.html",
        qc_trimmed=scratch_dir + "01-analysis/06-trimmed-multiqc/multiqc.html"
    output:
        trimmed=directory(output_dir + "01-reads-quality-trimmed/"),
        qc_initial=output_dir + "01-reads-quality-trimmed-multiqc-initial.html",
        qc_trimmed=output_dir + "01-reads-quality-trimmed-multiqc-trimmed.html"
    params:
        trimmed=scratch_dir + "01-analysis/04-deduplicated"
    shell:
        '''
        mkdir -p {output.trimmed}
        cp {params.trimmed}/*_paired* {output.trimmed}/
        cp {input.qc_initial} {output.qc_initial}
        cp {input.qc_trimmed} {output.qc_trimmed}
        '''
