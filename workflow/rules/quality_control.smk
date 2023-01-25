rule initial_fastqc:
    input:
        input_dir + "{file}.fastq.gz"
    output:
        html=scratch_dir + "01-analysis/01-initial-fastqc/{file}.html",
        zip=scratch_dir + "01-analysis/01-initial-fastqc/{file}_fastqc.zip"
    params:
    log:
        scratch_dir + "03-log/01-initial-fastqc/{file}.log"
    wrapper:
        "v1.5.0/bio/fastqc"

rule initial_multiqc:
    input:
        expand(scratch_dir + "01-analysis/01-initial-fastqc/{file}_fastqc.zip", file=FILES)
    output:
        scratch_dir + "01-analysis/02-initial-multiqc/multiqc.html"
    params:
    log:
        scratch_dir + "03-log/02-initial-multiqc.log"
    wrapper:
        "v1.5.0/bio/multiqc"

rule trimmomatic:
    input:
        r1=input_dir + "{sample}.R1.fastq.gz",
        r2=input_dir + "{sample}.R2.fastq.gz"
    output:
        r1=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_paired.fastq.gz",
        r2=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_paired.fastq.gz",
        r1_unpaired=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R1_unpaired.fastq.gz",
        r2_unpaired=scratch_dir + "01-analysis/03-trimmomatic/{sample}_R2_unpaired.fastq.gz"
    params:
        trimmer=["ILLUMINACLIP:{}:2:30:7".format(adapters), "LEADING:2", "TRAILING:2", "SLIDINGWINDOW:4:2$
        extra="",
        compression_level="-9"
    log:
        scratch_dir + "03-log/03-trimmomatic/{sample}.log"
    threads:
        32
    wrapper:
        "v1.5.0/bio/trimmomatic/pe"

rule trimmed_fastqc:
    input:
        scratch_dir + "01-analysis/03-trimmomatic/{sample}_{dir}_{pair}.fastq.gz"
    output:
        html=scratch_dir + "01-analysis/04-trimmed-fastqc/{sample}_{dir}_{pair}_fastqc.html",
        zip=scratch_dir + "01-analysis/04-trimmed-fastqc/{sample}_{dir}_{pair}_fastqc.zip"
    params: ""
    log:
        scratch_dir + "03-log/04-trimmed-fastqc/{sample}_{dir}_{pair}.log"
    threads:
        2
    wrapper:
        "v1.5.0/bio/fastqc"

rule trimmed_multiqc:
    input:
        expand(scratch_dir + "01-analysis/04-trimmed-fastqc/{sample}_{dir}_{pair}_fastqc.zip", sample=SAM$
    output:
        scratch_dir + "01-analysis/05-trimmed-multiqc/multiqc.html"
    params: ""
    log:
        scratch_dir + "03-log/05-initial-multiqc.log"
    wrapper:
        "v1.21.2/bio/multiqc"
