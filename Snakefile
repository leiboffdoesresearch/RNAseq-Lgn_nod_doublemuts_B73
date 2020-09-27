configfile: "config.yaml"

container: "docker://continuumio/miniconda3:4.4.10"

rule all:
    input:
        expand('featureCounts/{sample}.txt', sample=config["samples"])

rule fastp:
    input:
        R1 = 'compressed_reads/{sample}_R1.fastq.gz',
        R2 = 'compressed_reads/{sample}_R2.fastq.gz',
    output:
        R1 = 'trimmed_reads/{sample}_R1.fastq.gz',
        R2 = 'trimmed_reads/{sample}_R2.fastq.gz',
        html_report = 'fastp_report/{sample}.html',
        json_report = 'fastp_report/{sample}.json'
    conda:
        "envs/fastp.yaml"
    threads: 4 #how many?
    shell: 
        "fastp -w {threads} \ "
        "-i {input.R1} -I {input.R2} \ "
        "-o {output.R1} -O {output.R2} \ "
        "-h {output.html_report} -j {output.json_report}"

rule hisat2:
    input:
        R1 = 'trimmed_reads/{sample}_R1.fastq.gz',
        R2 = 'trimmed_reads/{sample}_R2.fastq.gz',
        genome = config["genome_info"]["genome"] #genome index
    output:
        bam = 'aligned_reads/{sample}.bam',
        rep = 'hisat_report/{sample}_summary.txt'
    conda:
        "envs/hisat2.yaml"
    threads: 4 # how many?
    params: 
        min_i = config["hisat_params"]["min_i"], #min intronlen
        max_i = config["hisat_params"]["max_i"] #max intronlen
    shell:
        "hisat2 --min-intronlen {params.min_i} --max-intronlen {params.max_i} \ "
        "-p {threads} \ "
        "-x {input.genome} \ "
        "-1 {input.R1} -2 {input.R2} \ "
        "2> {output.rep} | samtools view -bS - \ "
        "> {output.bam}"

rule featureCounts_unionExon:
    input:
        bam = 'aligned_reads/{sample}.bam',
        anno = config["genome_info"]["gtf"]
        #count_by
        #summarize_by
    output:
        'featureCounts/{sample}.txt'
    conda:
        "envs/subread.yaml"
    threads: 4 # how many?
    shell:
        "featureCounts -T {threads} \ "
        "-a '{input.anno}' \ "
        "-t exon -g gene_id \ "
        "-o {output} \ "
        "{input.bam}"