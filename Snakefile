#load info from config file
#config file contains sample name info, genome file info, and hisat parameters
#consider adding featureCounts summary factors as config parameters
configfile: "config.yaml"

#use minimal debian containerized environment with conda
#useful for OS standardization
container: "docker://continuumio/miniconda3:4.5.11"

#this rule looks for all the final files
#drives the back propagation of all intermediate rules
#elimintates the need for specifying inputs or outputs on command line
rule all:
    input:
        expand('featureCounts/{sample}.txt', sample=config["samples"])

#fastp in paired-end mode
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
        "fastp -w {threads} "
        "-i {input.R1} -I {input.R2} "
        "-o {output.R1} -O {output.R2} "
        "-h {output.html_report} -j {output.json_report}"

#hisat2 in paired end mode
rule hisat2:
    input:
        R1 = 'trimmed_reads/{sample}_R1.fastq.gz',
        R2 = 'trimmed_reads/{sample}_R2.fastq.gz'
    output:
        bam = 'aligned_reads/{sample}.bam',
        report = "hisat2_report/{sample}_log.txt"
    conda:
        "envs/hisat2.yaml"
    threads: 10 # how many?
    params: 
        genome = config["genome_info"]["genome"], #genome index
        min_i = config["hisat_params"]["min_i"], #min intronlen
        max_i = config["hisat_params"]["max_i"] #max intronlen
    shell:
        "PERL5LIB=/quobyte/Leiboff_Lab/Lgn_nod_RNAseq/.snakemake/conda/f58acdc0/lib/5.30.3/Perl; "
        "hisat2 --min-intronlen {params.min_i} --max-intronlen {params.max_i} "
        "-p {threads} "
        "-x {params.genome} "
        "-1 {input.R1} -2 {input.R2} "
        "--summary-file {output.report} "
        "| samtools view -Sb -F 4 -o {output.bam}"

#Union-Exon RNA alignment counting
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
        "featureCounts -T {threads} "
        "-a '{input.anno}' "
        "-t exon -g gene_id "
        "-o {output} "
        "{input.bam}"