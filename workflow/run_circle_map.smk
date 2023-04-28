configfile: "workflow/config.yml"

def get_samples(input):
    samples = {}
    with open(input,'r') as f:
        for line in f:
            sample,*file = line.strip("\n").split("\t")
            samples[sample] = file
    return samples

samples = get_samples(config['samples'])

def get_fastq_files(wildcards):
    return samples[wildcards.sample]

rule all:
    input:
        lambda wildcards: expand("report/{sample}_summary.html", sample=samples.keys()),
        
rule fastqc:
    input:
        get_fastq_files,
    output:
        directory("FastQC/{sample}")
    group:
        "precessing"
    threads:
        config['threads']
    resources:
        mem_gb=config['Normal_mem']
    singularity:
        "workflow/softwares/fastqc.sif"
    shell:
        "mkdir -p {output} &&"
        "fastqc -t {threads} -o {output} {input}"

rule fastp:
    input:
        get_fastq_files,
        "FastQC/{sample}"
    output:
        temp("Fastp/{sample}.clean_1.fq.gz"),
        temp("Fastp/{sample}.clean_2.fq.gz"),
        protected("Fastp/{sample}.clean.html"),
        protected("Fastp/{sample}.clean.json")
    group:
        "precessing"
    threads:
        config['threads']
    resources:
        mem_gb=config['Normal_mem']
    singularity:
        "workflow/softwares/fastp.sif"
    shell:
        "fastp -w {threads} -i {input[0]} -o {output[0]} -I {input[1]} -O {output[1]} -h {output[2]} -j {output[3]}"

rule bwa:
    input:
        "Fastp/{sample}.clean_1.fq.gz",
        "Fastp/{sample}.clean_2.fq.gz",
    output:
        temp("alignData/{sample}_circle.sam")
    group:
        "precessing"
    threads:
        config['threads']
    resources:
        mem_gb=config['Normal_mem']
    singularity:
        "workflow/softwares/bwa.sif"
    shell:
        "bwa mem -q -t {threads} {config[reference]} {input[0]} {input[1]} > {output}"

rule coordinate_sort_index:
    input:
        "alignData/{sample}_circle.sam"
    output:
        protected("preData/sorted_{sample}_circle.bam")
    group:
        "precessing"
    threads:
        config['threads']
    resources:
        mem_gb=config['Normal_mem']
    singularity:
        "workflow/softwares/samtools.sif"
    shell:
        "samtools sort -@ {threads} -o {output} {input} &&"
        "samtools index -@ {threads} {output}"

rule name_sort:
    input:
        "preData/sorted_{sample}_circle.bam"
    output:
        temp("preData/qname_{sample}_circle.bam")
    group:
        "precessing"
    threads:
        config['threads']
    resources:
        mem_gb=config['Normal_mem']
    singularity:
        "workflow/softwares/samtools.sif"
    shell:
        "samtools sort -n -@ {threads} -o {output} {input}"

rule extract_circle_reads:
    input:
        "preData/qname_{sample}_circle.bam"
    output:
        temp("preData/{sample}_circular_read_candidates.bam")
    group:
        "precessing"
    threads:
        config['threads']
    resources:
        mem_gb=config['Normal_mem']
    singularity:
        "workflow/softwares/CircleMapPlus.sif"
    shell:
        "circle_map++ ReadExtractor -t {threads} -i {input} -o {output}"

rule sort_circle_reads_index:
    input:
        "preData/{sample}_circular_read_candidates.bam"
    output:
        temp("preData/{sample}_circular_read_candidates.sort.bam")
    group:
        "precessing"
    threads:
        config['threads']
    resources:
        mem_gb=config['Normal_mem']
    singularity:
        "workflow/softwares/samtools.sif"
    shell:
        "samtools sort -@ {threads} -o {output} {input} &&"
        "samtools index -@ {threads} {output}"

rule circle_map_realign:
    input:
        "preData/{sample}_circular_read_candidates.sort.bam",
        "preData/qname_{sample}_circle.bam",
        "preData/sorted_{sample}_circle.bam",
    output:
        protected("FinallyData/{sample}_circle_site.bed")
    group:
        "fin"
    singularity:
        "workflow/softwares/CircleMapPlus.sif"
    threads:
        config['Big_threads']
    resources:
        mem_gb=config['Big_mem']
    shell:
        "circle_map++ Realign -t {threads} -i {input[0]} -qbam {input[1]} -sbam {input[2]} -fasta {config[reference]} -o {output}"

rule total_mappings:
    input:
        "preData/sorted_{sample}_circle.bam"
    output:
        temp("mappingState/{sample}.mapping.total.tsv")
    group:
        "downstream"
    threads:
        1
    resources:
        mem_gb=1
    script:
        "scripts/get_mapping.py"

rule cal_ecc_num:
    input:
        "FinallyData/{sample}_circle_site.bed"
    output:
        temp("BaseState/{sample}.eccNum.txt")
    group:
        "downstream"
    threads:
        1
    resources:
        mem_gb=1
    script:
        "scripts/cal_ecc_num.py"

rule calGC:
    input:
        ecc="FinallyData/{sample}_circle_site.bed"
    output:
        temp("BaseState/{sample}.gcContents.txt")
    group:
        "downstream"
    threads:
        1
    resources:
        mem_gb=1
    singularity:
        "workflow/softwares/Julia.sif"
    script:
        "scripts/calGC.jl"

rule calLength:
    input:
        "FinallyData/{sample}_circle_site.bed"
    output:
        temp("BaseState/{sample}.length.txt")
    group:
        "downstream"
    threads:
        1
    resources:
        mem_gb=1
    script:
        "scripts/calLength.py"

rule makeReport:
    input:
        mpTotal = "mappingState/{sample}.mapping.total.tsv",
        eccN = "BaseState/{sample}.eccNum.txt",
        gc = "BaseState/{sample}.gcContents.txt",
        lgth = "BaseState/{sample}.length.txt",
        ecc = "FinallyData/{sample}_circle_site.bed",
        fp = "Fastp/{sample}.clean.json"
    output:
        protected("report/{sample}_summary.html")
    group:
        "downstream"
    threads:
        1
    resources:
        mem_gb=1
    script:
        "scripts/makeReport.py"
