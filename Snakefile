#################################################################################################################
# Look for mutations in Arabidopsis thaliana RNA to determine whether stressful conditions induce higher
# mutation rates in the RNA
#
# Author - Kobie Kirven
#################################################################################################################

###################
##### Inputs ######
###################
configfile: "config.yaml"
SAMPLE_NAMES = ["AA", "AB", "AC", "AD","AE", "AF"]
INTERMEDIATE_DIRECTORIES = ["160411_7001126F_0117_AHKWWWBCXX_1","160415_7001126F_0121_AHLFCKBCXX_2",
"160503_7001126F_0122_BHLCYKBCXX_3", "160517_7001126F_0125_AHMWNJBCXX_4", "160527_7001126F_0126_BHMN3LBCXX_5"]

###################
###### Rules ######
###################

rule all:
    input:
        expand(config["working_directory"] + "0_raw/trimmed/{wild}/fastq/{sample_name}_trimmed_fastqc.html", wild=INTERMEDIATE_DIRECTORIES, sample_name=SAMPLE_NAMES),
        #expand(config["working_directory"] + "0_raw/variants/{wild}/{sample_name}_variants.vcf", wild=INTERMEDIATE_DIRECTORIES, sample_name=SAMPLE_NAMES)
        expand(config["working_directory"] + "0_raw/alignment/{wild}/{sample_name}_sorted.bam", wild=INTERMEDIATE_DIRECTORIES, sample_name=SAMPLE_NAMES)

# Trim reads: Trim single-end reads with cutadapt
rule trim_reads:
    input:
        config["working_directory"] + "0_raw/{dir}/fastq/Sample_{sample}/{sample}_R1.fastq.gz"
    output:
        config["working_directory"] + "0_raw/trimmed/{dir}/fastq/{sample}_trimmed.fastq.gz"
    shell:
        "cutadapt -a AGATCGGAAGAG -q 30,30 -m 20 -o {output} {input}"


# Run FASTQC: get the FASTQC quality reports for the sequences 
rule run_fastqc:
    input:
        config["working_directory"] + "0_raw/trimmed/{dir}/fastq/{sample}_trimmed.fastq.gz"
    output:
        config["working_directory"] + "0_raw/trimmed/{dir}/fastq/{sample}_trimmed_fastqc.html"
    shell: 
        "fastqc {input}"


# Get reference transcriptome: Download the reference transcriptome 
rule get_reference_transcriptome:
    output:
        config["working_directory"] + "TAIR10_chr_all.fas.gz"
    shell:
        "cd " + config["working_directory"] + "&& " + "wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz"


# Build ref transcriptome: Build the botwite2 index for the reference transcriptome 
rule build_ref_transcriptome:
    input:
        config["working_directory"] + "TAIR10_chr_all.fas.gz"
    output:
        config["working_directory"] + "TAIR10_chr_all.fas.gz.1.bt2"
    shell:
        "bowtie2-build {input} {input}"

# Align reads: Align reads with bowtie2 
rule align_reads:
    input:
        bowtie = config["working_directory"] + "TAIR10_chr_all.fas.gz.1.bt2",
        reads = config["working_directory"] + "0_raw/trimmed/{dir}/fastq/{sample}_trimmed.fastq.gz",
        genome = config["working_directory"] + "TAIR10_chr_all.fas.gz"
    output:
        config["working_directory"] + "0_raw/alignment/{dir}/{sample}.bam"
    shell:
        "bowtie2 --very-sensitive-local -x {input.genome} -U {input.reads} | samtools view -bS - > {output}"

        
# Sort and Index: sort and index the alignments to the reference genome
rule sort_bam:
    input:
        config["working_directory"] + "0_raw/alignment/{dir}/{sample}.bam"
    output:
        config["working_directory"] + "0_raw/alignment/{dir}/{sample}_sorted.bam"
    shell:
        "samtools sort {input} > {output}"
        
# Index BAM: index BAM files
rule index_bam:
    input: 
        config["working_directory"] + "0_raw/alignment/{dir}/{sample}_sorted.bam"
    output:
        config["working_directory"] + "0_raw/alignment/{dir}/{sample}_sorted.bam.bai"
    shell:
        "samtools index {input}"

# Unzip reference: Unzip the reference genome
rule unzip_reference:
    input:
        config["working_directory"] + "{file}.fas.gz"
    output:
        config["working_directory"] + "{file}.fas"
    shell:
        "gunzip -c {input} > {output}"
    
    
# # Determine genotype likelyhoods: 
# rule determine_genotype_likelyhood:
#     input:
#         in_bam = config["working_directory"] + "0_raw/alignment/{dir}/{sample}_sorted.bam",
#         index = config["working_directory"] + "0_raw/alignment/{dir}/{sample}_sorted.bam.bai",
#         genome = config["working_directory"] + "TAIR10_chr_all.fas"
#     output:
#         config["working_directory"] + "0_raw/variants/{dir}/{sample}_genotypes.vcf"
#     shell:
#         "bcftools mpileup -Ov -f {input.genome} {input.in_bam} > {output}"


# # Call variants: Call variants in the data using bcftools
# rule call_variants:
#     input:
#         genotype = config["working_directory"] + "0_raw/variants/{dir}/{sample}_genotypes.vcf",
#         genome = config["working_directory"] + "TAIR10_chr_all.fas"
#     output:
#         config["working_directory"] + "0_raw/variants/{dir}/{sample}_variants.vcf"
#     shell:
#         "bcftools call --ploidy 2 -vm -Ov {input.genotype} > {output}"