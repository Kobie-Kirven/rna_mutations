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
        expand(config["working_directory"] + "0_raw/alignment/{wild}/{sample_name}.sam", wild=INTERMEDIATE_DIRECTORIES, sample_name=SAMPLE_NAMES)

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
        "Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"
    shell:
        "wget http://ftp.ensemblgenomes.org/pub/plants/release-53/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"


# Build ref transcriptome: Build the botwite2 index for the reference transcriptome 
rule build_ref_transcriptome:
    input:
        "Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"
    output:
        "Arabidopsis_thaliana.TAIR10.cdna.all.1.bt2"
    shell:
        "bowtie2-build Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz Arabidopsis_thaliana.TAIR10.cdna.all"

# Align reads: Align reads with bowtie2 
rule align_reads:
    input: 
        config["working_directory"] + "0_raw/trimmed/{dir}/fastq/{sample}_trimmed.fastq.gz"
    output:
        config["working_directory"] + "0_raw/alignment/{dir}/{sample}.sam"
    shell:
        "bowtie2 --very-sensitive-local -x Arabidopsis_thaliana.TAIR10.cdna.all -U {input} -S {output}"