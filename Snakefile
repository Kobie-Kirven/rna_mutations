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
        expand(config["working_directory"] + "0_raw/trimmed/{wild}/fastq/{sample_name}_trimmed.fastq.gz", wild=INTERMEDIATE_DIRECTORIES, sample_name=SAMPLE_NAMES)

rule trim_reads:
    #Trim single-end reads with cutadapt
    input:
        config["working_directory"] + "0_raw/{dir}/fastq/Sample_{sample}/{sample}_R1.fastq.gz"
    output:
        config["working_directory"] + "0_raw/trimmed/{dir}/fastq/{sample}_trimmed.fastq.gz"
    shell:
        "cutadapt -a AGATCGGAAGAG -q 30,30 -m 20 -o {output} {input}"

# rule run_fastqc:
#     "Run fastqc on trimmed reads"
#     input:

