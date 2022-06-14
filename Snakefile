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
SAMPLE_NAMES = ["AA", "AB", "AC", "AD","AE", "AF", "AM", "AN", "AO", "AP", "AQ", "AR"]
INTERMEDIATE_DIRECTORIES = ["160411_7001126F_0117_AHKWWWBCXX_1","160415_7001126F_0121_AHLFCKBCXX_2",
"160503_7001126F_0122_BHLCYKBCXX_3", "160517_7001126F_0125_AHMWNJBCXX_4", "160527_7001126F_0126_BHMN3LBCXX_5"]
group_names = ["normal", "salt_stress"]
###################
###### Rules ######
###################

prefixes = expand(config["working_directory"] + "0_raw/trimmed/{dir}/fastq/", dir=INTERMEDIATE_DIRECTORIES)

rule all:
    input:
        expand(config["working_directory"] + "0_raw/trimmed/{wild}/fastq/{sample_name}_trimmed_fastqc.html", wild=INTERMEDIATE_DIRECTORIES, sample_name=SAMPLE_NAMES),
        expand(config["working_directory"] + "0_raw/variants/{sample}_Aligned.sortedByCoord.dup.sorted.reformat.vcf", sample=SAMPLE_NAMES),

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

# Merge fasta files: Merge the multiple sequencing runs for each biological replicate
rule merge_fasta:
    input:
        [fileName + "{sample}_trimmed.fastq.gz" for fileName in prefixes]
    output:
        config["working_directory"] + "0_raw/trimmed/merged/{sample}_trimmed.fastq.gz"
    shell:
        "cat {input} > {output}"

# Get reference transcriptome: Download the reference transcriptome 
rule get_reference_transcriptome:
    output:
        config["working_directory"] + "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
    shell:
        "cd " + config["working_directory"] + "&& " + "wget http://ftp.ensemblgenomes.org/pub/plants/release-53/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"

# Get reference GTF: Download the reference GTF
rule get_reference_gtf:
    output:
        config["working_directory"] + "Arabidopsis_thaliana.TAIR10.53.gtf"
    shell:
        "cd " + config["working_directory"] + "&& " + "wget http://ftp.ensemblgenomes.org/pub/plants/release-53/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.53.gtf.gz && gunzip Arabidopsis_thaliana.TAIR10.53.gtf.gz"

# # Get chromosomes: Get the non-chloroplast and non-mitochondrial chromosomes from the reference transcriptome
# rule get_chromosomes:
#     input:
#         config["working_directory"] + "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
#     output: 
#         config["working_directory"] + "TAIR10_genomic.fa"
#     shell:
#         "python3 format_genome.py -i " + config["working_directory"] + "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz -o " + config["working_directory"] + "TAIR10_genomic.fa"


# Build STAR reference genome
rule build_STAR_reference:
    input:
        gtf = config["working_directory"] + "Arabidopsis_thaliana.TAIR10.53.gtf",
        fas = config["working_directory"] + "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
    output:
       config["working_directory"] + "STAR_genome/sjdbList.fromGTF.out.tab"
    shell:
        "STAR --runMode genomeGenerate --genomeDir " + config["working_directory"] + "STAR_genome --genomeFastaFiles {input.fas} --sjdbGTFfile {input.gtf} --sjdbOverhang 100 --runThreadN 4 --genomeSAindexNbases 12"

# Create sequence dictonary
rule create_sequence_dictionary:
    input:
        config["working_directory"] + "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
    output:
        config["working_directory"] + "Arabidopsis_thaliana.TAIR10.dna.toplevel.dict"
    shell:
        "java -jar /home/kjk6173/bin/picard.jar CreateSequenceDictionary R={input} O={output}"

# Create a fasta index file
rule create_fasta_index:
    input:
        config["working_directory"] + "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
    output:
        config["working_directory"] + "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai"
    shell:
        "samtools faidx {input}"

# Align reads to the reference transcriptome
rule align_reads:
    input:
        reads = config["working_directory"] + "0_raw/trimmed/merged/{sample}_trimmed.fastq.gz",
        genome = config["working_directory"] + "STAR_genome/sjdbList.fromGTF.out.tab",

    params:
        prefix = config["working_directory"] + "0_raw/alignment/star/{sample}_"
    output:
        out = config["working_directory"] + "0_raw/alignment/star/{sample}_Aligned.sortedByCoord.out.bam",

    shell:
        "STAR --genomeDir " + config["working_directory"] + "STAR_genome --readFilesIn {input.reads} --outFileNamePrefix {params.prefix} " + config["working_directory"] + "0_raw/alignment/star/{wildcards.sample}_ --outSAMtype BAM SortedByCoordinate --readFilesCommand 'gunzip -c' --runThreadN 8 --sjdbOverhang 100 --twopassMode Basic --outTmpDir /home/kjk6173/Desktop/rna_mutations/temp_{wildcards.sample}"

# Mark duplicates in the alignment file
rule mark_duplicates:
    input:
        config["working_directory"] + "0_raw/alignment/star/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        bam = config["working_directory"] + "0_raw/alignment/star/{sample}_Aligned.sortedByCoord.dup.bam",
        metrics = config["working_directory"] + "0_raw/alignment/star/{sample}_Aligned.sortedByCoord.dup.metrics"
    shell:
        "java -jar /home/kjk6173/bin/picard.jar MarkDuplicates --INPUT {input} --OUTPUT {output.bam} --METRICS_FILE {output.metrics} --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT"


# Sort the alignment file
rule sort_bam:
    input:
        config["working_directory"] + "0_raw/alignment/star/{sample}_Aligned.sortedByCoord.dup.bam"
    output:
        config["working_directory"] + "0_raw/alignment/star/{sample}_Aligned.sortedByCoord.dup.sorted.bam"
    shell:
        "java -jar /home/kjk6173/bin/picard.jar SortSam --INPUT {input} --OUTPUT {output} --SORT_ORDER coordinate --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT"

# Split N cigar reads
rule split_N_cigar:
    input:
        bam = config["working_directory"] + "0_raw/alignment/star/{sample}_Aligned.sortedByCoord.dup.sorted.bam",
        ref_dict = config["working_directory"] + "Arabidopsis_thaliana.TAIR10.dna.toplevel.dict",
        index = config["working_directory"] + "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai"
    output:
        config["working_directory"] + "0_raw/alignment/star/{sample}_Aligned.sortedByCoord.dup.sorted.split.bam"
    shell:
        "java -Xmx20g -jar /home/kjk6173/bin/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar SplitNCigarReads -I {input.bam} -O {output} -R /run/media/kjk6173/kobie_ext/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"


# Recalibrate base quality scores
rule recalibrate_base_quality:
    input:
        config["working_directory"] + "0_raw/alignment/star/{sample}_Aligned.sortedByCoord.dup.sorted.split.bam"
    output:
        config["working_directory"] + "0_raw/alignment/star/{sample}_Aligned.sortedByCoord.dup.sorted.table"
    shell:
        "java -jar /home/kjk6173/bin/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar BaseRecalibrator -I {input} -O {output} -R Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"


# Reformat the algnment based on the recalibration table
rule reformat_alignment:
    input:
        bam = config["working_directory"] + "0_raw/alignment/star/{sample}_Aligned.sortedByCoord.dup.sorted.split.bam",
        table = config["working_directory"] + "0_raw/alignment/star/{sample}_Aligned.sortedByCoord.dup.sorted.table"
    output:
        config["working_directory"] + "0_raw/alignment/star/{sample}_Aligned.sortedByCoord.dup.sorted.reformat.bam"
    shell:
        "java -jar /home/kjk6173/bin/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar ApplyBQSR -R /run/media/kjk6173/kobie_ext/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa --bqsr-recal-file {input.table} -I {input.bam} -O {output}"


# Call variants
rule call_variants:
    input:
        config["working_directory"] + "0_raw/alignment/star/{sample}_Aligned.sortedByCoord.dup.sorted.reformat.bam"
    output:
        config["working_directory"] + "0_raw/variants/{sample}_Aligned.sortedByCoord.dup.sorted.reformat.vcf"
    shell:
        "java -jar /home/kjk6173/bin/gatk-package-4.2.6.1-local.jar -T HaplotypeCaller -R /run/media/kjk6173/kobie_ext/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -I {input} -dontUseSoftClippedBases -stand_call_conf 20 -o {output}"