##############################################################
# Only get chromosomes 1-5 of the TAR10 genome
#
# Author: Kobie Kirven
##############################################################

# imports 
from Bio import SeqIO

dir_path = "/run/media/kjk6173/kobie_ext/"
fn = open(dir_path + "TAIR10_genomic.fa", "w")

for rec in SeqIO.parse(dir_path + "TAIR10_chr_all.fas", "fasta"):
    if rec.id in ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]:
        fn.write(">" + rec.id + "\n" + str(rec.seq) + "\n")
fn.close()

