library("phyloseq")

#set working directory
setwd("/home/thibault/CROQ/METABARCODING/data_for_phyloseq")

#import files
otufile <- "filtered_biom-w-taxo.tsv"
mapfile <- "sample-metadata.csv"
treefile <- "tree.nwk"
rs_file <- "dna-sequences.fasta"
qiimedata <- import_qiime(otufile, mapfile, treefile, rs_file)

