
############### Hasan Can Demirci #################

############## GGMSA by taking data from NCBI #####


# necessary libraries

library(ggmsa)
library(compbio4all)
library(rentrez)
library(Biostrings)
library(msa)


# creating data with accession numbers

TP53_table <- c("NP_000537.3", "Human", "TP53",
                 "NP_035770.2 ","Mus Musculus","TP53",
                 "NP_001258749.1","ZebraFish","TP53")




# convert the vector to matrix using matrix()

TP53_table_matrix <- matrix(TP53_table, byrow = T, nrow = 3)

# convert the matrix to a dataframe using data.frame()
# stringsAsFactors: used to determine the data change into vector or factor

TP53_table <- as.data.frame(TP53_table_matrix, stringsAsFactors = F)

# name columns of dataframe using names() function

colnames(TP53_table) <- c("ncbi.protein.accession", "species", "gene.name")

# get fasta sequences

TP53_list <- entrez_fetch_list(db = "protein",
                                id = TP53_table$ncbi.protein.accession,
                                rettype = "fasta")

# Go through list and remove FASTA header from each sequence
# fasta_cleaner : provide reading fasta format file as DNA, RNA, or Protein depends on your choice

TP53_vector <- sapply(TP53_list, fasta_cleaner, parse = F)
TP53_vector
# name the vector

names(TP53_vector) <- TP53_table$species

# create AAStringSet for protein alignment

TP53_vector_ss <- AAStringSet(TP53_vector)

# perform MSA (visualization)

alignment <- msa(TP53_vector_ss, "ClustalW")

# change class of alignment to AAMultipleAlignment

class(alignment) <- "AAMultipleAlignment"

# plot the alignment using ggmsa

ggmsa::ggmsa(alignment, start = 1, end = 40, char_width = 0.5, seq_name =T) +
  geom_seqlogo() +
  geom_msaBar()

