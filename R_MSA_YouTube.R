
############### Hasan Can Demirci #################

############## GGMSA by taking data from NCBI #####


# necessary libraries

library(ggmsa)
library(compbio4all)
library(rentrez)
library(Biostrings)
library(msa)


# creating data with accession numbers

CD150_table <- c("NP_001108225.1", "Human", "CD150",
                 "NP_001074356.2","Gallus gallsus","CD150",
                 "NP_031958.2","Mus Musculus","CD150")




# convert the vector to matrix using matrix()

CD150_table_matrix <- matrix(CD150_table, byrow = T, nrow = 3)

# convert the matrix to a dataframe using data.frame()
# stringsAsFactors: used to determine the data change into vector or factor

CD150_table <- as.data.frame(CD150_table_matrix, stringsAsFactors = F)

# name columns of dataframe using names() function

colnames(CD150_table) <- c("ncbi.protein.accession", "species", "gene.name")

# get fasta sequences

CD150_list <- entrez_fetch_list(db = "protein",
                                id = CD150_table$ncbi.protein.accession,
                                rettype = "fasta")

# Go through list and remove FASTA header from each sequence
# fasta_cleaner : provide reading fasta format file as DNA, RNA, or Protein depends on your choice

CD150_vector <- sapply(CD150_list, fasta_cleaner, parse = F)

# name the vector

names(CD150_vector) <- CD150_table$species

# create AAStringSet for protein alignment

CD150_vector_ss <- AAStringSet(CD150_vector)

# perform MSA (visualization)

alignment <- msa(CD150_vector_ss, "ClustalW")

# change class of alignment to AAMultipleAlignment

class(alignment) <- "AAMultipleAlignment"

# plot the alignment using ggmsa

ggmsa::ggmsa(alignment, start = 1, end = 40, char_width = 0.5, seq_name =T) +
  geom_seqlogo() +
  geom_msaBar()

