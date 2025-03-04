#     Loading the data "stringAsFactors" function ensures that text (character) data are not converted
#     into categorical factors automatically which could make string manipulations harder
getwd()
data <- read.csv("Master Proteomics.csv", stringsAsFactors = FALSE)
colnames(data)
library(dplyr)
library(stringr)


#     Installing biomaRt package to fetch accession numbers from Ensemble database
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

library(biomaRt)
library(biomartr)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://useast.ensembl.org")


# Extracting Gene names from Gene_Proteins
gene_name_extract <- function(Gene_Proteins) {
  if (grepl(";", Gene_Proteins)) {
    return(str_split(Gene_Proteins, ";")[[1]][1])
  } else {
    return(str_split(Gene_Proteins, "_")[[1]][1])
  }
}


# Extracting Gene descriptions from Gene_Proteins
gene_description_extract <- function(Gene_Proteins) {
  if (grepl("GN=", Gene_Proteins)) { 
    desc <- str_extract(Gene_Proteins, "GN=.*")  
    desc <- str_remove(desc, "GN=")  
    desc <- str_remove_all(desc, "GN=;*")  
    return(str_trim(desc))  
  } else if (grepl("_", Gene_Proteins)) {  
    return(str_split(Gene_Proteins, "_")[[1]][2])  
  } else {
    return(NA)  
  }
}


# Extracting Accession Numbers from Gene_Proteins
accession_extract <- function(Gene_Proteins) {
  accession <- str_extract_all(Gene_Proteins, "[OPQ][0-9][A-Z0-9]{3}[0-9]")  
  if (length(accession[[1]]) > 0) {  
    return(accession[[1]][1])  # Return the first match
  } else {
    return(NA)  
  }
}


# Inserting gene names, descriptions, and accessions into table
proteomics_data <- data %>%
  mutate(Gene_Names = sapply(Gene_Proteins, gene_name_extract),
         Gene_Descriptions = sapply(Gene_Proteins, gene_description_extract),
         Accession = sapply(Gene_Proteins, accession_extract))


#There are two issues at this point: 
# 1) Some rows have gene names but no accession number,
# 2) Some rows have accession numbers but no gene names
#BiomaRt can use gene names used to retrieve accessions and accessions can be used to retrieve gene names

#1 getting accession
get_accession_by_gene <- function(Gene_Names) {
  # 1. If the gene name is missing, empty, or just "GN=", return NA
  if (is.na(Gene_Names) || Gene_Names == "" || Gene_Names == "GN=") {
    return(NA)
  }
  
  # 2. Query Ensembl using biomaRt to find the corresponding accession number
  result <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),  # What to return; 'getBM' searches for accession using gene name
                  filters = "hgnc_symbol",  # Search by gene name
                  values = Gene_Names,  # The gene name to look up
                  mart = ensembl)  # Ensembl database connection
  
  # 3. If a match is found, return the first accession number
  if (nrow(result) > 0) {
    return(result$ensembl_gene_id[1])  # Take the first result
  } else {
    return(NA)  # If no match, return NA
  }
}

#2 getting gene
get_gene_by_accession <- function(accession) {
  if (is.na(accession) || accession == "") {
    return(NA)
  }
  
  result <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = accession,
                  mart = ensembl)
  
  if (nrow(result) > 0) {
    return(result$hgnc_symbol[1])  # Return the first match
  } else {
    return(NA)
  }
}

#mutate proteomics to add data
proteomics_data <- proteomics_data %>%
  mutate(
    Accession = ifelse(is.na(Accession) | Accession == "", 
                       sapply(Gene_Names, get_accession_by_gene), 
                       Accession),
    Gene_Names = ifelse(is.na(Gene_Names) | Gene_Names == "GN=", 
                        sapply(Accession, get_gene_by_accession), 
                        Gene_Names)
  )
