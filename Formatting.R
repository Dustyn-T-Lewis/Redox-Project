# Loading the data "stringAsFactors" function ensures that text (character) data are not converted
#     into categorical factors automatically which could make string manipulations harder
getwd()
data <- read.csv("Master Proteomics.csv", stringsAsFactors = FALSE)
colnames(data)

# Load Required Libraries
library(dplyr)
library(stringr)

# Install and Load `rentrez` for UniProt Queries
if (!requireNamespace("rentrez", quietly = TRUE)) {
  install.packages("rentrez")
}
library(rentrez)

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

# Extracting UniProt Accession Numbers from Gene_Proteins
accession_extract <- function(Gene_Proteins) {
  accession <- str_extract_all(Gene_Proteins, "[OPQ][0-9][A-Z0-9]{3}[0-9]")  
  if (length(accession[[1]]) > 0) {  
    return(accession[[1]][1])  # Return the first UniProt Accession match
  } else {
    return(NA)  
  }
}

# Inserting gene names, descriptions, and accessions into table
proteomics_data <- data %>%
  mutate(Gene_Names = sapply(Gene_Proteins, gene_name_extract),
         Gene_Descriptions = sapply(Gene_Proteins, gene_description_extract),
         Accession = sapply(Gene_Proteins, accession_extract))

# There are two issues at this point: 
# 1) Some rows have gene names but no accession number,
# 2) Some rows have accession numbers but no gene names
# UniProt can use gene names to retrieve accessions and accessions to retrieve gene names

# Function to Query UniProt for Accession Numbers by Gene Name
get_accession_by_gene <- function(Gene_Names) {
  if (is.na(Gene_Names) || Gene_Names == "" || Gene_Names == "GN=") {
    return(NA)
  }
  
  query <- paste(Gene_Names, "[GENE] AND Homo sapiens[ORGN]")
  result <- entrez_search(db = "protein", term = query)
  
  if (!is.null(result) && length(result$ids) > 0) {
    return(result$ids[1])  # Return the first UniProt Accession match
  } else {
    return(NA)
  }
}

# Function to Query UniProt for Gene Name by Accession Number
get_gene_by_accession <- function(accession) {
  if (is.na(accession) || accession == "") {
    return(NA)
  }
  
  result <- entrez_summary(db = "protein", id = accession)
  gene_name <- result$organism
  
  if (!is.null(gene_name)) {
    return(gene_name)  # Return the gene name
  } else {
    return(NA)
  }
}

# Mutating Proteomics Data Using Batches to Avoid API Limits
batch_size <- 100  # Adjust based on your dataset size
num_batches <- ceiling(nrow(proteomics_data) / batch_size)

for (i in 1:num_batches) {
  start <- (i - 1) * batch_size + 1
  end <- min(i * batch_size, nrow(proteomics_data))
  
  proteomics_data[start:end, "Accession"] <- ifelse(
    is.na(proteomics_data[start:end, "Accession"]) | proteomics_data[start:end, "Accession"] == "",
    sapply(proteomics_data[start:end, "Gene_Names"], get_accession_by_gene),
    proteomics_data[start:end, "Accession"]
  )
  
  proteomics_data[start:end, "Gene_Names"] <- ifelse(
    is.na(proteomics_data[start:end, "Gene_Names"]) | proteomics_data[start:end, "Gene_Names"] == "GN=",
    sapply(proteomics_data[start:end, "Accession"], get_gene_by_accession),
    proteomics_data[start:end, "Gene_Names"]
  )
}

# Save Processed Data to CSV
save_proteomics_data <- function(data, filename = "Processed_Proteomics_Data.csv") {
  write.csv(data, file = filename, row.names = FALSE)
  message(paste("✅ Data successfully saved to:", filename))
}

# Call the function to save your data
save_proteomics_data(proteomics_data)

View(proteomics_data)
