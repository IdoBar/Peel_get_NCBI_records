# Cleaning Reference Databases - 60+ species --------------------------------------------
# Cinthia Pietromonaco
# Test 1 - with species exported from ala within mbrc
# Rules: 1. complete sequence 2. herbarium voucher  3.ITS2  4. sequence length between 200-2500 5.Viridiplantae
#load packages
library(readxl)
library(stringr)
library(tidyr)
library(dplyr) 
library(readr)
library(purrr)

library(tidyverse) # this includes all of the packages above 
library(glue)

# Step 1: Search & Download Species List ----------------------------------
#upload list
#make sure spelling is ok and if there are synonyms for certain species.
taxa <- read_csv("records-2022-12-01.csv") %>% select(`Species Name`) %>% 
  mutate(spp_string=glue('"{`Species Name`}"[Organism] OR'))
taxa <- taxa[,2]
#species <- as.character(taxa$'Species Name')

#string 1 for ITS2 marker
searchdef <- '((internal transcribed spacer[All Fields]) OR (ITS[All Fields])) AND (plants[filter]) AND (200[SLEN] : 2500[SLEN]) AND (flowering plants[porgn]) AND'

#string for species (do you have the special double quotes in purpose or is t because it was copied-pasted from a web page or a word document?)
species.str <- '“SPECIES_NAME”[Organism] OR '

# lists greater than 60 species -------------------------------------------
add.taxa <- taxa %>%
  mutate(spp_string = str_replace(species.str,"SPECIES_NAME", `Species Name`))

#create empty vector to store output
genbank.search <- c()

for(i in seq(1, nrow(add.taxa), 30)) {
  tmp <- paste(add.taxa$spp_string[i:(i+29)], collapse = " ")
  # Remove trailing NA for groups with less than 60 species
  tmp <- gsub("(\\sNA){2,}$", "", tmp)
  # remove trailing " OR "
  tmp <- gsub(" OR $", "", tmp)
  tmp <- paste(tmp, sep ="", collapse = "")
  tmp <- paste('(', tmp, sep ="", collapse = "")
  tmp <- paste(tmp, ')', sep ="", collapse = "")
  tmp <- paste(searchdef, tmp)
  genbank.search <- c(genbank.search, tmp)
}
# 2: Download FASTA files from GenBank ---------------------------------------
library(rentrez)
# set NCBI API key - this allows yo to send more queries to Entrez without timing out - see more details at https://support.nlm.nih.gov/knowledgebase/article/KA-05318/en-us and https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us
set_entrez_key("33ec35e1973e4928beafd47bee86492be109")

#entrez_db_summary('nucleotide')	 #checking its the right database
#entrez_db_searchable("nucleotide") #checking what search terms we can use

#create loop to extract ids for each species sequence result
empty_data <- data.frame() #make empty data frame
prova_EF <- lapply(genbank.search, FUN = function(a) {
  r_search = entrez_search(db="nucleotide",
                           term=a,
                           retmax=5000)
  search_ids <- data.frame(r_search$ids) # Convert search results to a data frame
  empty_data <- rbind(empty_data, search_ids)
})

# You can do this much more efficiently with purrr:
prova_EF <- genbank.search %>% 
  map_dfr(.f = ~entrez_search(db="nucleotide",
                              term=.x,
                              retmax=5000))


# OLD VERSION -------------------------------------------------------------
#for (i in 1:nrow(genbank.search)){
#  r_search = entrez_search(db="nucleotide",
#                           term=genbank.search[i,],
#                           retmax=12000)
#  search_ids <- data.frame(r_search$ids) # Convert search results to a data frame
#  empty_data <- rbind(empty_data, search_ids)    # Append the rows of the search results to the empty data frame
#}

# -------------------------------------------------------------------------

#create loop to extract title information about each sequence
features <- data.frame()
prova_CP <- lapply(prova_EF[1:2], FUN = function(b) {
    apply(b, MARGIN = 1, FUN = function(c) {
    output = entrez_fetch(db="nucleotide", id=c, rettype="text")
    output = data.frame(ids = c, features = output)
    features = rbind(output, features)   
  })
})

data <- map_dfr(prova_CP, bind_rows)

#this will take awhile if you have a lot of ids. go take a coffee break :)
############################## stopped here!!

# Step 3: filtering ---------------------------------------------------------------
#nucleotide.summary <- entrez_summary(db="nucleotide",id=data$ids) #see variables to consider

with_voucher <- dplyr::filter(data, grepl('specimen_voucher', features))
#View(with_voucher)

full_seq <- dplyr::filter(with_voucher, !grepl('internal transcribed spacer 2, partial sequence', features))
full_seq <- dplyr::filter(full_seq, !grepl('internal transcribed spacer 2-like', features))

data2 <- data.frame(ids = character(), sequence_inf = character())

# Fetch the nucleotide sequence for each ID in full_seq$ids
for (n in full_seq$ids){
  sequence_inf = entrez_fetch(db="nucleotide", id=n, rettype="fasta")
  # Add the ID and sequence to the sequence_inf data frame
  sequence_inf <- data.frame(ids = n, sequence_inf = sequence_inf)
  # Add the row to the data2 data frame
  data2 <- rbind(data2, sequence_inf)
}

#View(data2)

genbank_sequences <- merge(data2, full_seq, by = "ids") %>% 
  drop_na()

# Step 4: separate and rename FASTA text accordingly ------------------------------
#accession
#View(genbank_sequences)
accession <- data.frame(accession = as.character(sapply(strsplit(genbank_sequences$sequence_inf, " "), "[", 1)))
accession$accession <- gsub('>','',accession$accession)

with_accession <- cbind(genbank_sequences, accession)

#View(with_accession)

#species name
species <- data.frame(species = as.character(word(with_accession$sequence_inf, 2,3, sep=" ")))
with_species <- cbind(with_accession, species)

#sequence
sequence <- data.frame(sequence = as.character(word(with_species$sequence_inf,-1)))
sequence <- sequence %>%
  mutate(sequence = str_sub(sequence, 9, -1))
with_sequence <- cbind(with_species, sequence)

#View(with_sequence)

#Step 5: Get taxa info -----------------------------------------------------------
species_list <- (unique(with_sequence$species)) #get unique species so it doesn't take as long
#create loop to extract each taxonomy
#need to have a good internet connection for this one.. if not it may need to be redone several times to retrieve all info without error
emptydata = list()

for (n in species_list) {
  taxa_inf = data.frame(myTAI::taxonomy(organism = n, #export information based on species
                                        db       = "ncbi",
                                        output   = "classification" ))
  taxa_inf$n <- n  # maybe you want to keep track of which iteration produced it?
  emptydata[[n]] <- taxa_inf
}


big_data <- data.table::rbindlist(emptydata)
#View(big_data)

big_data <- big_data[big_data$rank %in% c('kingdom', 'phylum', 'order', 'family', 'genus', 'species'),] #export taxa wanted
big_data <- select(big_data, -id)
big_data <- dplyr::mutate(big_data,#create link
                          link = dplyr::case_when(rank == "kingdom" ~"tax=k:",
                                                  rank == "phylum" ~"p:",
                                                  rank == "order" ~"o:",
                                                  rank == "family" ~"f:",
                                                  rank == "genus" ~"g:",
                                                  rank == "species" ~"s:"))
big_data$link_name <- paste(big_data$link, big_data$name, sep = "")

with_sequence <- with_sequence %>% 
  rename("n" = "species")

total <- merge(with_sequence,big_data, by="n")
#change the link_name value for species to the one from the dataset because myTAI may export synonyms instead
total$link_name <- ifelse(total$rank == "species", paste0("s:", total$n), total$link_name)

indiv <- split(total, total$ids) #make a nested list with each species

result_list <- list() #make empty list
#create loop to make link for each sp. in nested list
for (i in 1:length(indiv)) {
  result <- paste(indiv[[i]][["link_name"]], sep = "", collapse = ",")
  # result <- str_replace(result," ", "_")
  result_list <- c(result_list, result)
}

#View(result_list)
taxa_link <- unlist(result_list)
result_df <- taxa_link %>%
  as.data.frame(taxa_link, stringsAsFactors = FALSE) %>%
  mutate(n = str_remove(str_extract(taxa_link, "s:.*"), "s:")) %>%
  mutate(n = str_replace(n, "_", " ")) %>%
  rename_at(1,~"taxa_link") 

total <- merge(total,result_df, by="n")

total$taxa_link <- paste(total$accession, total$taxa_link, sep = ";")
total$taxa_link <- paste(total$taxa_link, ";", sep = "")
total$taxa_link <- str_replace(total$taxa_link," ", "_")

end <-  total[!duplicated(total$ids), ]
# Step 6 - Create FASTA file ----------------------------------------------
library(seqinr)
write.fasta(as.list(end$sequence),names = end$taxa_link, 
            file.out = "/Users/cinthiapietromonaco/Desktop/its2.somethingthatmakessense.date.tax.fasta")

