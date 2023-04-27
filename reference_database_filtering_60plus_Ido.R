# Cleaning Reference Databases - 60+ species --------------------------------------------
# Cinthia Pietromonaco & Ido Bar
# Test 1 - with species exported from ala within mbrc
# Rules: 1. complete sequence 2. herbarium voucher  3.ITS2  4. sequence length between 200-2500 5.Viridiplantae
#install.packages("pak") # best package ever to install other packages, see https://pak.r-lib.org/
#library(pak)
# install packages for parallel processing, see https://furrr.futureverse.org/articles/progress.html
#pkg_install(c("furrr", "progressr")) 
# install utility packages that make life easy
# pacman provides the easiest way to load packages in one line, see https://trinker.github.io/pacman/vignettes/Introduction_to_pacman.html
#pkg_install(c("pacman", "janitor")) 
# install bioinformatics-related packages
#pkg_install(c("rentrez", "seqinr", "myTAI", "taxize", "taxonomizr"))

# load custom functions that Ido wrote 
devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "Util.R")

#load packages
pacman::p_load(tidyverse, rentrez, janitor, glue, furrr, progressr, seqinr, taxonomizr, tictoc)

# set NCBI API key - this allows yo to send more queries to Entrez without timing out - see more details at https://support.nlm.nih.gov/knowledgebase/article/KA-05318/en-us and https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us
set_entrez_key("33ec35e1973e4928beafd47bee86492be109")

# library(readxl)
# library(stringr)
# library(tidyr)
# library(dplyr) 
# library(readr)
# library(purrr)
# 
# library(tidyverse) # this includes all of the packages above 
# library(glue)
# 
# library(furrr)
# library(progressr)

# Step 1: Search & Download Species List ----------------------------------
#string 1 for ITS2 marker
searchdef <- '((internal transcribed spacer[All Fields]) OR (ITS[All Fields])) AND (plants[filter]) AND (200[SLEN] : 2500[SLEN]) AND (flowering plants[porgn]) AND'

species_name <- "SPECIES_NAME"
# species.str <- '"SPECIES_NAME"[Organism] OR '



#upload list
#make sure spelling is ok and if there are synonyms for certain species.
taxa <- read_csv("input_data/records-2023-03-16.csv") %>% select(`Species Name`) %>% 
  mutate(genebank.search=glue('{searchdef} "{`Species Name`}"[Organism]'))
# taxa <- taxa[,2]
# create test dataset
test_set <- head(taxa, n = 50)

test_search_string <- paste0(searchdef, " (",
                            paste(glue('"{test_set$`Species Name`}"[Organism]'), collapse = " OR "), ")")

# try to retrieve all search results for all species at once!!
res <- entrez_search(db="nuccore",
                     term=test_search_string,
                     use_history = TRUE)
res$web_history

#chunk searches query strings into groups of 50 species
char_limit <- 2048
species_chunks <- split(taxa$`Species Name`, ceiling(seq_along(taxa$`Species Name`)/50))
chunked_search_strings <- lapply(species_chunks, function(chunk) {
  paste0(searchdef, " (",
         paste(glue('"{taxa$`Species Name`}"[Organism]'), collapse = " OR "), ")")
})

chunked_search_strings <- data.frame(chunked_search_strings)

# Step 2: Download FASTA files from GenBank ---------------------------------------
#entrez_db_summary('nucleotide')	 #checking its the right database
#entrez_db_searchable("nucleotide") #checking what search terms we can use

  
# create a function to extract the title, accession and sequence from an Entrez id

get_NCBI_info <- function(res_id){
  # res_id = 1441206483
  # fetch this entry as FASTA
  res_fasta <- entrez_fetch(db='nuccore', id = res_id, rettype = c("fasta"))
  
  # res_features <- entrez_fetch(db="nuccore", id=res_id, rettype="text") 
  # fetch this entry's summary info
  res_summary <- entrez_summary(db='nuccore',id = res_id) %>% 
    extract_from_esummary(., c("title", "extra", "taxid", 'subtype', 'organism', 
                               'gi'), simplify = FALSE)
  # res_summary[[1]]
  summary_table <- res_summary %>% map_dfr(.f = as_tibble)
  # summary_table$extra[2]
  
  # create a table with the transcript info
  sequence_info <- tibble(id = res_id, accession=res_summary$extra,
                        title= res_summary$title, taxid = res_summary$taxid, 
                        organism = res_summary$organism, gi = res_summary$gi,
                        fasta = res_fasta, feature = res_summary$subtype)
  return(sequence_info)
}

# use web history!! (see https://docs.ropensci.org/rentrez/articles/rentrez_tutorial.html#get-a-web_history-object-from-entrez_search-or-entrez_link)

get_NCBI_info_from_web_history <- function(web_history_obj, rec_start, chunk_size=100){
  # rec_start = 1
  # chunk_size=50
  # web_history_obj = res$web_history
  # fetch this entry as FASTA
  res_fasta <- entrez_fetch(db='nuccore', rettype = c("fasta")  , # id = res_id, 
                            web_history = web_history_obj,
                            retmax=chunk_size, retstart=rec_start)
  # convert into a vector
  fasta_vector <- str_split(res_fasta, pattern = ">") %>% unlist() %>% 
    stringi::stri_remove_empty_na() %>% paste0(">", .)
  # res_features <- entrez_fetch(db="nuccore",  rettype="text", # id=res_id,
  #                              web_history = web_history_obj , 
  #                              retmax=chunk_size, retstart=rec_start) 
  # fetch this entry's summary info
  res_summary <- entrez_summary(db='nuccore',web_history = web_history_obj, #  id = res_id, 
                                retmax=chunk_size, retstart=rec_start) %>% 
    extract_from_esummary(., c("title", "extra", "taxid", 'subtype', 'organism', 
                               'gi', 'uid'), simplify = FALSE)
  
  # res_summary[[1]]
  summary_table <- res_summary %>% 
    map_dfr(.f = ~as_tibble(.x) %>% mutate(across(everything(), .fns = as.character)))
  # summary_table$extra[2]
  
    # create a table with the transcript info
  if (!"subtype" %in% names(summary_table)) summary_table$subtype <- ""
  sequence_info <- summary_table %>% rename(feature = subtype) %>% 
    mutate(fasta = fasta_vector, taxid = as.integer(taxid)) # also process accession out of the extra column
  # sequence_info$organism
  return(sequence_info)
}


# use furrr and progressr to process these in parallel, see https://furrr.futureverse.org/articles/progress.html
plan(multisession, workers = min(4, availableCores())) # adjust cores based on your computer (but too many will cause the server to reject the requests)

# run in parallel (with progress bar)
# rerun command if it fails (see https://purrr.tidyverse.org/reference/insistently.html)
rate <- rate_delay(0.3, max_times = 10) # introduce a delay between queries to not to overload the server
insistent_retrieve_info <- insistently(get_NCBI_info_from_web_history, rate = rate, quiet = FALSE)
#insistent_retrieve_info <- insistently(get_res_then_NCBI_info_from_web_history, rate = rate, quiet = FALSE)
# process all ids from all species
# test_input <- head(taxa, 50)
chunk_size=100


# # process results

test_res <- seq(1,res$count,chunk_size) %>% 
  map_dfr(~get_NCBI_info_from_web_history(web_history_obj = res$web_history,
                                           rec_start = .x, chunk_size = chunk_size))


#NCBI recommends that users post no more than three URL requests per second and limit large jobs to either weekends or between 9:00 PM and 5:00 AM Eastern time during weekdays.
results_data <- chunked_search_strings$chunked_search_strings[1] %>% 
  imap_dfr(.f = ~{
    index = .y
    species_name = .x
    query <- glue('{searchdef} "{species_name}"[Organism]')
    res <- entrez_search(db="nuccore",
                         term=query,
                         retmax=5000)
    LogMsg(glue("Processing search query for '{species_name}' ({index}/{nrow(test_input)}), please wait..."))
    with_progress({
      p <- progressor(steps = length(res$ids))
      
      results_table <- res$ids %>% 
        future_map_dfr(.f = ~{
          id_res <- insistent_retrieve_info(.x)
          p()
          return(id_res)
        }) 
    })
    # LogMsg(glue("Finished processing search query for '{species_name}', waiting 1 seconds before the next one"))
    Sys.sleep(0.5) # introduce a delay between queries to not to overload the server
    results_table %>% mutate(species_name = species_name, query = query)
  
})

# save the results
save(results_data, file = filedate(filename = "ITS_NCBI_data", ext = ".RData"))

# load the results (so you can continue straight from here)
recent_data_file <- recent_file(".", "ITS_NCBI_data")
load(recent_data_file)

# perform your filters...
# Step 3: filtering ---------------------------------------------------------------
#nucleotide.summary <- entrez_summary(db="nucleotide",id=data$ids) #see variables to consider
# results_data = all_results_table
filtered_data_with_voucher <- results_data %>% dplyr::filter(grepl('specimen_voucher', feature)) %>% 
  dplyr::filter(!grepl('internal transcribed spacer 2, partial sequence|internal transcribed spacer 2-like', 
                       perl = TRUE,  title))

# Step 5: Get taxa info -----------------------------------------------------------
# BTW, this can be done on your original species list from your CSV file
species_list <- (unique(filtered_data_with_voucher$taxid))
# prepare taxonomizr database (will download the db locally)
dir.create("tax_db")
tax_db <- prepareDatabase('tax_db/accessionTaxa.sql', getAccessions=FALSE)
# tax_db <- normalizePath('C:/Users/idoid/Griffith University/QX Disease in SRO - Documents/General/Niki_Nenadic/Niki_PhD/RNA-Seq/TransPi/Annotation_analysis/tax_db/accessionTaxa.sql')
# normalize taxa
taxonomy_table <- getTaxonomy(species_list, tax_db,
                              desiredTaxa = c("superkingdom", "kingdom", "phylum", "class", "clade", 
                                              "order", "family", "genus",
                                              "species"))  %>% 
  as.data.frame() %>% 
  mutate(taxid=rownames(.)) %>% 
  as_tibble() %>% remove_empty()


# This is where I got so far, I can then use the taxonomy classification to create the fasta headers

# Step 6: save sequences as FASTA ------------------------------

# save the sequences as fasta (this uses the original fasta headers)
fasta_output <- filedate(filename = "filtered_ITS_sequences", ext = ".fasta", outdir = "fasta")
walk(filtered_data_with_voucher$fasta, .f = ~cat(.x, file = fasta_output, append = TRUE))



























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
# with taxonomizr package you can do this work offline (after it downloads the entire taxonomy database locally)
# see https://github.com/sherrillmix/taxonomizr
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




# testing on 27/04/2023 ---------------------------------------------------
species_chunk_size <- 50

process_taxa <- function(taxa_list){
  # create test dataset
  #taxa_subset <- head(taxa, n = 50)
  taxa_search_string <- paste0(searchdef, " (",
                               paste(glue('"{taxa_list}"[Organism]'), collapse = " OR "), ")")
  # try to retrieve all search results for all species at once!!
  res <- entrez_search(db="nuccore",
                       term=taxa_search_string,
                       use_history = TRUE)
  return(res) 
}

taxa <- read_csv("input_data/records-2023-03-16.csv") %>% select(`Species Name`) %>% 
  mutate(genebank.search=glue('{searchdef} "{`Species Name`}"[Organism]')) 

# use furrr and progressr to process these in parallel, see https://furrr.futureverse.org/articles/progress.html
plan(multisession, workers = min(4, availableCores())) # adjust cores based on your computer (but too many will cause the server to reject the requests)

# run in parallel (with progress bar)
# rerun command if it fails (see https://purrr.tidyverse.org/reference/insistently.html)
rate <- rate_delay(0.3, max_times = 10) # introduce a delay between queries to not to overload the server
insistent_retrieve_info <- insistently(get_NCBI_info_from_web_history, rate = rate, quiet = FALSE)
# insistent_retrieve_info <- insistently(process_taxa, rate = rate, quiet = FALSE)

history_records_chunk_size=100
taxa_chunk_size <- 50

# taxa_chunks <- seq(1, nrow(taxa), species_chunk_size) 
# split a dataframe into subsets by size
split_df <- function(df, group_size){
  split(df, gl(ceiling(nrow(df)/group_size), group_size, nrow(df)))
}
# taxa_subset <- taxa[1:1000,]
# taxa_chunks <- split_df(taxa_subset, 50)
taxa_chunks <- split_df(taxa, taxa_chunk_size)
# iter_sleep <- 0.5 # delay between chunks

tic() # start timer (run together with the code that follows until and including the toc() at the end to stop the timer)
LogMsg(glue("Processing search queries for {nrow(taxa)} species (in {length(taxa_chunks)} chunks of {taxa_chunk_size} taxa), please wait..."))
with_progress({
      p <- progressor(steps = length(taxa_chunks))
      results_table <- taxa_chunks %>% 
        future_imap_dfr(.f = ~{
          # making sure that global variables are available within the future (parallel) function
          # taxa_chunk_size = taxa_chunk_size 
          # taxa_subset = taxa_subset
          # index = as.integer(.y) #chunks of list being processed
          taxa_chunk <- .x
          # taxa_chunk <- taxa_chunks[[1]]
         # LogMsg(glue("Processing search queries ({index}/{length(taxa_chunks)}), please wait..."))
          chunk_res <- process_taxa(taxa_list = taxa_chunk$`Species Name`)
          res <- seq(1, chunk_res$count, history_records_chunk_size) %>% 
            map_dfr(~insistent_retrieve_info(web_history_obj = chunk_res$web_history,
                                                    rec_start = .x, chunk_size = history_records_chunk_size))
          p()
          # LogMsg(glue("Finished processing taxa chunk ({index*taxa_chunk_size}/{nrow(taxa_subset)})"))
          return(res)
        }) 
    
    # Sys.sleep(iter_sleep) # introduce a delay between queries to not to overload the server
   # results_table %>% mutate(species_name = species_name, query = query)
  })
LogMsg(glue("Finished processing taxa chunks. Overall time was:"))
toc() # stop timer
# save the results
save(results_table, file = filedate(filename = "ITS_NCBI_data_all_taxa_results", ext = ".RData"))