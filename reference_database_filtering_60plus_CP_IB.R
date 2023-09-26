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
packages <- c("tidyverse", "janitor", "glue", "furrr", "progressr", "seqinr", 
              "taxonomizr", "tictoc", 'allenbaron/rentrez')
pak::pak(packages) # install packages

pacman::p_load(char = basename(packages))

# set NCBI API key - this allows yo to send more queries to Entrez without timing out - see more details at https://support.nlm.nih.gov/knowledgebase/article/KA-05318/en-us and https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us
set_entrez_key("33ec35e1973e4928beafd47bee86492be109")
# Step 1: Load data and requirements --------------------------------------
#string 1 for ITS2 marker
searchdef <- '((internal transcribed spacer[All Fields]) OR (ITS[All Fields])) AND (plants[filter]) AND (200[SLEN] : 2500[SLEN]) AND (flowering plants[porgn]) AND'

species_name <- "SPECIES_NAME"
taxa <- read_csv("input_data/records-2023-09-26.csv")[15050:17000,] %>% select(`Species Name`) %>% 
  mutate(genebank.search=glue('{searchdef} "{`Species Name`}"[Organism]')) 

# split a dataframe into subsets by size (now included in my Utils.R file that is loaded at the beginning)
taxa_chunk_size <- 50
history_records_chunk_size=100

split_df <- function(df, group_size) {
  num_rows <- nrow(df)
  num_chunks <- ceiling(num_rows / group_size)
  rows_needed <- num_chunks * group_size
  
  if (rows_needed > num_rows) {
    last_row <- df[num_rows, ]
    df <- rbind(df, last_row)
    num_rows <- nrow(df)
  }
  
  chunks <- split(df, rep(1:num_chunks, each = group_size, length.out = num_rows))
  return(chunks)
}

taxa_chunks <- split_df(taxa, taxa_chunk_size)
# Step 2: Load functions --------------------------------------------------
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
get_NCBI_info_from_web_history <- function(web_history_obj, rec_start, chunk_size=100){
  # rec_start = 1
  # chunk_size=50
  # web_history_obj = res$web_history
  # fetch this entry as FASTA
  res_fasta <- entrez_fetch(db='nuccore', rettype = c("fasta")  , # id = res_id, 
                            web_history = web_history_obj,
                            retmax=chunk_size, retstart=rec_start)
  # convert into a vector
  fasta_vector <- str_split(res_fasta, pattern = ">") %>% unlist() %>% .[-1] %>% 
    paste0(">", .)
  # res_features <- entrez_fetch(db="nuccore",  rettype="text", # id=res_id,
  #                              web_history = web_history_obj , 
  #                              retmax=chunk_size, retstart=rec_start) 
  # fetch this entry's summary info
  res_summary <- entrez_summary(db='nuccore', web_history = web_history_obj, #  id = res_id, 
                                retmax=chunk_size, retstart=rec_start) %>% 
    extract_from_esummary(., c("title", "extra", "taxid", 'subtype', 'organism', 
                               'gi', 'uid'), simplify = FALSE)
  
  # res_summary[[1]]
  summary_table <- res_summary %>% 
    map_dfr(.f = ~as_tibble(.x) %>% mutate(across(everything(), .fns = as.character)))
  # summary_table$extra[2]
  
  # create a table with the transcript info
  if (!"subtype" %in% names(summary_table)) summary_table$subtype <- ""
  if (!"taxid" %in% names(summary_table)) summary_table$taxid <- ""
  sequence_info <- summary_table %>% rename(feature = subtype) %>% 
    mutate(fasta = fasta_vector, taxid = as.integer(taxid)) # also process accession out of the extra column
  # sequence_info$organism
  return(sequence_info)
}
# use furrr and progressr to process these in parallel, see https://furrr.futureverse.org/articles/progress.html
plan(multisession, workers = min(4, availableCores())) # adjust cores based on your computer (but too many will cause the server to reject the requests)

# run in parallel (with progress bar)
# rerun command if it fails (see https://purrr.tidyverse.org/reference/insistently.html)
rate <- rate_delay(1, max_times = 5) # introduce a delay between queries to not to overload the server
insistent_retrieve_info <- insistently(get_NCBI_info_from_web_history, rate = rate, quiet = FALSE)

# Step 3: Process ---------------------------------------------------------
tic() # start timer (run together with the code that follows until and including the toc() at the end to stop the timer)
LogMsg(glue("Processing search queries for {nrow(taxa)} species (in {length(taxa_chunks)} chunks of {taxa_chunk_size} taxa), please wait..."))
with_progress({
  p <- progressor(steps = length(taxa_chunks))
  results_table10 <- taxa_chunks %>%
    future_imap_dfr(.f = ~{
      tryCatch({
        # making sure that global variables are available within the future (parallel) function
        # taxa_chunk_size = taxa_chunk_size 
        # taxa_subset = taxa_subset
        # index = as.integer(.y) #chunks of list being processed
        taxa_chunk <- .x
        # taxa_chunk <- taxa_chunks[[1]]
        # LogMsg(glue("Processing search queries ({index}/{length(taxa_chunks)}), please wait..."))
        chunk_res <- process_taxa(taxa_list = taxa_chunk$`Species Name`)
        res <- seq(0, chunk_res$count - 1, history_records_chunk_size) %>%
          map_dfr(~insistent_retrieve_info(web_history_obj = chunk_res$web_history,
                                           rec_start = .x, chunk_size = history_records_chunk_size))
        p()
        # LogMsg(glue("Finished processing taxa chunk ({index*taxa_chunk_size}/{nrow(taxa_subset)})"))
        return(res)
      }, error = function(e) {
        # Handle the error (e.g., print an error message)
        cat("Error:", conditionMessage(e), "\n")
        return(NULL)  # Return a placeholder value or NULL to indicate an error
      })
    })
})
LogMsg(glue("Finished processing taxa chunks. Overall time was:"))
toc() # stop timer
# Step 4: Merge/Filter/Format ---------------------------------------------

data_frames <- list(
  results_table, 
  results_table2, 
  results_table3, 
  results_table4, 
  results_table5, 
  results_table6, 
  results_table7, 
  results_table8, 
  results_table9.1,
  results_table9.2,
  results_table10
)

# Merge the data frames by the first 8 columns
merged_data <- Reduce(function(x, y) merge(x, y, by = 1:8, all = TRUE), data_frames)

filtered_data_with_voucher <- merged_data %>% dplyr::filter(grepl('specimen_voucher', feature)) %>% 
  dplyr::filter(!grepl('internal transcribed spacer 2, partial sequence|internal transcribed spacer 2-like', 
                       perl = TRUE,  title))

accession <- data.frame(accession = as.character(sapply(strsplit(filtered_data_with_voucher$fasta, " "), "[", 1)))
accession$accession <- gsub('>','',accession$accession)

with_accession <- cbind(filtered_data_with_voucher, accession)

sequence <- data.frame(sequence = as.character(word(with_accession$fasta,-1)))
sequence <- sequence %>%
  mutate(sequence = str_sub(sequence, 9, -1))
with_sequence <- cbind(with_accession, sequence)

taxa <- read_csv("input_data/records-2023-09-26.csv") %>% 
  select(organism = `Species Name`, Kingdom:Genus)

with_taxa <- merge(with_sequence, taxa, by = "organism")


with_link <- with_taxa %>%
  mutate(organism = gsub(" ", "_", organism),
         link_name = gsub(" ", "", paste(sub("\\.1$", "", accession), "; tax=k:", Kingdom,
                                         ",p:", Phylum, ",c:", Class, ",o:", Order,
                                         ",f:", Family, ",g:", Genus, ",s:", organism, ";")))


# Step 5 - Create FASTA file ----------------------------------------------
library(seqinr)
write.fasta(as.list(with_link$sequence),names = with_link$link_name, file.out = "/Users/cinthiapietromonaco/Desktop/its2.didijustdothis.tax.fasta")
