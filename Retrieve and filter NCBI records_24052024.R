# Semi-Automated Retrieval of ITS2 Sequences from GenBank: Generating a Reference Database for Queensland Plant Species
# Cinthia Pietromonaco, Ido Bar, Alison Peel
# Critera for search- 1. Complete sequence 2. Herbarium voucher  3.ITS2 gene 4. Sequence length between 200-2500 

# Set paths
NCBI_key <- "33ec35e1973e4928beafd47bee86492be109"
species_list <- "/Users/cinthiapietromonaco/Downloads/data.csv" #species list downloaded from ALA spatial.

# Load custom functions created by Ido
devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "Util.R")

# Install and load packages
packages <- c("devtools", "tidyverse", "janitor", "glue", "furrr", "progressr", "seqinr", 
              "taxonomizr", "tictoc", 'rentrez')
pacman::p_load(char = basename(packages))

# Set NCBI API key - this allows you to send more queries to Entrez without timing out - see more details at https://support.nlm.nih.gov/knowledgebase/article/KA-05318/en-us and https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us
set_entrez_key(NCBI_key)

# Define search parameters
searchdef <- '((internal transcribed spacer[All Fields]) OR (ITS[All Fields])) AND (plants[filter]) AND (200[SLEN] : 2500[SLEN]) AND (flowering plants[porgn]) AND'

# Load and preprocess data
taxa <- species_list %>%
  filter(Kingdom == "Plantae") %>% 
  slice(1:1000) %>%
  select(`Species Name`) %>%
  mutate(genebank.search = glue('{searchdef} "{`Species Name`}"[Organism]'))

# Function to split dataframe into chunks
split_df <- function(df, group_size) {
  split(df, ceiling(seq_along(df[[1]]) / group_size))
}

# Split data into chunks
taxa_chunks <- split_df(taxa, 30)

# Function to process taxa
process_taxa <- function(taxa_list) {
  taxa_search_string <- paste0(searchdef, " (", paste(glue('"{taxa_list}"[Organism]'), collapse = " OR "), ")")
  res <- entrez_search(db = "nuccore", term = taxa_search_string, use_history = TRUE)
  return(res)
}

# Function to retrieve NCBI info
get_NCBI_info_from_web_history <- function(web_history_obj, rec_start, chunk_size = 100) {
  res_fasta <- entrez_fetch(db = 'nuccore', rettype = "fasta", web_history = web_history_obj, retmax = chunk_size, retstart = rec_start)
  fasta_vector <- str_split(res_fasta, pattern = ">") %>% unlist() %>% .[-1] %>% paste0(">", .)
  res_summary <- entrez_summary(db = 'nuccore', web_history = web_history_obj, retmax = chunk_size, retstart = rec_start) %>%
    extract_from_esummary(., c("title", "extra", "taxid", 'subtype', 'organism', 'gi', 'uid'), simplify = FALSE)
  summary_table <- map_dfr(res_summary, ~as_tibble(.x) %>% mutate(across(everything(), as.character)))
  if (!"subtype" %in% names(summary_table)) summary_table$subtype <- ""
  if (!"taxid" %in% names(summary_table)) summary_table$taxid <- ""
  sequence_info <- summary_table %>% rename(feature = subtype) %>%
    mutate(fasta = fasta_vector, taxid = as.integer(taxid))
  return(sequence_info)
}

# Set up parallel processing 
plan(multisession, workers = min(2, availableCores()))
rate <- rate_delay(1, max_times = 5)
insistent_retrieve_info <- insistently(get_NCBI_info_from_web_history, rate = rate, quiet = FALSE)

# Process data with progress bar
tic()
LogMsg(glue("Processing search queries for {nrow(taxa)} species (in {length(taxa_chunks)} chunks of {30} taxa), please wait..."))

with_progress({
  p <- progressor(steps = length(taxa_chunks))
  results_table <- future_imap_dfr(taxa_chunks, .f = ~{
    tryCatch({
      chunk_res <- process_taxa(.x$`Species Name`)
      res <- seq(0, chunk_res$count - 1, by = 100) %>%
        map_dfr(~insistent_retrieve_info(web_history_obj = chunk_res$web_history, rec_start = .x, chunk_size = 100))
      p()
      return(res)
    }, error = function(e) {
      cat("Error:", conditionMessage(e), "\n")
      return(NULL)
    })
  })
})
LogMsg("Finished processing taxa chunks. Overall time was:")
toc()

# Filter and format the results
filtered_data_with_voucher <- results_table %>%
  filter(grepl('specimen_voucher', feature)) %>%
  filter(!grepl('internal transcribed spacer 2, partial sequence|internal transcribed spacer 2-like|ITS2 \\(partial\\)', title, perl = TRUE)) %>% 
  filter(grepl('internal transcribed spacer 2|ITS2', title))

accession <- data.frame(accession = gsub('>', '', as.character(sapply(strsplit(filtered_data_with_voucher$fasta, " "), "[", 1))))
sequence <- data.frame(sequence = str_sub(word(filtered_data_with_voucher$fasta, -1), 9, -1))
with_sequence <- cbind(filtered_data_with_voucher, accession, sequence)

# Section below is for missing taxonomy
# Manually download taxdump.tar.gz file from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/
# can also use function below if connection is stable from taxonomizr package.
# prepareDatabase('accessionTaxa.sql')

# Prepare the NCBI taxonomy database manually
nodes_file <- "/Users/cinthiapietromonaco/Library/CloudStorage/OneDrive-GriffithUniversity/Metabarcoding (Shared CP)/Reference Database (Shared CP)/Peel_get_NCBI_records/NCBI Taxonomy Database/taxdump/nodes.dmp"
names_file <- "/Users/cinthiapietromonaco/Library/CloudStorage/OneDrive-GriffithUniversity/Metabarcoding (Shared CP)/Reference Database (Shared CP)/Peel_get_NCBI_records/NCBI Taxonomy Database/taxdump/names.dmp"

# Create the taxonomy database 
db <- "accessionTaxa.sql"
read.nodes.sql(nodes_file, db)
read.names.sql(names_file, db)

# Function to get taxonomic information
get_taxonomic_info <- function(species_name) {
  tax_ids <- getId(species_name, db)
  if (length(tax_ids) == 0) {
    return(data.frame(Kingdom = NA, Phylum = NA, Order = NA, Family = NA, Genus = NA))
  }
  tax_info <- getTaxonomy(tax_ids, db)
  colnames(tax_info) <- c("Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax_info <- tax_info[, c("Superkingdom", "Phylum", "Order", "Family", "Genus")]
  return(tax_info)
}

# Replace any taxonomy missing 
taxa <- species_list %>%
  select(organism = `Species Name`, Kingdom, Phylum, Order, Family, Genus) %>%
  mutate(first_two_words = word(organism, 1, 2))

taxa <- taxa %>%
  replace_na(list(Kingdom = "NA", Phylum = "NA", Order = "NA", Family = "NA", Genus = "NA"))

with_sequence <- with_sequence %>%
  mutate(first_two_words = word(organism, 1, 2))

with_taxa <- merge(with_sequence, taxa, by = "first_two_words", all.x = TRUE) %>%
  mutate(organism_full = organism.x)

missing_taxa <- with_taxa %>%
  filter(is.na(Kingdom) | is.na(Phylum) | is.na(Order) | is.na(Family) | is.na(Genus) |
           Kingdom == "NA" | Phylum == "NA" | Order == "NA" | Family == "NA" | Genus == "NA") %>%
  select(organism_full) %>%
  distinct()

taxonomic_info <- missing_taxa %>%
  rowwise() %>%
  mutate(tax_info = list(tryCatch(get_taxonomic_info(organism_full), 
                                  error = function(e) return(data.frame(Kingdom = NA, Phylum = NA, Order = NA, Family = NA, Genus = NA))))) %>%
  unnest_wider(tax_info) %>%
  select(organism_full, Superkingdom, Phylum, Order, Family, Genus) %>% 
  rename(Kingdom = Superkingdom)

for (i in 1:nrow(taxonomic_info)) {
  species_name <- taxonomic_info$organism_full[i]
  
  with_taxa <- with_taxa %>%
    mutate(
      Kingdom = ifelse(is.na(Kingdom) | Kingdom == "NA" & organism_full == species_name, taxonomic_info$Kingdom[i], Kingdom),
      Phylum = ifelse(is.na(Phylum) | Phylum == "NA" & organism_full == species_name, taxonomic_info$Phylum[i], Phylum),
      Order = ifelse(is.na(Order) | Order == "NA" & organism_full == species_name, taxonomic_info$Order[i], Order),
      Family = ifelse(is.na(Family) | Family == "NA" & organism_full == species_name, taxonomic_info$Family[i], Family),
      Genus = ifelse(is.na(Genus) | Genus == "NA" & organism_full == species_name, taxonomic_info$Genus[i], Genus)
    )
}


# Create FASTA file with taxonomic information
with_link <- with_taxa %>%
  mutate(organism_full = gsub(" ", "_", organism.x),
         link_name = glue("{sub('\\\\.1$', '', accession)}; tax=k:{Kingdom},p:{Phylum},o:{Order},f:{Family},g:{Genus},s:{organism_full};"), sequence = sub("^ast", "", sequence))

write.fasta(as.list(with_link$sequence), names = with_link$link_name, file.out = "/Users/cinthiapietromonaco/Desktop/its2.ffdietfromlit.2024-01.fasta")
