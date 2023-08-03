taxa <- read_csv("input_data/records-2023-03-16.csv") %>% select(`Species Name`) %>% 
  mutate(genebank.search=glue('{searchdef} "{`Species Name`}"[Organism]')) 

# split a dataframe into subsets by size (now included in my Utils.R file that is loaded at the beginning)
split_df <- function(df, group_size){
  split(df, gl(ceiling(nrow(df)/group_size), group_size, nrow(df)))
}
# taxa_subset <- taxa[1:1000,]
# taxa_chunks <- split_df(taxa_subset, 50)
history_records_chunk_size=50
taxa_chunk_size <- 50
taxa_chunks <- split_df(taxa, taxa_chunk_size)

# use furrr and progressr to process these in parallel, see https://furrr.futureverse.org/articles/progress.html
plan(multisession, workers = min(4, availableCores()-1)) # adjust cores based on your computer (but too many will cause the server to reject the requests)
# run in parallel (with progress bar)
# rerun command if it fails (see https://purrr.tidyverse.org/reference/insistently.html)
rate <- rate_delay(0.3, max_times = 10) # introduce a delay between queries to not to overload the server
insistent_retrieve_info <- insistently(get_NCBI_info_from_web_history, rate = rate, 
                                       quiet = FALSE)

history_records_chunk_size=50
taxa_chunk_size <- 50


taxa_chunk <- taxa_chunks[[280]]
tail(taxa_chunk)
chunk_res <- process_taxa(taxa_list = taxa_chunk$`Species Name`)

process_taxa_chunk <- function(taxa_chunk){
  chunk_res <- process_taxa(taxa_list = taxa_chunk$`Species Name`)
  
  if (chunk_res$count %% history_records_chunk_size == 1) {
    history_records_chunk_size <- history_records_chunk_size+2
  }
  history_chunk_indices <- seq(1, chunk_res$count, history_records_chunk_size)
  
 # tic() # start timer (run together with the code that follows until and including the toc() at the end to stop the timer)
  # LogMsg(glue("Processing search queries for {nrow(taxa_chunk)} species (in {length(history_chunk_indices)} records of web_history), please wait..."))
  with_progress({
    p <- progressor(steps = length(history_chunk_indices))
    chunk_results <- history_chunk_indices %>% 
      future_map_dfr(.f = ~{
        # making sure that global variables are available within the future (parallel) function
        # taxa_chunk_size = taxa_chunk_size 
        # taxa_subset = taxa_subset
        # index = as.integer(.y) #chunks of list being processed
        
        # taxa_chunk <- taxa_chunks[[1]]
        # LogMsg(glue("Processing search queries ({index}/{length(taxa_chunks)}), please wait..."))
        
        res <-  insistent_retrieve_info(web_history_obj = chunk_res$web_history,
                                        rec_start = .x, chunk_size = history_records_chunk_size)
        p()
        # LogMsg(glue("Finished processing taxa chunk ({index*taxa_chunk_size}/{nrow(taxa_subset)})"))
        return(res)
      }) 
    
    # Sys.sleep(iter_sleep) # introduce a delay between queries to not to overload the server
    # results_table %>% mutate(species_name = species_name, query = query)
  })
  # LogMsg(glue("Finished processing taxa chunk. Overall time was:"))
  # toc() # stop timer
  return(chunk_results)
}
# test one chunk
test_results <- process_taxa_chunk(taxa_chunks[[219]])


# try running a large set
history_records_chunk_size=50
taxa_chunk_size <- 50
taxa_chunks <- split_df(taxa, taxa_chunk_size)
taxa_chunks2 <- taxa_chunks[51:60]
tic()
all_results <-  taxa_chunks %>% 
  imap_dfr(.f = ~{
    
    LogMsg(glue("Processing search queries for {taxa_chunk_size} species (total {taxa_chunk_size*as.integer(.y)}/{nrow(taxa)}), please wait..."))
    process_taxa_chunk(taxa_chunk=.x)
    })
toc()