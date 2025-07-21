# Multi-Pool Peptide Voting Analysis
# Processes 5 dataset pools through three search engines

library(dplyr)

# Set base directory
# setwd("/Users/juntaozhao/R/PepTrust")
setwd("D:/R/PepTrust")

# Voting function
get_vote <- function(pep_comet, pep_msgf, pep_xtandem) {
  peptides <- c(pep_comet, pep_msgf, pep_xtandem)
  peptides <- peptides[!is.na(peptides)]
  
  if (length(peptides) == 0) {
    return(list(peptide = "", vote = 0))
  }
  
  freq <- table(peptides)
  max_vote <- max(freq)
  
  if (max_vote == 1) {
    selected_peptide <- sample(peptides, 1)
    return(list(peptide = selected_peptide, vote = 0))
  }
  
  winner <- names(freq[freq == max_vote])[1]
  return(list(peptide = winner, vote = max_vote))
}

# # Create output directory
# dir.create("combined_results/voting", recursive = TRUE, showWarnings = FALSE)

# Process each pool
for (pool_num in 1:5) {
  cat("Processing Pool", pool_num, "\n")
  
  # Load data
  data_comet <- read.csv(paste0("searched_results/comet/pool", pool_num, ".csv"))
  data_msgf <- read.csv(paste0("searched_results/msgf/pool", pool_num, ".csv"))
  data_xtandem <- read.csv(paste0("searched_results/tandem/pool", pool_num, ".csv"))
  
  # Prepare for merging
  comet <- data_comet %>% select(scan, Peptide) %>% rename(peptide_comet = Peptide)
  msgf <- data_msgf %>% select(scan, Peptide) %>% rename(peptide_msgf = Peptide)
  xtandem <- data_xtandem %>% select(scan, Peptide) %>% rename(peptide_xtandem = Peptide)
  
  # Merge and vote
  result <- comet %>%
    full_join(msgf, by = "scan") %>%
    full_join(xtandem, by = "scan") %>%
    rowwise() %>%
    mutate(
      vote_result = list(get_vote(peptide_comet, peptide_msgf, peptide_xtandem)),
      plain_peptide = vote_result$peptide,
      vote = vote_result$vote
    ) %>%
    ungroup() %>%
    select(scan, plain_peptide, vote)
  
  # Save results
  write.csv(result, paste0("combined_results/voting_02/pool_", pool_num, ".csv"), row.names = FALSE)
}

cat("All pools processed successfully!\n")