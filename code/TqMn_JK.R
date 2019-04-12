# TaqMan data

library(tidyverse)
library(readxl)

# Load TaqMan data
Taqman_results_V1 <- read_excel("data/raw/Taqman_results.xlsx", 
                             sheet = "Taqman_visit1")

Taqman_results_V5 <- read_excel("data/raw/Taqman_results.xlsx", 
                             sheet = "Taqman_visit5")

# Select data of interest
stool_tqmn_V1 <- select(Taqman_results_V1,   
                             starts_with("Bacterial"), starts_with("AGE"), starts_with("CTX"), starts_with("KPC"), 
                             starts_with("NDM"), starts_with("SHV"), starts_with("TEM"), starts_with("CMY"), starts_with("STUDY"))

stool_tqmn_V5 <- select(Taqman_results_V5,   
                        starts_with("Bacterial"), starts_with("AGE"), starts_with("CTX"), starts_with("KPC"), 
                        starts_with("NDM"), starts_with("SHV"), starts_with("TEM"), starts_with("CMY"), starts_with("STUDY"))
