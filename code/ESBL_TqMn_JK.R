# Load libraries
library(tidyverse)
library(readxl)
library(xtable)

# Load data.  This data is from TrEAT TD trial and contains observations from 366
# Where do the 366 come from? TrEAT only randomized 339. 
# The numbers given are CT/Cq values i.e. # of cycles of PCR necessary in order to 
# reach threshold of detection. We deem 35 or less to indicate the genetic material 
# present (value is low enough that we deem it to be present)
# set working directory to gut_MPH

######################################
#### Both TaqMan & Phenotypic data ###
######################################

treat <- read_excel("data/raw/TrEAT_Merge_ESBL_2018.09.13_v2.XLSX")

str(treat)
# Get a conception of the data set:
ls()  # contents of current work space
class(treat) # dataframe
dim(treat)  
print(object.size(treat), units = "Mb") # just over 1 Mb

# Average TLUS by site (using dplyr). Remove NA's from columns
#   and include a column with number of instances per group
treat %>% 
  group_by(SITE) %>% 
  summarise(mean_TLUS = mean(TLUS, na.rm = TRUE),
            n = n())


########################
#### Phenotypic data ###
########################

# Let's narrow the data we want to observe:
Slim_treat <- treat %>%
  select(STUDY_ID_TRUNC, AGE,
         SEX, starts_with("ESBL_V"))

str(Slim_treat)  # PROBLEM:   R sees phenotypic resistance status (ESBL_V1, ESBL_V5) as
# character data instead of numeric/integer.  
# Must convert these data to numeric/integer. 

# Create new df with ESBL_V1/V5 variables as integer instead of character/factor:
# Code as: 
# NS -> NA
# Negative -> 0
# Positive -> 1
# N/A -> NA

int_ESBL_V1V5 <- Slim_treat %>%
  mutate(int_V1 = ifelse(ESBL_V1 == "NS", NA, 
                  ifelse(ESBL_V1 == "Negative", 0, 
                  ifelse(ESBL_V1 == "Positive", 1, 
                  ifelse(ESBL_V1 == "N/A", NA, NA))))) %>%
  mutate(int_V5 = ifelse(ESBL_V5 == "NS", NA, 
                  ifelse(ESBL_V5 == "Negative", 0, 
                  ifelse(ESBL_V5 == "Positive", 1, 
                  ifelse(ESBL_V5 == "N/A", NA, NA)))))



str(int_ESBL_V1V5)

# Find number of missing in int_V1 variable:
sum(is.na(int_ESBL_V1V5$int_V1)) # 134/366
# missing in int_V5:
sum(is.na(int_ESBL_V1V5$int_V5)) # 196/366

# More visual methods: 
na_V1 <- is.na(int_ESBL_V1V5$int_V1)
table(na_V1)

int_ESBL_V1V5 %>%
  select(int_V1) %>% 
  group_by(int_V1) %>% 
  tally()

int_ESBL_V1V5 %>%
  select(int_V5) %>% 
  group_by(int_V5) %>% 
  tally()

# Summary: 
# Using phenotypic testing we find that there are 12/366 observations that were positive
# among non-pathogenic E. coli on V1 (visit 1) for ESBL and 17/170 positive on 
# V5 (visit day 5)




########################
#### TaqMan data #######
########################
# Select the variables of interest:
#note: this treat dataset only contains TaqMan data for V1; must go to HumiChip dataset 
# to get V5 TaqMan data. (thus the V1 in the variable name)
card_stool_tqmn_V1 <- select(treat,   
                       starts_with("Bacterial"), starts_with("AGE"), starts_with("CTX"), starts_with("KPC"), 
                       starts_with("NDM"), starts_with("SHV"), starts_with("TEM"), starts_with("CMY"), starts_with("STUDY"))

##### Problem: #####
# All these values in this dataframe are coded as factors.  We want them to be
# numeric type. 
# First, what values do these variables contain in terms of our analysis? There are 3 types:
# number(Cq value), "Undetermined" and NA. We interpret Undetermined as meaning 
# there was no bacteria detected using the TaqMan real-time PCR containing the gene 
# coding for the type of ESBL resistance specified. NA means missing data. 

# Let's make these character vectors into doubles. 
# Since the Cq values in this data are unique numbers it's not feasible to use an ifelse 
# statement for each one of them when converting them into numeric type data. Instead, 
# let's designate each "Undetermined" as above the Cq threshold for detection (Undetermined
# designates that during testing we did not detect anything, but it's not a missing/NA 
# value).  We will keep "NA" as NA while leaving the rest as is. In this process we will
# recode Undetermined and Cq numbers from factor to double.



# Use the scoped forms of mutate and transmutate such as mutate_if in order to change
# from one type of variable to another type (thanks: https://dplyr.tidyverse.org/reference/mutate_all.html)
# Also, using Ryan's code, "Undetermined" or "NA" will be NA while the Cq values will
# become doubles. 
tqmn_na_int <- card_stool_tqmn_V1 %>%  # Not the best method since leaves both undetected
  # values and missing values as NA. 
  mutate_at(vars(-STUDY_ID_TRUNC), funs(as.double))  # note: introduces NAs by coercion,
# including all Ryan's code for interpreting the taq values (there are a few "Indeterminate"
# values which we coerce to NA)

# Ryan's code: 
tqmn_na_int_Cq <- card_stool_tqmn_V1 %>%
  # Change 'Undertermined' to max threshold (40)
  mutate_at(vars(-STUDY_ID_TRUNC), funs(ifelse(. == "Undetermined", 40, .))) %>%
  # Convert all columns to numeric (will coerce 'Indeterminate' to NA)
  mutate_at(vars(-STUDY_ID_TRUNC), as.numeric) %>%
  # Convert all values under 35 to 1 and over 35 to 0
  mutate_at(vars(-STUDY_ID_TRUNC), funs(ifelse(. < 35, 1, 0)))

# Note: funs() has been soft deprecated as of dplyr 0.8.0. Instead, use list()
# funs(name = f(.))   now use: list(name = ~f(.)) 
tqmn_na_int_Cq_list <- card_stool_tqmn_V1 %>%
  # Change 'Undertermined' to max threshold (40)
  mutate_at(vars(-c(STUDY_ID_TRUNC, Bacterial_16s_STOOL, Bacterial_16s_CARD)), 
                 list(~ifelse(. == "Undetermined", 40, .))) %>%
  # Convert all columns to numeric (will coerce 'Indeterminate' to NA)
  mutate_at(vars(-c(STUDY_ID_TRUNC, Bacterial_16s_STOOL, Bacterial_16s_CARD)), as.numeric) %>%
  # Convert all values under 35 to 1 and over 35 to 0
  mutate_at(vars(-c(STUDY_ID_TRUNC, Bacterial_16s_STOOL, Bacterial_16s_CARD)), 
            list(~ifelse(. < 35, 1, 0)))


# Now obtain proportion of positive probes among total population within each probe category.

# We want to look at the correspondence bw Card vs Stool.  
# First collect all Card data from all types (i.e. TEM, KPC, etc.) into a list and all
# stool data from all types into a single list and then, for those pairs for which we have
# data for both card and stool, perform the McNemar test. 

CARD_list <- unite(tqmn_na_int_Cq_list, all_CARD_data, c(CTX_CARD, KPC_CARD, 
                                                         NDM_CARD, SHV_CARD, TEM_CARD, 
                                                         CMY_CARD), remove = TRUE)  # Not
# quite the result intended: gave each observation as all 6 variable values together

CARD_list <- unite(tqmn_na_int_Cq_list, all_CARD_data, c(CTX_CARD, KPC_CARD, 
                                                         NDM_CARD, SHV_CARD, TEM_CARD, 
                                                         CMY_CARD), remove = FALSE) # Nope;
# just retains the original 6 variable columns. 
STOOL_list <- unite(tqmn_na_int_Cq_list, all_STOOL_data, c(CTX_STOOL, KPC_STOOL, NDM_STOOL, SHV_STOOL, TEM_STOOL, CMY_STOOL) ), sep = "_", remove = TRUE)

# Try using rbind:
CARD_list <-tqmn_na_int_Cq_list %>%
  rbind("CTX_CARD", "KPC_CARD", "NDM_CARD", "SHV_CARD", "TEM_CARD", "CMY_CARD") # Nope;
# yields 16 columns and 372 observations, appending these 6 string to the bottom of 
# each variable.

CARD_list <- rbind(select(tqmn_na_int_Cq_list, "CTX_CARD", "KPC_CARD", "NDM_CARD", 
                          "SHV_CARD", "TEM_CARD", "CMY_CARD")) # Nope; gives df w/ these.

# Use match:
CTX_match <- match(tqmn_na_int_Cq_list$CTX_CARD, tqmn_na_int_Cq_list$CTX_STOOL)
CTX_match <- (factor(tqmn_na_int_Cq_list$CTX_CARD) %in% factor(tqmn_na_int_Cq_list$CTX_STOOL))
sum(CTX_match)
str(CTX_match)  # This gives 100% matching which we know is not true by looking at the data. 

# Use count:
C_S_corr_CTX <- tqmn_na_int_Cq_list %>%
  mutate(CTX_C_S = paste(CTX_CARD, CTX_STOOL, sep = ",")) %>%
  count(CTX_C_S)

# To make tables of data: 
C_S_corr_CTX_tab <- xtable(C_S_corr_CTX)
print.xtable(C_S_corr_CTX_tab, type="html", file="Card_stool_tables.html")


CTX_matrix <- matrix(c(188, 9, 8, 13), nrow = 2, ncol = 2, byrow = FALSE, dimnames = 
list(c("CTX_S_neg", "CTX_S_pos"), c("CTX_C_neg", "CTX_C_pos")))

# Apply McNemar test to CTX card vs stool data: 
mcnemar.test(CTX_matrix, correct = TRUE)  # Yields p-value = 1, which would tell us that
# probability of seeing this distribution of values would be expected 100% of the time 
# assuming the test hypothesis that number of discordants of each type (1,0 and 0,1) are
# equal, i.e. total row 1 = total column 1 and total row 2 = total column 2
mcnemar.test(CTX_matrix, correct = FALSE)




# CTX
tqmn_na_int_Cq_list %>%
  select(CTX_STOOL) %>% 
  group_by(CTX_STOOL) %>% 
  tally()
tqmn_na_int_Cq_list %>%
  select(CTX_CARD) %>% 
  group_by(CTX_CARD) %>% 
  tally()

# KPC
tqmn_na_int_Cq_list %>%
  select(KPC_STOOL) %>% 
  group_by(KPC_STOOL) %>% 
  tally()
tqmn_na_int_Cq_list %>%
  select(KPC_CARD) %>% 
  group_by(KPC_CARD) %>% 
  tally()

# NDM
tqmn_na_int_Cq_list %>%
  select(KPC_STOOL) %>% 
  group_by(KPC_STOOL) %>% 
  tally()
tqmn_na_int_Cq_list %>%
  select(KPC_CARD) %>% 
  group_by(KPC_CARD) %>% 
  tally()

tqmn_na_int_Cq_list %>%
  select(TEM_STOOL) %>% 
  group_by(TEM_STOOL) %>% 
  tally()

pos_ESBL <- function()
  variable <- vars(-STUDY_ID_TRUNC)
  group_by(variable)
  
##### Next goal: ##### 
# Iterate through all rows of each column and recode each value less than/equal
# to 35 as 1, each value > 35 as 0, and each NA remains NA. 
# Cq value of <= 35 means present. 
presence <- function(cut, dat, na.rm = TRUE){
  cut(dat)
}


presence(cut = cut(tqmn_na_int, all_vars(is_integer(.)), c(0, 35, 100)))
presence(cut = cut(chocolate, c(0, 35, 100)))

tqmn_na_int_pos <- tqmn_na_int %>%
  mutate_if(all_vars(is_integer(.)), ifelse(type))

sums <- summarise_all(tqmn_na_int, sum, na.rm = TRUE)




##### For loops:
  x <- c("a", "b", "c", "d")
for(i in 1:7){
  print(x[i])
}

for(i in rownames(tqmn_na_int))
  print(tqmn_na_int[i, "CTX_STOOL"], na.print = FALSE)

##### Attempt 1
probe_test1 <- function(variable) {
  for(row in seq_len(nrow(tqmn_na_int))) {
  p <- tqmn_na_int[row, variable]
  ifelse(p > 35, 0, ifelse(p <= 35, 1, NA))
  data.frame(p)
  }
}

##### Attempt 2
probe_test2 <- function(variable) {
  for(row in seq_len(nrow(tqmn_na_int))) {
    p <- tqmn_na_int[row, variable]
    mutate(tqmn_na_int, probe_pos = ifelse(p > 35, 0, ifelse(p <= 35, 1, NA)))
  }
}   # Testing with: probe_test2("CTX_STOOL") yields nothing. 

##### Attempt 3
probe_test3 <- function(variable) {
  for(row in seq_len(nrow(tqmn_na_int))) {
    p <- tqmn_na_int[row, variable]
  }  
    mutate(tqmn_na_int, probe_pos = ifelse(p > 35, 0, ifelse(p <= 35, 1, NA)))
}    # Testing with: df_test <- probe_test3("CTX_STOOL") yields df with probe_pos column but
# all are NAs. 

##### Attempt 3.1
probe_test3.1 <- function(variable) {
  for(row in seq_len(nrow(variable))) {
    p <- nrow(variable)
  }  
  mutate(tqmn_na_int, probe_pos = ifelse(p > 35, 0, ifelse(p <= 35, 1, NA)))
} # Testing with: df_test <- probe_test3(tqmn_na_int) yields df with probe_pos column but all 
# are 0s. Duh, bc the number of rows is p and it's greater than 35. 


##### Attempt 4
probe_test4 <- function(variable) {
  for(row in seq_len(nrow(tqmn_na_int))) {
    p <- tqmn_na_int[row, variable]
  }  
  mutate(tqmn_na_int, probe_pos = ifelse(p > 35, 0, ifelse(p <= 35, 1, NA)))
  tibble(p, row)
}   

##### Attempt 5
probe_test5 <- function(x) {
  nc <- ncol(x)
  r <- nrow(x)
  for(i in 1:nc) {
    for(ii in 1:r) {
      mutate(i, probe_pos = ifelse(ii > 35, 0, ifelse(ii <= 35, 1, NA)))    
    }
  }
}    # Testing w/ probe_test5(tqmn_na_int) yielded error since mutate can't take just an integer
# class object.

##### Attempt 6
probe_test6 <- function(x) {
  nc <- ncol(x)
  r <- nrow(x)
  for(i in 1:nc) {
    for(ii in 1:r) {
      mutate(x, probe_pos = ifelse(ii > 35, 0, ifelse(ii <= 35, 1, NA)))    
    }
  }
}  # Renders NULL, likely bc you didn't specify i in first for loop. 


##### Attempt 7
probe_test7 <- function(x) {
  nc <- ncol(x)
  r <- nrow(x)
  clnms <- colnames(x)
  rnms <- list(1:366)
  matx <- matrix(nrow = r, ncol = nc, byrow = FALSE)
  for(i in 1:nc) {
    for(ii in 1:r) {
      matx[i] <- transmute(x, probe_pos = ifelse(ii > 35, 0, ifelse(ii <= 35, 1, NA)))    
    }
  }
}   # Testing w/ df <- probe_test7(tqmn_na_int) yields NULL. 


df <- probe_test7(tqmn_na_int)

g <- card_stool_tqmn$CTX_STOOL
as.numeric(levels(g))[as.integer(f)]
str(g)
g > 30
g[5:50]
g



# To find the correlation or to perform the McNemar test b/w ESBL_V1 and 
# ISOLATES_POSITIVE_V1 they must be numeric. 
# Change ESBL_V1 into numeric. 
num_ESBL <- as.numeric(treat$ESBL_V1) # Now NS coded as 3, 
# Negative coded as 2, positive coded 4 and N/A as 1. 
# However, 
# the as.numeric function has now converted my dataframe into a numeric vector.
# We want Negative = 0 Positive = 1 and N/A as NA, NS as NA.
ESBL_V1_unc <- unclass(num_ESBL)
binary_ESBL_V1 <- ifelse(num_ESBL$ESBL_V1 == 3, NA, 
                         ifelse(ESBL_V1 = 1, NA, 
                                ifelse(ESBL_V1 = 2, 0, 
                                       ifelse(ESBL_V1 = 4, 1, NA)))) # does not work


#### APPENDIX/Self-reference for coding ####

# Some failed attempts to convert a factor vector to numeric using case_when:
int_ESBL_V1V5_c <- Slim_treat %>%
  case_when(
    ESBL_V1 == "NS" ~ 0, 
    ESBL_V1 == "Negative" ~ 0, 
    ESBL_V1 == "Positive" ~ 1, 
    ESBL_V1 == "N/A" ~ NA_integer_
  )

levels(Slim_treat$ESBL_V1)

int_ESBL_V1V5_d <- Slim_treat %>%
  case_when(
    ESBL_V1 == level 3 ~ 0, 
    ESBL_V1 == 2 ~ 0, 
    ESBL_V1 == 4 ~ 1, 
    ESBL_V1 == 1 ~ NA_integer_
  )

int_ESBL_V1V5_e <- ESBL_V1_unc %>%
  case_when(
    ESBL_V1 == 3 ~ 0, 
    ESBL_V1 == 2 ~ 0, 
    ESBL_V1 == 4 ~ 1, 
    ESBL_V1 == 1 ~ NA_integer_
  )

int_ESBL_V1V5_e2 <- ESBL_V1_unc %>%
  case_when(
    ESBL_V1 = 3 ~ 0, 
    ESBL_V1 = 2 ~ 0, 
    ESBL_V1 = 4 ~ 1, 
    ESBL_V1 = 1 ~ NA_integer_
  )

int_ESBL_V1V5_c2 <- Slim_treat %>%
  case_when(
    .$ESBL_V1 == "NS" ~ 0, 
    .$ESBL_V1 == "Negative" ~ 0, 
    .$ESBL_V1 == "Positive" ~ 1, 
    .$ESBL_V1 == "N/A" ~ NA_integer_,
    TRUE ~ as.numeric(as.character(.$ESBL_V1))
  )

new_ESBL_V1 <- Slim_treat %>%
  transmute()



Corr_ESBL_c <- cor(int_ESBL_V1V5, treat$ISOLATES_POSITIVE_V1, use = "complete.obs")
Corr_ESBL_c
int_ESBL_V1V5

no_na_V1 <- V1 %>% 
  ifelse("N/A", remove(ESBL_V1), ESBL_V1)
no_na_V1 <- V1 %>% 
  ifelse("N/A", remove(ESBL_V1), ESBL_V1)
no_na_V1 <- V1 %>% 
  remove(N/A)


