# Load libraries
library(tidyverse)
library(readxl)

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
# NS -> 0
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
sum(is.na(int_ESBL_V1V5$int_V1))


int_ESBL_V1V5 %>%
  select(int_V1) %>% 
  group_by(int_V1) %>% 
  tally()

# Using phenotypic testing we find that there are 12/366 observations that were positive
# among non-pathogenic E. coli on V1 (visit 1) for ESBL. 
# However, as noted below, there are 55 NAs among this data. Do these NAs indicate not 
# tested, not enough sample, or nothing grew? 


# How many NAs are in the int_V1 vector/variable?
na_V1 <- is.na(int_ESBL_V1V5$int_V1)
table(na_V1)  # So there are 55 TRUE which means 55 NAs



########################
#### TaqMan data #######
########################
# Select the variables of interest:
card_stool_tqmn_V1 <- select(treat, 
                       starts_with("Bacterial"), starts_with("AGE"), starts_with("CTX"), starts_with("KPC"), 
                       starts_with("NDM"), starts_with("SHV"), starts_with("TEM"), starts_with("CMY"), starts_with("STUDY"))



##### Problem: #####
# All these values in this dataframe are coded as factors.  We want them to be
# numeric or integer type. 
# First, what values do these variables contain in terms of our analysis? There are 3 types:
# number(Cq value), "Undetermined" and NA. We interpret Undetermined and NA as meaning 
# there was no bacteria detected using the TaqMan real-time PCR containing the gene 
# coding for the type of ESBL resistance specified. 
# Let's make these character vectors into integer. 
# Since the Cq values in this data are unique numbers it's not feasible to use an ifelse 
# statement for each one of them when converting them into numeric type data. Instead, 
# let's designate each Undetermined or NA as NA while leaving the rest as is.  Then, we 
# can convert them to integer. 


# Use the scoped forms of mutate and transmutate such as mutate_if in order to change
# from one type of variable to another type (thanks: https://dplyr.tidyverse.org/reference/mutate_all.html)
tqmn_na_int <- card_stool_tqmn_V1 %>%
  mutate_at(vars(c(-STUDY_ID_TRUNC)), funs(as.double))  # note: introduces NAs by coercion, including all 

# Ryan's code for interpreting the taq values
tqmn_na_int <- card_stool_tqmn_V1 %>%
  # Change 'Undertermined' to max threshold (40)
  mutate_at(vars(-STUDY_ID_TRUNC), funs(ifelse(. == "Undetermined", 40, .))) %>%
  # Convert all columns to numeric (will coerce 'Indeterminate' to NA)
  mutate_at(vars(-STUDY_ID_TRUNC), as.numeric) %>%
  # Convert all values under 35 to 1 and over 35 to 0
  mutate_at(vars(-STUDY_ID_TRUNC), funs(ifelse(. < 35, 1, 0)))


# Check:
str(tqmn_na_int)
tqmn_na_int$CTX_STOOL[10:50]


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
}    # Testing with: df_test <- probe_test3("CTX_STOOL") yields df with probe_pos column but all 
# are NAs. 

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

##### Missing Data Analysis #####


# Check for missing data in the ESBL_V1 column:
is.na(Slim_treat[, "ESBL_V1"])  # This result gives all FALSE which
# doesn't make sense based on what you see when looking at dataframe, namely "Undetermined"
# and just blank. R must see them as... characters. (tested this below)
str(Slim_treat$ESBL_V1)   # PROBLEM:  R sees 
str(treat)
pure_Slim_treat <- remove_missing(Slim_treat)  # I am not certain 
# which values R is assuming are missing so this is dangerous move, besides the fact that
# I've removed over half the rows (192/366). 
# If we only look at ESBL_V1, ISOLATES_POSITIVE_V1
V1 <- select(treat, ESBL_V1, ISOLATES_POSITIVE_V1)

remove_missing(V1)  # 79 rows w/ missing values removed

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


