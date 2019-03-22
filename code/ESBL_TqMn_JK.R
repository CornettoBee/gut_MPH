# Load libraries
library(tidyverse)
install.packages("readxl")
library(readxl)

# Load data.  This data is from TrEAT TD trial and contains observations from 366
# Where do the 366 come from? TrEAT only randomized 339. 
# The numbers given are CT/Cq values i.e. # of cycles of PCR necessary in order to 
# reach threshold of detection. We deem 35 or less to indicate the genetic material 
# present (this value is low enough that it )
treat <- read_excel("data/raw/TrEAT_Merge_ESBL_2018.09.13_v2.XLSX")

# Get a conception of the data set:
ls()  # contents of current work space
class(TrEAT_Merge_ESBL_2018.09.13_v2) # dataframe
dim(TrEAT_Merge_ESBL_2018.09.13_v2)  
object.size(TrEAT_Merge_ESBL_2018.09.13_v2)  # almost one MB
summary(TrEAT_Merge_ESBL_2018.09.13_v2)  # not very helpful

# Select relevant data
S_TrEAT_Merge_ESBL <- select(TrEAT_Merge_ESBL_2018.09.13_v2, STUDY_ID_TRUNC, AGE,   
                             SEX, ends_with("STOOL"), ends_with("CARD"), ends_with("V1"), ends_with("V5")) 

# Check for missing data in the Aeromonas_STOOL column:
is.na(Slim_TrEAT_Merge_ESBL[, "Aeromonas_STOOL"])  # This result gives all FALSE which
# doesn't make sense based on what you see when looking at dataframe, namely "Undetermined"
# and just blank. R must see them as... factors! (tested this below)
is.character(Slim_TrEAT_Merge_ESBL[, "Aeromonas_STOOL"]) # FALSE
is.factor(Slim_TrEAT_Merge_ESBL[, "Aeromonas_STOOL"])    # TRUE

# Let's find the average TLUS (time to last unformed stool i.e. time to cure) indexed
# /grouped by site:
> tapply(TrEAT_Merge_ESBL_2018.09.13_v2$TLUS, TrEAT_Merge_ESBL_2018.09.13_v2$SITE, mean)
# This command does work but sites 48, 61 and 92 all render NA. There is at least one
# NA in the data for site 48 so perhaps we just need to tell R how to handle NAs. 
# Let's make a new dataframe with only SITE and TLUS then remove NAs from it. 
SITE_TLUS <- select(TrEAT_Merge_ESBL_2018.09.13_v2, SITE, TLUS)
class(SITE_TLUS)
str(SITE_TLUS)
tapply(SITE_TLUS$TLUS, SITE_TLUS$SITE, mean) # same problem, NAs nullify any attempt
# to get the mean.  Tried using na.omit, na.exclude, na.pass but none of those work

no_na <- !is.na(SITE_TLUS) # Just gives whether true or false NA

mean(no_na) # But taking mean will give us proportion of TRUE since TRUE and FALSE are
# coded as 1 and 0, respectively. This way we can see how many NAs there are and a better
# idea of what effect removal will have on the data. (here: 99% are TRUE i.e. few NAs)
# We will see below exactly what number of NAs are removed.

# To remove the few missing values
pure_SITE_TLUS <- remove_missing(SITE_TLUS)  # (4 rows removed)

# Then apply the mean to each TLUS as indexed/grouped by SITE:
tapply(pure_SITE_TLUS$TLUS, pure_SITE_TLUS$SITE, mean)  # we see that SITE 76 with 
# much higher mean time (57hrs) than the other 4 sites (14-17hrs). 

######

# Let's narrow the data we want to observe:
Slim_TrEAT_Merge_ESBL <- select(TrEAT_Merge_ESBL_2018.09.13_v2, STUDY_ID_TRUNC, AGE, 
                                SEX, ends_with("V1"), ends_with("V5")) 

Corr_ESBL <- Slim_TrEAT_Merge_ESBL %>%
  cor(ESBL_V1, ISOLATES_POSITIVE_V1)

pure_Slim_TrEAT_Merge_ESBL <- remove_missing(Slim_TrEAT_Merge_ESBL)  # I am not certain 
# which values R is assuming are missing so this is dangerous move, besides the fact that
# I've removed over half the rows (192/366). 
# If we only look at ESBL_V1, ISOLATES_POSITIVE_V1
V1 <- select(TrEAT_Merge_ESBL_2018.09.13_v2, ESBL_V1, ISOLATES_POSITIVE_V1)

remove_missing(V1)  # 79 rows w/ missing values removed

# To find the correlation or to perform the McNemar test b/w ESBL_V1 and 
# ISOLATES_POSITIVE_V1 they must be numeric. 
# Change ESBL_V1 into numeric. 
num_ESBL_V1 <- as.numeric(TrEAT_Merge_ESBL_2018.09.13_v2$ESBL_V1) # Now NS coded as 3, 
# Negative coded as 2, positive coded 4 and N/A as 1. 
# However, 
# the as.numeric function has now converted my dataframe into a numeric vector.
# We want Negative = 0 Positive = 1 and N/A as NA, NS as NA.
ESBL_V1_unc <- unclass(num_ESBL_V1)
binary_ESBL_V1 <- ifelse(num_ESBL_V1$ESBL_V1 == 3, NA, 
                         ifelse(ESBL_V1 = 1, NA, 
                                ifelse(ESBL_V1 = 2, 0, 
                                       ifelse(ESBL_V1 = 4, 1, NA)))) # does not work

binary_ESBL_V1 <- ifelse(num_ESBL_V1[ESBL_V1] == 3, NA, 
                         ifelse(ESBL_V1 = 1, NA, 
                                ifelse(ESBL_V1 = 2, 0, 
                                       ifelse(ESBL_V1 = 4, 1, NA))))# does not work
######### 

num_ESBL_V1 <- Slim_TrEAT_Merge_ESBL %>%
  mutate(numbers_V1 = ifelse(ESBL_V1 == "NS", 0, 
                             ifelse(ESBL_V1 == "Negative", 0, 
                                    ifelse(ESBL_V1 == "Positive", 1, 
                                           ifelse(ESBL_V1 == "N/A", NA, NA)))))  # Success!

test_mcN_ESBL <- select(num_ESBL_V1, numbers_V1, ISOLATES_POSITIVE_V1)

no_ESBL_V1 <- Slim_TrEAT_Merge_ESBL$ESBL_V1 %>%
  count()

num_ESBL_V1$numbers_V1 %>%
  sum(na.rm = TRUE)   # gives 12

d <- num_ESBL_V1$numbers_V1 # then take sum:

sum(d, na.rm = TRUE)  # gives 12 ; so there are 12 positives
# Using phenotypic testing we find that there are 12/366 observations that were positive
# among non-pathogenic E. coli on V1 (visit 1) for ESBL. 
# However, as noted below, there are 55 NAs among this data. Do these NAs indicate not 
# tested, not enough sample, or nothing grew? 

num_ESBL_V1$ISOLATES_POSITIVE_V1 %>%
  sum(na.rm = TRUE)   # gives 20 because sometimes a positive results in multiple positive
# isolates

e <- num_ESBL_V1$ISOLATES_POSITIVE_V1
sum(e, na.rm = TRUE)  # gives 20 

# How many NAs are in the numbers_V1 vector/variable?
na_V1 <- is.na(num_ESBL_V1$numbers_V1)
table(na_V1)  # So there are 55 TRUE which means 55 NAs

# Select the variables of interest:
CRD_STL_TqMn <- select(TrEAT_Merge_ESBL_2018.09.13_v2, 
                       starts_with("Bacterial"), starts_with("AGE"), starts_with("CTX"), starts_with("KPC"), starts_with("NDM"), starts_with("SHV"), starts_with("TEM"), starts_with("CMY"), starts_with("STUDY"))

# Problem: all these values in this dataframe are coded as factors.  We want them to be
# numeric

f <- CRD_STL_TqMn$Bacterial_16s_STOOL
as.numeric(levels(f))[as.integer(f)]

as.numeric(levels(f))[f]
g <- CRD_STL_TqMn$CTX_STOOL

as.numeric(levels(g))[as.integer(f)]
str(g)
g > 30

# Some failed attempts to convert a factor vector to numeric using case_when:
num_ESBL_V1_c <- Slim_TrEAT_Merge_ESBL %>%
  case_when(
    ESBL_V1 == "NS" ~ 0, 
    ESBL_V1 == "Negative" ~ 0, 
    ESBL_V1 == "Positive" ~ 1, 
    ESBL_V1 == "N/A" ~ NA_integer_
  )

levels(Slim_TrEAT_Merge_ESBL$ESBL_V1)

num_ESBL_V1_d <- Slim_TrEAT_Merge_ESBL %>%
  case_when(
    ESBL_V1 == level 3 ~ 0, 
    ESBL_V1 == 2 ~ 0, 
    ESBL_V1 == 4 ~ 1, 
    ESBL_V1 == 1 ~ NA_integer_
  )

num_ESBL_V1_e <- ESBL_V1_unc %>%
  case_when(
    ESBL_V1 == 3 ~ 0, 
    ESBL_V1 == 2 ~ 0, 
    ESBL_V1 == 4 ~ 1, 
    ESBL_V1 == 1 ~ NA_integer_
  )

num_ESBL_V1_e2 <- ESBL_V1_unc %>%
  case_when(
    ESBL_V1 = 3 ~ 0, 
    ESBL_V1 = 2 ~ 0, 
    ESBL_V1 = 4 ~ 1, 
    ESBL_V1 = 1 ~ NA_integer_
  )

num_ESBL_V1_c2 <- Slim_TrEAT_Merge_ESBL %>%
  case_when(
    .$ESBL_V1 == "NS" ~ 0, 
    .$ESBL_V1 == "Negative" ~ 0, 
    .$ESBL_V1 == "Positive" ~ 1, 
    .$ESBL_V1 == "N/A" ~ NA_integer_,
    TRUE ~ as.numeric(as.character(.$ESBL_V1))
  )

new_ESBL_V1 <- Slim_TrEAT_Merge_ESBL %>%
  transmute()



Corr_ESBL_c <- cor(num_ESBL_V1, TrEAT_Merge_ESBL_2018.09.13_v2$ISOLATES_POSITIVE_V1, use = "complete.obs")
Corr_ESBL_c
num_ESBL_V1

no_na_V1 <- V1 %>% 
  ifelse("N/A", remove(ESBL_V1), ESBL_V1)
no_na_V1 <- V1 %>% 
  ifelse("N/A", remove(ESBL_V1), ESBL_V1)
no_na_V1 <- V1 %>% 
  remove(N/A)


