###############################################################################
# Project title: "Do egg hormones have fitness consequences in wild birds? 
# A systematic review and meta-analysis"

# Script from Alfredo Sanchez-Tojar - modified by LM - 28 Jan 2021

###############################################################################
# Description of script and Instructions
###############################################################################

# This script is to import the results of a systematic review
# on egg hormones and fitness.

###############################################################################

###############################################################################
# Packages needed
###############################################################################

# load pacakges
pacman::p_load(dplyr, revtools)

# cleaning up
rm(list=ls())


###############################################################################
# Functions needed
###############################################################################

# none


###############################################################################
# Importing reference data
###############################################################################

# The literature searches were conducted on the 27th January 2021
# Timespan = all years. Remaining settings = default.


################################
# WEB OF SCIENCE: general search

# importing the .bib files
wos <- read_bibliography("data/raw_data/systematic_review/02_search_results/Web of Science - 355 - 27 Jan 2021.bib")

# reducing fields to the minimum number of fields
# so that all databases have the same columns. Also, these fields
# are the important ones for the screening (see below).
reducing.fields.wos <- c("label","title","author","journal","issn","volume",
                         "number","pages","year","publisher","doi","abstract") #number to issue, 

wos.red <- wos[,reducing.fields.wos]


########################
# SCOPUS: general search

# importing the .bib files
scopus <- read_bibliography("data/raw_data/systematic_review/02_search_results/Scopus list - 1778 - 28 Jan 2021.bib")

# reducing fields to the minimum number of fields
# so that all databases have the same columns. Also, these fields
# are the important ones for the screening (see below).
reducing.fields.scopus <- c("label","title","author","journal","issn","volume",
                            "number","pages","year","publisher","doi","abstract") #number to issue,

scopus.red <- scopus[,reducing.fields.scopus]

###############################################################################
# Full reference data: before deduplication
###############################################################################

# building the full reference list before deduplication
full.ref.data <- rbind(wos.red,
                       scopus.red)

write.csv(full.ref.data,"data/raw_data/systematic_review/02_search_results/deduplic.csv",row.names=FALSE)
full.ref.data <- read.table("data/raw_data/systematic_review/02_search_results/deduplic.csv",
                           header=T,sep=",")


###############################################################################
# Full reference data: after deduplication
###############################################################################

# # searching duplicates using revtools: 
search.duplicated <- find_duplicates(data = full.ref.data,
                                     match_variable = "title",
                                      group_variable = NULL,
                                      match_function = "fuzzdist",
                                      method = "fuzz_m_ratio",
                                      remove_punctuation = T,
                                      threshold = 0.1) # 0.1 = 3578, 0.2 = 3558 I choose 0.2
# 
# # extracing duplicates
screening.ref.data <- extract_unique_references(full.ref.data,search.duplicated)
write.csv(screening.ref.data,"data/raw_data/systematic_review/02_search_results/search_unique_references_extracted.csv",row.names=FALSE)
screening.ref.data <- read.table("data/raw_data/systematic_review/02_search_results/search_unique_references_extracted.csv",
                                 header=T,sep=",")


###############################################################################
# Formatting data for RAYYAN QCRI
###############################################################################

# Choose only the fields needed for creating a .csv file importable by: https://rayyan.qcri.org

# Example of a valid .csv file. The fields are the following:
# key,title,authors,journal,issn,volume,issue,pages,year,publisher,url,abstract
names.rayyan <- c("key","title","authors","journal","issn","volume","issue",
                  "pages","day", "month", "year","publisher","pmc_id", "pubmed",
                  "url","abstract", "notes")
names.rayyan
names(screening.ref.data)

# Standardizing fields according to rayyan.example

# What's different between the two?
setdiff(names.rayyan,names(screening.ref.data))
setdiff(names(screening.ref.data), names.rayyan)


# Excluding variables that are not needed
#screening.ref.data$n_duplicates <- NULL
#screening.ref.data$origin <- NULL

#adding variables that are not included
screening.ref.data$pmc_id <- "" #simply adding an empty column since we did not extracted.
screening.ref.data$pubmed <- "" #simply adding an empty column since we did not extracted.
screening.ref.data$notes <- "" #simply adding an empty column since we did not extracted.
screening.ref.data$day <- "" #simply adding an empty column since we did not extracted.
screening.ref.data$month <- "" #simply adding an empty column since we did not extracted.

# what's different now?
setdiff(names.rayyan, names(screening.ref.data))
setdiff(names(screening.ref.data), names.rayyan)


screening.ref.data.rayyan <- plyr::rename(screening.ref.data, 
                                          c("label" = "key", "author" = "authors", 
                                            "number" = "issue", "doi" = "url"))
names(screening.ref.data.rayyan)

# Reorder
screening.ref.data.rayyan <- screening.ref.data.rayyan[, names.rayyan]


# finding authors with missing initial(s) as that causes an error when importing into rayyan
table(grepl(",  ", screening.ref.data.rayyan$authors, fixed = T))

for(i in 1:nrow(screening.ref.data.rayyan)){
  
  if(grepl(",  ", screening.ref.data.rayyan$authors[i], fixed = T)){
    
    print(i)
  }
  
}

# manual fixes (i might differ for the same data)
# none


###############################################################################
# Creating output
###############################################################################

write.csv(screening.ref.data.rayyan[order(screening.ref.data.rayyan$title),],
          "data/raw_data/systematic_review/02_search_results/search_unique_references_rayyan.csv",row.names=FALSE)
#remember to manually remove the quotes for the column names only in the .csv file


# saving versions used for reproducibility purposes
sink("data/raw_data/systematic_review/02_search_results/deduplicating_Rpackages_session.txt")
sessionInfo()
sink()
