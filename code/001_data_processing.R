################################################################################
# Authors: 
#
# Lucia Mentesana (luciamentesana@gmail.com) wrote most of the original code
#
# Alfredo Sanchez-Tojar (alfredo.tojar@gmail.com) revised & adjusted the code, 
# wrote some original code (nostly for the analytical part) and prepared the 
# code for publication

# Script first created in Feb 2023

################################################################################
# Description of script and Instructions
################################################################################

# Script used to process the meta-analytic data for the following meta-analysis:

# Do egg hormones have fitness consequences in wild birds? A systematic review 
# and meta-analysis

# by Lucia Mentesana, Michaela Hau, Pietro B. D'Amelio, Nicolas M. Adreani, and
# Alfredo Sanchez-Tojar. 

# preprint available at: https://doi.org/10.1101/2024.10.29.620852


# This script:

# 1) Calculates the final sample sizes for rows that reported statistical tests
# and those that reported means, LSM (least squared means) and absolute values

# 2) Estimates effect sizes of the relationship between yolk hormones and 
# fitness traits, and their sampling variances

# 3) Prepares the data frame for analysis (i.e., checking all variables within 
# each factor)


# NOTE: corticosterone is the main glucocorticoid present in birds. Throughout
# the code we use glucocorticoids and corticosterone interchangeably. 


################################################################################
# Packages needed
################################################################################

# install.packages("pacman") #if not already installed
# pacman::p_load(metafor, plyr, stringr, rotl, ape, orchaRd, dplyr, 
#                devtools, tidyverse, patchwork, R.rsp, emmeans, ggplot2, maps, 
#                phytools, flextable, wesanderson)
pacman::p_load(metafor,plyr,stringr)

# cleaning up
rm(list = ls())


################################################################################
# Functions:
################################################################################

# none

################################################################################
# Importing the raw data 
################################################################################

meta.final <- read.csv("data/raw_data/meta-analysis_egg_hormones_raw_data.csv")


################################################################################
# 0) Organizing variables
################################################################################

# exploring and transforming variables
names(meta.final)

# As factors:
meta.final$Include <- as.factor(meta.final$Include)
meta.final$StudyID <- as.factor(meta.final$StudyID)
meta.final$GroupID <- as.factor(meta.final$GroupID)
meta.final$LaboratoryID <- as.factor(meta.final$LaboratoryID)
meta.final$PopulationID <- as.factor(meta.final$PopulationID)
meta.final$EffectID <- as.factor(meta.final$EffectID)
meta.final$Year <- as.factor(meta.final$Year)
meta.final$Lab_PI <- as.factor(meta.final$Lab_PI)
meta.final$Species <- as.factor(meta.final$Species)
meta.final$Developmental_mode <- as.factor(meta.final$Developmental_mode)
meta.final$Focal_clutch_number <- as.factor(meta.final$Focal_clutch_number)
meta.final$Study_area <- as.factor(meta.final$Study_area)
meta.final$Year_collection <- as.factor(meta.final$Year_collection)
meta.final$Study_type <- as.factor(meta.final$Study_type)
meta.final$Egg_n_collected <- as.factor(meta.final$Egg_n_collected)
meta.final$Time_collection <- as.factor(meta.final$Time_collection)
meta.final$Egg_sampl_method <- as.factor(meta.final$Egg_sampl_method)
meta.final$Egg_store <- as.factor(meta.final$Egg_store)
meta.final$Egg_homog <- as.factor(meta.final$Egg_homog)
meta.final$Egg_homog_name <- as.factor(meta.final$Egg_homog_name) #only contains NA's because no information was provided in the original articles
meta.final$Site_horm_measured <- as.factor(meta.final$Site_horm_measured)
meta.final$Horm_measured <- as.factor(meta.final$Horm_measured)
meta.final$Horm_extracted_colum <- as.factor(meta.final$Horm_extracted_colum)
meta.final$Method_horm_meas <- as.factor(meta.final$Method_horm_meas)
meta.final$AntibodyKit <- as.factor(meta.final$AntibodyKit)
meta.final$Horm_units <- as.factor(meta.final$Horm_units)
meta.final$Horm_Var_Type <- as.factor(meta.final$Horm_Var_Type)
meta.final$Exp_mother_egg <- as.factor(meta.final$Exp_mother_egg)
meta.final$Method_mother_exp <- as.factor(meta.final$Method_mother_exp)
meta.final$Type_mother_exp <- as.factor(meta.final$Type_mother_exp)
meta.final$Horm_mother_exp <- as.factor(meta.final$Horm_mother_exp)
meta.final$Horm_mother_units <- as.factor(meta.final$Horm_mother_units)
meta.final$Horm_mother_lab <- as.factor(meta.final$Horm_mother_lab)
meta.final$Horm_mother_cateogrical_dose <- as.factor(meta.final$Horm_mother_cateogrical_dose)
meta.final$Time_mother_exp <- as.factor(meta.final$Time_mother_exp)
meta.final$Egg_order_exp <- as.factor(meta.final$Egg_order_exp)
meta.final$Time_egg_exp <- as.factor(meta.final$Time_egg_exp)
meta.final$Method_egg_exp <- as.factor(meta.final$Method_egg_exp)
meta.final$Site_egg_exp <- as.factor(meta.final$Site_egg_exp)
meta.final$Horm_egg_exp <- as.factor(meta.final$Horm_egg_exp)
meta.final$Horm_egg_units <- as.factor(meta.final$Horm_egg_units)
meta.final$Horm_egg_lab <- as.factor(meta.final$Horm_egg_lab)
meta.final$Horm_egg_cateogrical_dose <- as.factor(meta.final$Horm_egg_cateogrical_dose)
meta.final$Solvent_egg_used_exp <- as.factor(meta.final$Solvent_egg_used_exp)
meta.final$Fitness_trait <- as.factor(meta.final$Fitness_trait)
meta.final$Fitness_mother <- as.factor(meta.final$Fitness_mother)
meta.final$Fitness_offs <- as.factor(meta.final$Fitness_offs)
meta.final$Off_life_stage <- as.factor(meta.final$Off_life_stage)
meta.final$Off_age <- as.factor(meta.final$Off_age)
meta.final$Off_sex <- as.factor(meta.final$Off_sex)
meta.final$Blind_data_recording <- as.factor(meta.final$Blind_data_recording)
meta.final$Partial_or_selective_data_rep <- as.factor(meta.final$Partial_or_selective_data_rep)
meta.final$TestStat_N_level <- as.factor(meta.final$TestStat_N_level)
meta.final$TestStat_Type <- as.factor(meta.final$TestStat_Type)
meta.final$TestDirection <- as.factor(meta.final$TestDirection)
meta.final$Exp_N_level <- as.factor(meta.final$Exp_N_level)
meta.final$Shared_control <- as.factor(meta.final$Shared_control)
meta.final$Shared_control_ID <- as.factor(meta.final$Shared_control_ID)
meta.final$Var_type <- as.factor(meta.final$Var_type)
meta.final$ExtractionID <- as.factor(meta.final$ExtractionID)
meta.final$Extraction_check_observer <- as.factor(meta.final$Extraction_check_observer)
meta.final$Least_square_means <- as.factor(meta.final$Least_square_means)
meta.final$Cross_fostering <- as.factor(meta.final$Cross_fostering)
meta.final$Incubation_prior <- as.factor(meta.final$Incubation_prior)
meta.final$Effect_size_type <- as.factor(meta.final$Effect_size_type)


# As numeric:
meta.final$Mean_clutch_size <- as.numeric(meta.final$Mean_clutch_size)
meta.final$SD_clutch_size <- as.numeric(meta.final$SD_clutch_size)
meta.final$SE_clutch_size <- as.numeric(meta.final$SE_clutch_size)
meta.final$Max_clutch_number <- as.numeric(meta.final$Max_clutch_number)
meta.final$Lat_study_area <- as.numeric(meta.final$Lat_study_area)
meta.final$Long_study_area <- as.numeric(meta.final$Long_study_area)
meta.final$Elev_study_area <- as.numeric(meta.final$Elev_study_area)
meta.final$Horm_mother_applied_dose <- as.numeric(meta.final$Horm_mother_applied_dose)
meta.final$Horm_mother_mean_nat <- as.numeric(meta.final$Horm_mother_mean_nat)
meta.final$Horm_mother_dose_exp <- as.numeric(meta.final$Horm_mother_dose_exp)
meta.final$Horm_mother_SD <- as.numeric(meta.final$Horm_mother_SD)
meta.final$Percent_egg_exp <- as.numeric(meta.final$Percent_egg_exp)
meta.final$Horm_egg_dose_exp <- as.numeric(meta.final$Horm_egg_dose_exp)
meta.final$Horm_egg_SD <- as.numeric(meta.final$Horm_egg_SD) #only contains NA's because no information was provided in the original articles
meta.final$Off_average_nesting_time <- as.numeric(meta.final$Off_average_nesting_time)
meta.final$Off_relative_age <- as.numeric(meta.final$Off_relative_age)
meta.final$TestStat_N_off <- as.numeric(meta.final$TestStat_N_off) 
meta.final$TestStat_NumberOfPredictors <- as.numeric(meta.final$TestStat_NumberOfPredictors) 
meta.final$TestStat_Value <- as.numeric(meta.final$TestStat_Value) 
meta.final$TestStat_df <- as.numeric(meta.final$TestStat_df) 
meta.final$Control_Mean <- as.numeric(meta.final$Control_Mean) 
meta.final$Control_Var <- as.numeric(meta.final$Control_Var) 
meta.final$Control_N_clutch <- as.numeric(meta.final$Control_N_clutch) 
meta.final$Control_N_off <- as.numeric(meta.final$Control_N_off) 
meta.final$Exp_Mean <- as.numeric(meta.final$Exp_Mean) 
meta.final$Exp_Var <- as.numeric(meta.final$Exp_Var) 
meta.final$Exp_N_clutch <- as.numeric(meta.final$Exp_N_clutch) 
meta.final$Exp_N_off <- as.numeric(meta.final$Exp_N_off) 


# Variables  that should be numeric but there is text on them (for example, 
# when researchers manipulated two hormones simultaneously and so the amount of 
# hormone experimentally added to the bird was reported separately for each 
# hormone)

# Leaving them as characters for the time being.
meta.final$Horm_mean <- as.character(meta.final$Horm_mean)
meta.final$Horm_Var <- as.character(meta.final$Horm_Var)
meta.final$Egg_n_exp <- as.character(meta.final$Egg_n_exp)
meta.final$Horm_egg_applied_dose <- as.character(meta.final$Horm_egg_applied_dose)
meta.final$Horm_egg_mean_nat <- as.character(meta.final$Horm_egg_mean_nat)
meta.final$TestStat_N_clutch <- as.character(meta.final$TestStat_N_clutch)
meta.final$TestStat_p <- as.character(meta.final$TestStat_p) 

names(meta.final)
summary(meta.final)

################################################################################
# 1) Obtaining final sample sizes
################################################################################

################################################################################
# A) For papers reporting test statistics and correlations

meta.final_stats <- droplevels(subset(meta.final, 
                                      meta.final$Effect_size_type != "Mean" & 
                                        meta.final$Effect_size_type != "Absolute" & 
                                        meta.final$Effect_size_type != "LSM"))

table(meta.final_stats$Effect_size_type)

# exploring TestDirection
table(meta.final_stats$TestDirection)
meta.final_stats[,c("StudyID","Fitness_trait",
                    "TestStat_N_level","TestStat_N_clutch","TestStat_N_off",
                    "TestStat_Type","TestStat_Value","TestStat_df",
                    "TestDirection")]

table(meta.final_stats$TestStat_Type)
table(meta.final_stats$TestStat_Type,meta.final_stats$TestDirection)


################################################################################
# clutch level
meta.final_stats_clutch <- droplevels(subset(meta.final_stats, 
                                             meta.final_stats$TestStat_N_level == "clutch"))

# Please, be aware that although there are some NA's in the variable
# meta.final_stats$TestStat_N_level, those corresponds to data points that we
# anyway cannot include in the final dataset, and therefore, although excluded
# by this type of sub-setting, this will not affect our results

meta.final_stats_clutch$stats_n_total <- meta.final_stats_clutch$TestStat_N_clutch


# There are 2 papers for which we do not have the number of clutches but we have 
# the number of offspring. We will therefore use the number of offspring as the
# sample size
nrow_na <- which(is.na(meta.final_stats_clutch$stats_n_total)) 
meta.final_stats_clutch$stats_n_total[nrow_na] <- meta.final_stats_clutch$TestStat_N_off[nrow_na]


# There are other 12 studies for which we do not have the number of clutches nor
# the number of offspring, but we do have the number of degrees of freedom. We 
# originally planned to get the N from the degrees of freedom as N = DF + number
# of predictors, even for those df that have decimal points (i.e. those from 
# mixed models), which means that those sample sizes will not be integers. Since 
# we do not have number of predictors for some effect sizes, we finally decided 
# to use N = DF + 1, which is a conservative approach

# create a new variable for the calculation of N based on df
meta.final_stats_clutch$stats_n_df <- (meta.final_stats_clutch$TestStat_df) + 1

nrow_unclear <- which((meta.final_stats_clutch$stats_n_total == "unclear")) 
meta.final_stats_clutch$stats_n_total[nrow_unclear] <- meta.final_stats_clutch$stats_n_df[nrow_unclear]

# names(meta.final_stats_clutch)
nrow(meta.final_stats_clutch)

# keep the final data base that has complete information: test and N
meta.final_stats_clutch <- droplevels(subset(meta.final_stats_clutch, 
                                             meta.final_stats_clutch$stats_n_total != "NA"))
nrow(meta.final_stats_clutch)

# # quick exploration
# summary(meta.final_stats_clutch[,c("Fitness_trait","Fitness_mother",
#                                    "Fitness_offs","TestStat_N_off",
#                                    "TestStat_N_clutch","TestStat_N_level",
#                                    "TestStat_df","stats_n_total",
#                                    "Partial_or_selective_data_rep")])


# create three more columns that will not contain information, but are needed 
# for the rows that have information on mean control and experimental groups
# and one study that reported raw values

meta.final_stats_clutch$control_n_total <- "NA"
meta.final_stats_clutch$exp_n_total <- "NA"
meta.final_stats_clutch$raw_n <- "NA"
meta.final_stats_clutch$datasubset <- "stats_clutch"

# names(meta.final_stats_clutch)


################################################################################
# offspring level
meta.final_stats_offspring <- droplevels(subset(meta.final_stats, 
                                                meta.final_stats$TestStat_N_level == "offspring"))

meta.final_stats_offspring$stats_n_total <- meta.final_stats_offspring$TestStat_N_off


# There are 3 papers for which we do not have the number of offspring but we 
# have the number of clutches. We will therefore the number of clutches as the
# sample size
nrow_na_off <- which(is.na(meta.final_stats_offspring$stats_n_total)) 
meta.final_stats_offspring$stats_n_total[nrow_na_off] <- meta.final_stats_offspring$TestStat_N_clutch[nrow_na_off]


# There are other 6 studies for which we do not have the number of clutches or 
# the number of offspring, but we do have the degrees of freedom. We will 
# therefore get the N from the degrees of freedom N = DF + number of predictors
# For those values that have df from mixed models (and have decimals) we will 
# round up the values. Since we do not have number of predictors for many papers,
# we will consider it to be 1 to avoid excluding these data points

# create a new variable for the calculation of N based on df
meta.final_stats_offspring$stats_n_df <- (meta.final_stats_offspring$TestStat_df) + 1

nrow_na_df <- which(is.na(meta.final_stats_offspring$stats_n_total)) 
# this works because the approach above (nrow_na_off) did not fill in any of the NA's we are targetting here
meta.final_stats_offspring$stats_n_total[nrow_na_df] <- meta.final_stats_offspring$stats_n_df[nrow_na_df]


# There are 2 papers for which final sample size is 2. This occurred because df
# were reported for Chi-square tests, where df are calculated from the sample 
# size but from the number of groups that have been compared (i.e. n groups - 1)
# Since this cannot be used for effect size calculations, we are excluding 
# them from the final db
meta.final_stats_offspring <- droplevels(subset(meta.final_stats_offspring, 
                                                meta.final_stats_offspring$stats_n_total != 2))

meta.final_stats_offspring$stats_n_total


nrow(meta.final_stats_offspring)

# quick exploration
# summary(meta.final_stats_offspring[,c("Fitness_trait","Fitness_mother",
#                                    "Fitness_offs","TestStat_N_off",
#                                    "TestStat_N_clutch","TestStat_N_level",
#                                    "TestStat_df","stats_n_total",
#                                    "Partial_or_selective_data_rep")])

#names(meta.final_stats_clutch)

# Adding the corresponding columns that will be needed to merge all 5 data bases
# after effect sizes are calculated from the different sources
meta.final_stats_offspring$control_n_total <- "NA"
meta.final_stats_offspring$exp_n_total <- "NA"
meta.final_stats_offspring$raw_n <- "NA"
meta.final_stats_offspring$datasubset <- "stats_offspring"


################################################################################
# B) For papers reporting means, absolute values and LSM from control and 
# experimental (exp) groups

meta.final_means <- droplevels(subset(meta.final, meta.final$Effect_size_type == "Mean" |
                                        meta.final$Effect_size_type == "Absolute" | 
                                        meta.final$Effect_size_type == "LSM"))
table(meta.final_means$Effect_size_type)


################################################################################
# clutch level
meta.final_means_clutch <- droplevels(subset(meta.final_means, 
                                             meta.final_means$Exp_N_level == "clutch"))


################################################################################
# offspring level
meta.final_means_offspring <- droplevels(subset(meta.final_means, 
                                                meta.final_means$Exp_N_level == "offspring"))
nrow(meta.final_means_offspring)

# first create a column that contains information with no value but has 
# the same name as the one used for studies that reported test statistics
meta.final_means_clutch$stats_n_total <- "NA"
meta.final_means_offspring$stats_n_total <- "NA"

# add the columns needed to have all data bases with the same order
meta.final_means_clutch$stats_n_df <- "NA"

meta.final_means_offspring$stats_n_df <- "NA"


################################################################################
# Clutch level
meta.final_means_clutch$control_n_total <- meta.final_means_clutch$Control_N_clutch
meta.final_means_clutch$exp_n_total <- meta.final_means_clutch$Exp_N_clutch

# There are 2 rows for which there is no information at the clutch level, but 
# there is at the offspring level. I will use the info for offspring
nrow_na_clutch_c <- which(is.na(meta.final_means_clutch$control_n_total)) 
meta.final_means_clutch$control_n_total[nrow_na_clutch_c] <- meta.final_means_clutch$Control_N_off[nrow_na_clutch_c]


nrow_na_clutch_exp <- which(is.na(meta.final_means_clutch$exp_n_total)) 
meta.final_means_clutch$exp_n_total[nrow_na_clutch_exp] <- meta.final_means_clutch$Exp_N_off[nrow_na_clutch_exp]


# I now subset the db to get those rows that have complete info: means and N
meta.final_means_clutch <- droplevels(subset(meta.final_means_clutch, 
                                             meta.final_means_clutch$control_n_total != "NA" &
                                               meta.final_means_clutch$exp_n_total != "NA"))
nrow(meta.final_means_clutch)

# # quick exploration
# summary(meta.final_means_clutch[,c("Fitness_trait","Fitness_mother",
#                                    "Fitness_offs",
#                                    "Control_N_off","control_n_total","Control_N_clutch",
#                                    "Exp_N_off","exp_n_total","Exp_N_clutch",
#                                    "Partial_or_selective_data_rep")])


meta.final_means_clutch$raw_n <- "NA"
meta.final_means_clutch$datasubset <- "means_clutch"


################################################################################
# Offspring level
meta.final_means_offspring$control_n_total <- meta.final_means_offspring$Control_N_off
meta.final_means_offspring$exp_n_total <- meta.final_means_offspring$Exp_N_off

# There are 21 rows for which there is no info on sample size for offspring but 
# there is for clutch. I will therefore consider the N of clutches for those 
# rows
nrow_na_off_c <- which(is.na(meta.final_means_offspring$control_n_total)) 
meta.final_means_offspring$control_n_total[nrow_na_off_c] <- meta.final_means_offspring$Control_N_clutch[nrow_na_off_c]

nrow_na_off_exp <- which(is.na(meta.final_means_offspring$exp_n_total))
meta.final_means_offspring$exp_n_total[nrow_na_off_exp] <- meta.final_means_offspring$Exp_N_clutch[nrow_na_off_exp]


# I will subset the data base to have complete information in each row
meta.final_means_offspring <- droplevels(subset(meta.final_means_offspring, 
                                                meta.final_means_offspring$control_n_total != "NA" &
                                                  meta.final_means_offspring$exp_n_total != "NA"))

nrow(meta.final_means_offspring)

# # quick exploration
# summary(meta.final_means_offspring[,c("Fitness_trait","Fitness_mother",
#                                    "Fitness_offs",
#                                    "Control_N_off","control_n_total","Control_N_clutch",
#                                    "Exp_N_off","exp_n_total","Exp_N_clutch",
#                                    "Partial_or_selective_data_rep")])


meta.final_means_offspring$raw_n <- "NA"
meta.final_means_offspring$datasubset <- "means_offspring"

# names(meta.final_means_offspring)


### We also have 4 rows with information from 1 paper that contains the 
# relationship between hormone level and clutch size with absolute value. This 
# information appears as "Absolute" in the effect_size_type column and sample 
# sizes appear in the "test_stat" section. I will therefore subset this 
# information and add the 4 rows to the general data base.

meta.final_abs <- droplevels(subset(meta.final, 
                                    meta.final$StudyID == "rayyan-123815683"))

# Adding the corresponding columns that will be needed to merge all 5 data bases
# after effect sizes are calculated from the different sources
meta.final_abs$stats_n_total <- "NA"
meta.final_abs$stats_n_df <- "NA"
meta.final_abs$control_n_total <- "NA"
meta.final_abs$exp_n_total <- "NA"
meta.final_abs$datasubset <- "abs"


meta.final_abs$raw_n <- meta.final_abs$TestStat_N_clutch

#names(meta.final_stats_clutch)


################################################################################
# Combine all 5 data bases
################################################################################
meta.final.n <- rbind(meta.final_stats_clutch, meta.final_stats_offspring, 
                      meta.final_means_clutch, meta.final_means_offspring,
                      meta.final_abs)
# names(meta.final.n)
nrow(meta.final.n)

meta.final.n$stats_n_total <- as.numeric(meta.final.n$stats_n_total)
meta.final.n$control_n_total <-as.numeric(meta.final.n$control_n_total)
meta.final.n$exp_n_total <- as.numeric(meta.final.n$exp_n_total)

# exploring TestDirection
table(meta.final.n$Effect_size_type, meta.final.n$TestDirection)


################################################################################
# 2) Estimating effect sizes and their variances
################################################################################
# Functions to convert statistics into our effect size of interest and its 
# sampling variance

# To do this, we need to check which was the statistical analysis done to then 
# apply the different transformations. Values should be transformed into 
# correlations, either Pearson's or biserial correlations.

# create a new column where we will place the transformed values
meta.final.n$cor <- "NA"
meta.final.n$cor_var <- "NA"

# also create to other columns that are needed for the transformation of 
# means (to biserial cor) and absolute values
meta.final.n$sd_control <- "NA"
meta.final.n$sd_exp <- "NA"
meta.final.n$n_control_surv <- "NA"
meta.final.n$n_control_no_surv <- "NA"
meta.final.n$n_exp_surv <- "NA"
meta.final.n$n_exp_no_surv <- "NA"        

# Time to create different data bases based on the statistics used to then
# transform the data

################################################################################
# Convert chi-squared to Pearson's r

# The equation we are using is for a classical Chi-square test. As far as we
# know, there aren't transformations specifically designed for:
# Chi-square (mixed model), Chi-square (generalized linear model), Chi-square 
# (negative binomial mixed model), nor Wald-Chi square (linear model).
# We are aware of this but will nonetheless used the original transformation for
# these (k = 14 cases in our dataset)

meta.final.n_chi <- droplevels(subset(meta.final.n, 
                                      meta.final.n$Effect_size_type == "Chi-square" |
                                        meta.final.n$Effect_size_type == "Chi-square (generalized linear model)" |
                                        meta.final.n$Effect_size_type == "Chi-square (negative binomial mixed model)" |
                                        meta.final.n$Effect_size_type == "Chi-square (mixed model)" |
                                        meta.final.n$Effect_size_type == "Wald-Chi square (linear model)"))

meta.final.n_chi$stats_n_total <- as.numeric(meta.final.n_chi$stats_n_total)

# Function to convert chi-squared into Pearson's r (Nakagawa and Cuthill 2007, Table 2)
# chi.to.r <- function(chi, n){
#  r <- sqrt(chi/n)}
meta.final.n_chi$cor <- sqrt((meta.final.n_chi$TestStat_Value)/
                               (meta.final.n_chi$stats_n_total))


# Important: Chi-square tests are statistics that do not have a sign (i.e., cor
# values are always positive). Hence, we need to check directionality of effect. 
# We obtained this information from the text, the figures and/or after contacting
# the authors to obtain this information. A positive test direction indicates that 
# the fitness proxy increases after an increase of hormone concentration. 
# The control group indicates lower hormone concentration and the experimental 
# group indicates higher hormone concentration. 
# In those cases where we lack information on directionality, we wrote down 
# 'unclear'. That means that, even if we do have the effect size, since we 
# do not know its directionality, this information unfortunately cannot be used. 

# Exploring signs: 
table(meta.final.n_chi$Effect_size_type, meta.final.n_chi$TestDirection)

table(meta.final.n_chi$cor, meta.final.n_chi$TestDirection)
# 10 rows have a negative directionality that at the moment have a positive cor
# value. We need to modify the sign of these 10 rows. Important: one of these
# rows has a 'negative' fitness trait (hatching failure). 
# 5 rows with 'unclear' directionality that cannot be included in the final 
# analysis

# To adjust the sign of the effect size
meta.final.n_chi$cor <- ifelse(meta.final.n_chi$TestDirection == "negative",
                               meta.final.n_chi$cor * (-1),
                               meta.final.n_chi$cor)

# To exclude those effect sizes for which their sign is unclear:
meta.final.n_chi <- meta.final.n_chi[meta.final.n_chi$TestDirection != "unclear",]


# Function to obtain sampling variance of r (Borenstein et al. 2009, Equation 6.1)
# Vr <- function(r,N){
#  Vr <- ((1-r^2)^2)/(N-1)}

meta.final.n_chi$cor_var <- ((1 - (meta.final.n_chi$cor ^ 2)) ^ 2)/
  (meta.final.n_chi$stats_n_total - 1)

meta.final.n_chi$final_n <- meta.final.n_chi$stats_n_total
nrow(meta.final.n_chi)

################################################################################
# Convert F-test to Pearson's r (Lajeunesse, 2013; p201)

# The equation we are using is for a classical F-test one-way ANOVA. As far as 
# we know, there aren't transformations specifically designed for:
# F-test (generalized linear model), F-test (linear mixed model),
# F-test (MANOVA), nor F-test (repeated-measure models)
# We are aware of this but will nonetheless use the original transformation for
# these (k = 35 cases in our dataset)


meta.final.n_ftest <- droplevels(subset(meta.final.n, meta.final.n$Effect_size_type == "F-test (linear mixed model)" |
                                          meta.final.n$Effect_size_type == "F-test (generalized linear model)" |
                                          meta.final.n$Effect_size_type == "F-test (MANOVA)" |
                                          meta.final.n$Effect_size_type == "F-test (repeated-measure models)"))  

# F.to.r <- function(ftest,N){ (Lajeunesse, 2013; p201)
# r <- sqrt(ftest/(ftest+N-2))}
meta.final.n_ftest$cor <- sqrt((meta.final.n_ftest$TestStat_Value)/
                                 (meta.final.n_ftest$TestStat_Value + 
                                    meta.final.n_ftest$stats_n_total - 2))

# As mentioned above for the Chi-sqaure tests, F-test also do not have a sign 
# (i.e., cor # values are always positive). We therefore checked the data base
# and whenever necessary applied the same modifications are previously described.

# Exploring signs: 
table(meta.final.n_ftest$Effect_size_type, meta.final.n_ftest$TestDirection)

table(meta.final.n_ftest$cor, meta.final.n_ftest$TestDirection)
# 6 rows have a negative directionality that at the moment have a positive cor
# value. We need to modify the sign of these 6 rows.
# 19 rows with 'unclear' directionality that cannot be included in the final analysis

# To adjust the sign of the effect size
meta.final.n_ftest$cor <- ifelse(meta.final.n_ftest$TestDirection == "negative",
                                 meta.final.n_ftest$cor * (-1),
                                 meta.final.n_ftest$cor)

# To exclude those effect sizes for which their sign is unclear:
meta.final.n_ftest <- meta.final.n_ftest[meta.final.n_ftest$TestDirection != "unclear",]


# Calculating sampling variance for correlations
meta.final.n_ftest$cor_var <- ((1 - meta.final.n_ftest$cor ^ 2) ^ 2/
                                 (meta.final.n_ftest$stats_n_total - 1))

meta.final.n_ftest$final_n <- meta.final.n_ftest$stats_n_total
nrow(meta.final.n_ftest)


################################################################################
# Convert Linear mixed model to Pearson's r
# Unfortunately, we do not have the parameters required to do a transformation 
# from regression coefficients from generalized/linear (mixed) models. Thus, we
# cannot include the study that reported this estimate (n=1)


################################################################################
# Convert Spearman's rho to Pearson's r (Lajeunesse, 2013)

meta.final.n_spearman <- droplevels(subset(meta.final.n, 
                                           meta.final.n$Effect_size_type == "Spearman"))

# since n < 90, the formula is:
#spearman_to_pearson <- function(rho){ (Lajeunesse, 2013)
#  r <- 2*sin((spearman*pi/6)}

meta.final.n_spearman$cor <- 2 * sin(((meta.final.n_spearman$TestStat_Value * pi)/6))

meta.final.n_spearman$cor_var <- ((1 - meta.final.n_spearman$cor ^ 2) ^ 2/
                                    (meta.final.n_spearman$stats_n_total - 1))

# checking sign
table(meta.final.n_spearman$cor, meta.final.n_spearman$TestDirection)
# ok

meta.final.n_spearman$final_n <- meta.final.n_spearman$stats_n_total
nrow(meta.final.n_spearman)


################################################################################
# Convert t-values from a continuous predictor variable to Pearson's r
# (Nakagawa & Cuthill, 2007, p. 598, Equation 11)

# The equation we are using is for a classical independent t-test. As far as 
# we know, there aren't transformations specifically designed for:
# t-test (general linear model). We are aware of this but will nonetheless used 
# the original transformation for the only case in our dataset and later on run 
# a sensitivity analysis
meta.final.n_tt <- droplevels(subset(meta.final.n, 
                                     meta.final.n$Effect_size_type == "t-test (unpaired)" |
                                       meta.final.n$Effect_size_type == "t-test (general linear model)"))
meta.final.n_tt$Effect_size_type

# t.to.r <- function(t,N){ (Nakagawa & Cuthill, 2007, p. 598, Equation 11)
#  df <- N-2 # This assumes that there are no other fixed effects included in 
# the models from which t was obtained.
#  r <- t/sqrt(((t^2)+df))}

meta.final.n_tt$cor <- (meta.final.n_tt$TestStat_Value)/
  sqrt(((meta.final.n_tt$TestStat_Value ^ 2) + 
          (meta.final.n_tt$stats_n_total - 2)))

meta.final.n_tt$cor_var <- ((1 - meta.final.n_tt$cor ^ 2) ^ 2/
                              (meta.final.n_tt$stats_n_total - 1))

meta.final.n_tt$final_n <- meta.final.n_tt$stats_n_total
nrow(meta.final.n_tt)

# The t-test estimates do have a directionality. However, sometimes 
# authors used the experimental and not the control group as a reference. If
# this was the case, a positive cor value would indicate that the experiment has
# a negative effect on fitness traits, which is in the opposite direction as the
# criteria we are using in our analysis. Hence, we need to double check that indeed 
# the sign of the correlation reflects the criteria we are following. 
# We also need to check for potential 'unclear'. 

# Exploring signs: 
table(meta.final.n_tt$Effect_size_type, meta.final.n_tt$TestDirection)
# We do not have information on the directionality of the tests. All 6 rows have
# papers where this was not reported. This is an issue because, for example, in
# the paper rayyan-123815691 it seems that the reference was the experimental 
# group (see result of tarsus length on day 2, figure 2 - here authors report a
# negative t-test but the effect of the experiment on tarsus length is positive).
# From this paper we have 5 other t-tests where, contrary to what happen with 
# tarsus length, authors do not indicate the directionality of the t-test.

# Since we cannot confidently identify if the sign of the test indeed reflects the
# effect of the experiment, we took a conservative approach and decided NOT TO 
# INCLUDE THESE 6 ROWS INTO THE FINAL ANALYSIS.


################################################################################
# Convert Z-score to Pearson's r (Lajeunesse, 2013)

# The function we are using is for a classical Z-test. As far as 
# we know, there aren't transformations specifically designed for:
# z-test (mixed model). We are aware of this but will nonetheless used the
# original transformation for the two cases in our dataset

# For z-tests that are more complex (general model), we will run a normal z-test 
# and then run a sensitivity analysis.

meta.final.n_z <- droplevels(subset(meta.final.n, meta.final.n$Effect_size_type == "z-test" |
                                      meta.final.n$Effect_size_type == "z-test (mixed model)")) 

meta.final.n_z$Effect_size_type

# zscore_to_pearson <- function(z,n){ Lajeunesse, 2013)
#  r <- z/sqrt(n)}

meta.final.n_z$cor <- meta.final.n_z$TestStat_Value/
  sqrt(meta.final.n_z$stats_n_total)

meta.final.n_z$cor_var <- ((1 - meta.final.n_z$cor ^ 2) ^ 2/
                             (meta.final.n_z$stats_n_total - 1))

# z-test estimates do have directionality. however, as it happened above with
# t-tests, we need to confirm that the directionality of the estimate follows the
# criteria explained above (i.e., positive value means an increase in hormone
# concentration has a positive effect on fitness proxies)

# Exploring signs: 
table(meta.final.n_z$Effect_size_type, meta.final.n_z$TestDirection)
table(meta.final.n_z$cor, meta.final.n_z$TestDirection)
# Good. Cor sign and test directionality have matching information. 

meta.final.n_z$final_n <- meta.final.n_z$stats_n_total
nrow(meta.final.n_z)


################################################################################
# subset the data for Pearson, which is ready.

meta.final.n_pearson <- droplevels(subset(meta.final.n, 
                                          meta.final.n$Effect_size_type =="Pearson"))

meta.final.n_pearson$cor <- meta.final.n_pearson$TestStat_Value

meta.final.n_pearson$cor_var <- ((1 - meta.final.n_pearson$cor ^ 2) ^ 2/
                                   (meta.final.n_pearson$stats_n_total - 1))


# Exploring signs: 
table(meta.final.n_pearson$Effect_size_type, meta.final.n_pearson$TestDirection)
table(meta.final.n_pearson$cor, meta.final.n_pearson$TestDirection)

# Important: four values have opposite information. This is because these 4
# rows have a 'negative\ fitness trait (i.e., egg mortality). Test direction
# indicates the true direction, whereas the value indicates the reported estimate.
# For the rest, cor sign and test directionality information match.

# To adjust the sign of the effect size
meta.final.n_pearson_sign <- droplevels(subset(meta.final.n_pearson, 
                                               meta.final.n_pearson$Fitness_mother == "egg mortality" |
                                                 meta.final.n_pearson$Fitness_offs == "offspring mortality"))
meta.final.n_pearson_sign$cor <- (meta.final.n_pearson_sign$cor * (-1))

# now get the remaining database that has the correct signs, and then merge 
# both data bases.
meta.final.n_pearson_sign_ok <- droplevels(subset(meta.final.n_pearson, 
                                                  meta.final.n_pearson$Fitness_mother != "egg mortality" |
                                                    meta.final.n_pearson$Fitness_offs != "offspring mortality"))

meta.final.n_pearson_ok <- rbind(meta.final.n_pearson_sign, 
                                 meta.final.n_pearson_sign_ok)

table(meta.final.n_pearson_ok$cor, meta.final.n_pearson_ok$TestDirection)
# ok

meta.final.n_pearson_ok$final_n <- meta.final.n_pearson_ok$stats_n_total
nrow(meta.final.n_pearson_ok)


################################################################################
# Functions to convert means into Biserial correlations
################################################################################
# Calculating the biserial correlation coefficients. This is for those rows for
# which we have Mean data for each of the groups being compared. 

meta.final.n_mean <- droplevels(subset(meta.final.n, 
                                       meta.final.n$Effect_size_type == "Mean"))
nrow(meta.final.n_mean)


# The calculation requires the SD.
# We therefore have to calculate it.
# We subset those papers for which we already have SD.
meta.final.n_mean_sd <- droplevels(subset(meta.final.n_mean, 
                                          meta.final.n_mean$Var_type == "SD"))

# separate control and experimental sd
meta.final.n_mean_sd$sd_control <- meta.final.n_mean_sd$Control_Var
meta.final.n_mean_sd$sd_exp <- meta.final.n_mean_sd$Exp_Var


# do the same now but for those rows for which we have SE and SEM data.
meta.final.n_mean_se <- droplevels(subset(meta.final.n_mean, 
                                          meta.final.n_mean$Var_type == "SE" |
                                            meta.final.n_mean$Var_type == "SEM"))



################################################################################
# need to convert SE into SD.
# sd = se * sqrt(n)
meta.final.n_mean_se$control_n_total <- as.numeric(meta.final.n_mean_se$control_n_total)
meta.final.n_mean_se$sd_control <- (meta.final.n_mean_se$Control_Var) * 
  sqrt(meta.final.n_mean_se$control_n_total)

# Due to Taylor's law, a positive correlation between log(mean) and log(sd) is
# expected. We are exploring this as a way of finding potential outliers.
cor(log(meta.final.n_mean_se$Control_Mean),log(meta.final.n_mean_se$sd_control))
plot(log(meta.final.n_mean_se$Control_Mean), log(meta.final.n_mean_se$sd_control))
# Things look more or less as expected, but noisy


meta.final.n_mean_se$exp_n_total <- as.numeric(meta.final.n_mean_se$exp_n_total)
meta.final.n_mean_se$sd_exp <- (meta.final.n_mean_se$Exp_Var) * 
  sqrt(meta.final.n_mean_se$exp_n_total)

# Due to Taylor's law, a positive correlation between log(mean) and log(sd) is
# expected. We are exploring this as a way of finding potential outliers.
cor(log(meta.final.n_mean_se$Exp_Mean),log(meta.final.n_mean_se$sd_exp),
    use = "complete")
plot(log(meta.final.n_mean_se$Exp_Mean),log(meta.final.n_mean_se$sd_exp))
# Things look more or less as expected, but noisy


# There are 2 rows in the data base for which we do not have variance for the 
# experimental group because sample size was 1, and 1 for which we do not have 
# information on variance for the experimental group. Therefore, I subset the 
# data base to include only those with complete information
meta.final.n_mean_se <- droplevels(subset(meta.final.n_mean_se, 
                                          meta.final.n_mean_se$sd_exp != "NA"))
nrow(meta.final.n_mean_se)


################################################################################
# merge both databases
meta.final.n_mean_ok <- rbind(meta.final.n_mean_sd, meta.final.n_mean_se)
nrow(meta.final.n_mean_ok)


################################################################################
### To calculate biserial correlations we require means, sds, and ns for both
# the experimental and the control group. We have some papers for which the same
# control group is compared to several experimental groups, that is, we have
# some cases of shared-control non independence in our data set that we need to
# account for. We are going to use a simple but effective method to essentially
# reduced the weight that the control group would otherwise have in the analyses
# if we would ignore this non independence. Our method of choice in this case is
# to simply divide the sample size of the control group by the number of times
# that control group was compared to another group. That is, if a control group
# was compared to 2 experimental groups, we would divide the sample size of the 
# control group by 2.

# shared controls
meta.final.n_mean_ok_shared_c <- droplevels(subset(meta.final.n_mean_ok, 
                                                   meta.final.n_mean_ok$Shared_control == "yes"))
meta.final.n_mean_ok_shared_c$Shared_control_ID
str(meta.final.n_mean_ok_shared_c$control_n_total)

meta.final.n_mean_ok_shared_c$Shared_control_ID <- as.character(meta.final.n_mean_ok_shared_c$Shared_control_ID)

meta.final.n_mean_ok_shared_c <- ddply(meta.final.n_mean_ok_shared_c,
                                       c("Shared_control_ID"),
                                       transform, 
                                       control_n_total = control_n_total/
                                         (length(Shared_control_ID)))

# select the data base for non shared controls and then merge both data bases
meta.final.n_mean_ok_noshared_c <- droplevels(subset(meta.final.n_mean_ok, 
                                                     meta.final.n_mean_ok$Shared_control == "no")) 
nrow(meta.final.n_mean_ok_noshared_c)


# merge both data bases
meta.final.n_mean_ok <- rbind(meta.final.n_mean_ok_shared_c, 
                              meta.final.n_mean_ok_noshared_c)
summary(meta.final.n_mean_ok$control_n_total)


# Calculating the biserial correlation coefficients and their variance
db.biserial <- db.biserial <- escalc(measure = "RBIS",
                                     n2i = as.numeric(meta.final.n_mean_ok$control_n_total),
                                     n1i = as.numeric(meta.final.n_mean_ok$exp_n_total),
                                     m2i = as.numeric(meta.final.n_mean_ok$Control_Mean),
                                     m1i = as.numeric(meta.final.n_mean_ok$Exp_Mean),
                                     sd2i = as.numeric(meta.final.n_mean_ok$sd_control),
                                     sd1i = as.numeric(meta.final.n_mean_ok$sd_exp))

# yi is the correlation, vi the variance
# I will therefore add this information to the general db.
meta.final.n_mean_ok$cor <- db.biserial$yi
meta.final.n_mean_ok$cor_var <- db.biserial$vi


# Checking that the signs reflect the general criteria used: 
expected.sign <- ifelse(meta.final.n_mean_ok$Control_Mean > 
                          meta.final.n_mean_ok$Exp_Mean, -1, 1)

observed.sign <- ifelse(meta.final.n_mean_ok$cor > 0, 1, -1)

# Checking if they agree
table(expected.sign == observed.sign) 
# It look oks. The 6 cases that appear as FALSE correspond to those rows where
# control and experimental groups had the same mean values and cor is 0. 


meta.final.n_mean_ok$final_n <- as.numeric(meta.final.n_mean_ok$control_n_total) +
  as.numeric(meta.final.n_mean_ok$exp_n_total)

################################################################################
###  LSM
meta.final.n_lsm <- droplevels(subset(meta.final.n, 
                                      meta.final.n$Effect_size_type == "LSM"))
nrow(meta.final.n_lsm)
# In total, we have 17 rows with LSM. Unfortunately we cannot use because
# we are not aware of any function that allows to calculate effect sizes from LSM.


################################################################################
### Absolute value
# For transforming these values I use the "odds transformation" (Lajeunesse, 2013)
# For this I need to first transform the data base

meta.final.n_abs <- droplevels(subset(meta.final.n, 
                                      meta.final.n$Effect_size_type == "Absolute"))
nrow(meta.final.n_abs)


# The absolute values are generally traits associated with survival: hatching 
# success, fledgling success.
# We have four different types of reported information here:
# 1) absolute values (i.e., number of dead/alive individuals per control and 
# experimental treatment) 
# 2) Percentages
# 3) Percentages and also absolute number
# 4) Proportions

# This information appears in the column "variance_type_ comment"

# For the odds formula we need to obtain sample sizes of individuals that
# survived and those that didn't.
# We therefore need to create two new columns that contain this information.
# First we will create the new columns we need:
nrow(meta.final.n_abs)
  

# There are some rows (n = 4) for which we do not have information on number of
# clutches/offspring
# We will keep the rows for which we have complete information

meta.final.n_abs <- droplevels(subset(meta.final.n_abs, 
                                      meta.final.n_abs$control_n_total != "NA"))

nrow(meta.final.n_abs)


meta.final.n_abs$control_n_total <- as.numeric(meta.final.n_abs$control_n_total)
meta.final.n_abs$exp_n_total <- as.numeric(meta.final.n_abs$exp_n_total)


# I will work with different data bases depending on the type of information 
# that we have available.


################################################################################
# 1) Absolute value
meta.final.n_abs_value <- droplevels(subset(meta.final.n_abs, 
                                            meta.final.n_abs$Var_type_comment == "Absolute value; no variance"))
nrow(meta.final.n_abs_value)


meta.final.n_abs_value$n_control_surv <- meta.final.n_abs_value$Control_Mean
meta.final.n_abs_value$n_control_no_surv <- (meta.final.n_abs_value$control_n_total) - 
  (meta.final.n_abs_value$Control_Mean)

meta.final.n_abs_value$n_exp_surv <- meta.final.n_abs_value$Exp_Mean
meta.final.n_abs_value$n_exp_no_surv <- (meta.final.n_abs_value$exp_n_total) - 
  (meta.final.n_abs_value$Exp_Mean)


################################################################################
# 2) Percentage
meta.final.n_abs_perc <- droplevels(subset(meta.final.n_abs, 
                                           meta.final.n_abs$Var_type_comment == "Absolute percentage; no variance"))
nrow(meta.final.n_abs_perc)


# All percentages are positive (i.e., survival, fledgling, ...).
# The function I need to use for success is:
# (value reported * N total)/100


################################################################################
# control
meta.final.n_abs_perc$n_control_surv <- ((meta.final.n_abs_perc$Control_Mean * 
                                            meta.final.n_abs_perc$control_n_total))/100

# The function for those clutches/offspring that were not successful is:
meta.final.n_abs_perc$n_control_no_surv <- (meta.final.n_abs_perc$control_n_total) - 
  (meta.final.n_abs_perc$n_control_surv)


################################################################################
# exp
meta.final.n_abs_perc$n_exp_surv <- ((meta.final.n_abs_perc$Exp_Mean * 
                                        meta.final.n_abs_perc$exp_n_total))/100

# The function for those clutches/offspring that were not successful is:
meta.final.n_abs_perc$n_exp_no_surv <- (meta.final.n_abs_perc$exp_n_total) - 
  (meta.final.n_abs_perc$n_exp_surv)


################################################################################
# 3) Percentage and absolute value (it's the same as percentage, but here I can
# check that indeed all values are alright)
meta.final.n_abs_perc_and_abs <- droplevels(subset(meta.final.n_abs, 
                                                   meta.final.n_abs$Var_type_comment == "Absolute percentage (also abs value); no variance"))
nrow(meta.final.n_abs_perc_and_abs)


################################################################################
# control
meta.final.n_abs_perc_and_abs$n_control_surv <- ((meta.final.n_abs_perc_and_abs$Control_Mean * 
                                                    meta.final.n_abs_perc_and_abs$control_n_total))/100

# The function for those clutches/offspring that were not successful is:
meta.final.n_abs_perc_and_abs$n_control_no_surv <- (meta.final.n_abs_perc_and_abs$control_n_total) - 
  (meta.final.n_abs_perc_and_abs$n_control_surv)

################################################################################
# exp
meta.final.n_abs_perc_and_abs$n_exp_surv <- (meta.final.n_abs_perc_and_abs$Exp_Mean * 
                                               meta.final.n_abs_perc_and_abs$exp_n_total)/100

# The function for those clutches/offspring that were not successful is:
meta.final.n_abs_perc_and_abs$n_exp_no_surv <- (meta.final.n_abs_perc_and_abs$exp_n_total) - 
  (meta.final.n_abs_perc_and_abs$n_exp_surv)


################################################################################
# 4) Proportion
meta.final.n_abs_prop <- droplevels(subset(meta.final.n_abs, 
                                           meta.final.n_abs$Var_type_comment == "Absolute proportion; no variance"))
nrow(meta.final.n_abs_prop)


# All proportions are positive (i.e., survival).
# The function I need to use for success is:
# (value reported * N total)/1

################################################################################
# control
meta.final.n_abs_prop$n_control_surv <- (meta.final.n_abs_prop$Control_Mean * 
                                           meta.final.n_abs_prop$control_n_total)

# The function for those clutches/offspring that were not successful is:
meta.final.n_abs_prop$n_control_no_surv <- (meta.final.n_abs_prop$control_n_total) - 
  (meta.final.n_abs_prop$n_control_surv)

################################################################################
# exp
meta.final.n_abs_prop$n_exp_surv <- (meta.final.n_abs_prop$Exp_Mean * 
                                       meta.final.n_abs_prop$exp_n_total)

# The function for those clutches/offspring that were not successful is:
meta.final.n_abs_prop$n_exp_no_surv <- (meta.final.n_abs_prop$exp_n_total) - 
  (meta.final.n_abs_prop$n_exp_surv)


################################################################################
### Now I'll combined all subsets before calculating the final effect size

meta.final.n_abs_ok <- rbind(meta.final.n_abs_value, meta.final.n_abs_perc, 
                             meta.final.n_abs_perc_and_abs, meta.final.n_abs_prop)


################################################################################

### Here again, I have shared controls (n = 9), so we will use the same approach
# as used above for dealing with shared-control non-independence

# shared controls
meta.final.n_abs_ok_shared.c <- droplevels(subset(meta.final.n_abs_ok, 
                                                  meta.final.n_abs_ok$Shared_control == "yes"))
meta.final.n_abs_ok_shared.c$n_control_surv

# surv
meta.final.n_abs_ok_shared.c <- ddply(meta.final.n_abs_ok_shared.c,
                                      c("Shared_control_ID"),
                                      transform, 
                                      n_control_surv = n_control_surv /
                                        (length(Shared_control_ID)))
# not survived
meta.final.n_abs_ok_shared.c <- ddply(meta.final.n_abs_ok_shared.c,
                                      c("Shared_control_ID"),
                                      transform, 
                                      n_control_no_surv = n_control_no_surv /
                                        (length(Shared_control_ID)))

# not shared controls
meta.final.n_abs_ok_noshared.c <- droplevels(subset(meta.final.n_abs_ok, 
                                                    meta.final.n_abs_ok$Shared_control == "no"))
meta.final.n_abs_ok_noshared.c$n_control_surv

# I merge both data bases
meta.final.n_abs_ok_i <- rbind(meta.final.n_abs_ok_shared.c, 
                               meta.final.n_abs_ok_noshared.c)

# function obtained from Lajeunesse, 2013 is:
# r <- [(db$exp_surv * db$control_no_surv) - (db$control_surv*db$exp_no_surv)]/
#  sqrt[(db$exp_surv + db$control_surv)*(db$control_no_surv+db$exp_no_surv)*
#     (db$exp_surv+db$exp_no_surv)*(db$control_surv+db$control_no_surv)]

meta.final.n_abs_ok_i$cor <- ((meta.final.n_abs_ok_i$n_exp_surv * 
                                 meta.final.n_abs_ok_i$n_control_no_surv) - 
                                (meta.final.n_abs_ok_i$n_control_surv * 
                                   meta.final.n_abs_ok_i$n_exp_no_surv))/
  (sqrt((meta.final.n_abs_ok_i$n_exp_surv + 
           meta.final.n_abs_ok_i$n_control_surv) *
          (meta.final.n_abs_ok_i$n_control_no_surv + 
             meta.final.n_abs_ok_i$n_exp_no_surv) *
          (meta.final.n_abs_ok_i$n_exp_surv +
             meta.final.n_abs_ok_i$n_exp_no_surv) *
          (meta.final.n_abs_ok_i$n_control_surv +
             meta.final.n_abs_ok_i$n_control_no_surv)))

meta.final.n_abs_ok_i$cor_var <- ((1 - meta.final.n_abs_ok_i$cor ^ 2) ^ 2/
                                    ((meta.final.n_abs_ok_i$control_n_total +
                                        meta.final.n_abs_ok_i$exp_n_total) - 1))

# Checking that the signs reflect the general criteria used: 
expected.sign_abs <- ifelse(meta.final.n_abs_ok_i$Control_Mean > 
                              meta.final.n_abs_ok_i$Exp_Mean, -1, 1)

observed.sign_abs <- ifelse(meta.final.n_abs_ok_i$cor > 0, 1, -1)

# Checking if they agree
table(expected.sign_abs == observed.sign_abs) 
# It looks ok. Only one value appears as FALSE which corresponds to a row where
# control and experimental groups had the same mean values and cor is 0. 


meta.final.n_abs_ok_i$final_n <- as.numeric(meta.final.n_abs_ok_i$control_n_total) +
  as.numeric(meta.final.n_abs_ok_i$exp_n_total)


################################################################################
### I know merge all the data bases together (the ones that had stats, the one
# with mean values and the one with absolute values)
################################################################################

meta.final_ok <- rbind(meta.final.n_chi, meta.final.n_ftest, 
                       meta.final.n_spearman, 
                       # meta.final.n_tt, reminder: could not included because of lack of info on directionality
                       meta.final.n_z, meta.final.n_pearson_ok,
                       meta.final.n_mean_ok, meta.final.n_abs_ok_i)
nrow(meta.final_ok)
  

# Final touch, combining the two following characters into one so that the info
# about the data location is all in one single variable, then deleting the
# separate variables
meta.final_ok$DataSourceLocation <- paste(meta.final_ok$DataSource,
                                          sep = "; ",
                                          meta.final_ok$DataLocation)

meta.final_ok <- meta.final_ok[ , - which(names(meta.final_ok) %in% 
                                            c("DataSource","DataLocation"))]


################################################################################
# 3. Preparing the moderators for the models 
################################################################################

# Filling EffectID so that we can model unit-level heterogeneity (i.e. model
# within-study/residual variance)
meta.final_ok$EffectID <- 1:nrow(meta.final_ok)

# First check that the last variables we created have the correct format
# If this is not the case, modify it.
meta.final_ok$stats_n_df <- as.numeric(meta.final_ok$stats_n_df) 
meta.final_ok$raw_n <- as.numeric(meta.final_ok$raw_n) 
meta.final_ok$cor <- as.numeric(meta.final_ok$cor) 
meta.final_ok$cor_var <- as.numeric(meta.final_ok$cor_var) 
meta.final_ok$sd_control <- as.numeric(meta.final_ok$sd_control) 
meta.final_ok$sd_exp <- as.numeric(meta.final_ok$sd_exp) 
meta.final_ok$n_control_surv <- as.numeric(meta.final_ok$n_control_surv) 
meta.final_ok$n_exp_surv <- as.numeric(meta.final_ok$n_exp_surv) 
meta.final_ok$n_exp_no_surv <- as.numeric(meta.final_ok$n_exp_no_surv) 
meta.final_ok$final_n <- as.numeric(meta.final_ok$final_n) 

### We will reorder the names of the fitness traits so they appear in a 
# chronological order.
meta.final_ok$Fitness_mother <- factor(meta.final_ok$Fitness_mother, 
                                       levels = c("clutch size", 
                                                  "hatching success",
                                                  "hatching probability",
                                                  "hatching failure",
                                                  "egg mortality",
                                                  "maternal longevity",
                                                  "offspring recruit",
                                                  "maternal life time reproductive success"))

# Same for offspring fitness traits
meta.final_ok$Fitness_offs <- factor(meta.final_ok$Fitness_offs, 
                                     levels = c("mass", 
                                                "tarsus",
                                                "wing",
                                                "head length",
                                                "head-bill length",
                                                "culmen", 
                                                "flipper length", 
                                                "gape width",
                                                "beak flank width", 
                                                "structural body size (mass and tarsus)",
                                                "structural body size (mass and bill length)",
                                                "growth (mass)", 
                                                "growth (PC1)", 
                                                "growth (PC2)", 
                                                "growth (tarsus)",                           
                                                "growth rate", 
                                                "growth rate (flipper length)",
                                                "growth rate (mass)", 
                                                "growth rate (tarsus)", 
                                                "growth rate (ulna)", 
                                                "mass gain", 
                                                "fledgling number", 
                                                "fledgling probability",
                                                "fledgling success",  
                                                "fledgling success (fledglings/hatchings)",  
                                                "pre-fledgling survival probability",
                                                "recruitment", 
                                                "clutch size",
                                                "hatching number", 
                                                "offspring mortality", 
                                                "offspring survival",
                                                "offspring survival years"))


################################################################################
# We now need to create a new variable that contains information on the hormone 
# measured. At the moment, the data base has this information but separated in 
# 3 different columns (i.e., info for correlative vs experimental studies, and 
# within experimental, mother and egg). We will merge this information into one.

# For this, we first create the new column where we will put all the info. 
meta.final_ok$Hormone_measured_unique <- "NA"

# We will subset the data base to get those rows that contain the required info.
meta.final_ok$Horm_measured <- as.character(meta.final_ok$Horm_measured)
meta.final_ok$Horm_mother_exp <- as.character(meta.final_ok$Horm_mother_exp)
meta.final_ok$Horm_egg_exp <- as.character(meta.final_ok$Horm_egg_exp)

# hormones measured in correlative studies
meta.final_ok_cor <- droplevels(subset(meta.final_ok, 
                                       meta.final_ok$Horm_measured != "NA"))
meta.final_ok_cor$Hormone_measured_unique <- meta.final_ok_cor$Horm_measured

meta.final_ok_cor_yes <- droplevels(subset(meta.final_ok_cor, 
                                           meta.final_ok_cor$Study_type == "correlational"))

# hormones measured in experimental studies - mothers
meta.final_ok_exp_mother <- droplevels(subset(meta.final_ok, 
                                              meta.final_ok$Horm_mother_exp != "NA"))
meta.final_ok_exp_mother$Hormone_measured_unique <- meta.final_ok_exp_mother$Horm_mother_exp

# hormones measured in experimental studies - eggs
meta.final_ok_exp_egg <- droplevels(subset(meta.final_ok, 
                                           meta.final_ok$Horm_egg_exp != "NA"))
meta.final_ok_exp_egg$Hormone_measured_unique <- meta.final_ok_exp_egg$Horm_egg_exp

# I now merge all data bases
meta.final_ok_ok <- rbind(meta.final_ok_cor_yes, 
                          meta.final_ok_exp_mother,
                          meta.final_ok_exp_egg)
nrow(meta.final_ok_ok)


# IMPORTANT: in the pre-registration we wrote "Furthermore, models including 
# categorical moderators will require a minimum of 5 data points per moderator 
# level, or else the model(s) will be run without that specific moderator level 
# unless the level can be meaningfully merged with another moderator level." 
# Flutamide has only 1 row with information. Furthermore, Flutamide is not a 
# hormone, but a blocker. And we did not have predictions in our pre-registrations 
# about blockers. Hence, this data point will not be included in the final data base. 
table(meta.final_ok_ok$Hormone_measured_unique)

meta.final_ok_ok <- droplevels(subset(meta.final_ok_ok, 
                                      meta.final_ok_ok$Hormone_measured_unique != "flutamide"))

nrow(meta.final_ok_ok)


# Important note: There are 6 rows that have a correlation value = 0. We checked these
# rows for potential mistakes, but the values are correct. The values occur
# because both control and experimental groups have the same mean or absolute
# value (variances are different though).

################################################################################
# Second, I create a column that groups together hormone types into bigger 
# categorical variables. 
meta.final_ok_ok$Hormone_measured_general <- "NA"

meta.final_ok_ok$Hormone_measured_general <- meta.final_ok_ok$Hormone_measured_unique

# Androgens are represented by: androstenedione, testosterone, DHT, and 
# androstenedione and testosterone.
meta.final_ok_ok$Hormone_measured_general <- str_replace(meta.final_ok_ok$Hormone_measured_general,
                                                         "androstenedione; dihydrotestosterone and testosterone",
                                                         "androgens")
meta.final_ok_ok$Hormone_measured_general <- str_replace(meta.final_ok_ok$Hormone_measured_general,
                                                         "androstenedione and testosterone", 
                                                         "androgens")
meta.final_ok_ok$Hormone_measured_general <- str_replace(meta.final_ok_ok$Hormone_measured_general,
                                                         "androstenedione", 
                                                         "androgens")
meta.final_ok_ok$Hormone_measured_general <- str_replace(meta.final_ok_ok$Hormone_measured_general,
                                                         "testosterone", 
                                                         "androgens")
meta.final_ok_ok$Hormone_measured_general <- str_replace(meta.final_ok_ok$Hormone_measured_general,
                                                         "DHT", 
                                                         "androgens")

# I transform the variables into factors:
meta.final_ok_ok$Hormone_measured_unique <- as.factor(meta.final_ok_ok$Hormone_measured_unique)
meta.final_ok_ok$Hormone_measured_general <- as.factor(meta.final_ok_ok$Hormone_measured_general)


################################################################################
# In our pre-registration, we stated that for some hypothesis the consequences of
# fitness should occur at the offspring and/or mother level. Hence, I need to split 
# the general data base into two: one for offspring and one for mother. Like this, 
# I will have 3 data bases (off, mother, off&mother) that I will use depending on
# the predictions made.

# I start with the offspring.
meta.final_ok_ok_off <- droplevels(subset(meta.final_ok_ok, 
                                          meta.final_ok_ok$Fitness_trait != "mother"))
nrow(meta.final_ok_ok_off)


meta.final_ok_ok$cor

meta.final_ok_ok_mom <- droplevels(subset(meta.final_ok_ok, 
                                          meta.final_ok_ok$Fitness_trait == "mother"))
nrow(meta.final_ok_ok_mom)



################################################################################
# Saving datasets for further analyses
################################################################################

# To save the final data base, run:
write.csv(meta.final_ok_ok,
          'data/processed_data/meta-analysis_egg_hormones_processed_data_full.csv', 
          row.names = F)

write.csv(meta.final_ok_ok_off,
          'data/processed_data/meta-analysis_egg_hormones_processed_data_offspring.csv', 
          row.names = F)

write.csv(meta.final_ok_ok_mom,
          'data/processed_data/meta-analysis_egg_hormones_processed_data_mom.csv', 
          row.names = F)
