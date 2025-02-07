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

# Script used to analyse the meta-analytic data for the following meta-analysis:

# Do egg hormones have fitness consequences in wild birds? A systematic review 
# and meta-analysis

# by Lucia Mentesana, Michaela Hau, Pietro B. D'Amelio, Nicolas M. Adreani, and
# Alfredo Sanchez-Tojar. 

# preprint available at: https://doi.org/10.1101/2024.10.29.620852

# Note that for all figures – they were exported as .svg format using R's export 
# function. Then, we used the free vector graphics editor Inkscape to enhance 
# their appearance. Specifically, we adjusted the size and font of the legends 
# and modified the colors of the data. Finally, we exported each figure in .png 
# format, which is the version we submitted to the journal.

################################################################################
# Packages needed
################################################################################

# install.packages("pacman") #if not already installed
pacman::p_load(rotl,
               ape,
               ggplot2,
               ggcorrplot,
               metafor,
               orchaRd,
               wesanderson,
               patchwork,
               tidyverse,
               ggtree)

# cleaning up
rm(list = ls())


################################################################################
# Functions: from https://github.com/Yefeng0920/heterogeneity_ecoevo/tree/main/function
################################################################################

# # function to calculate heterogeneity
# h.calc <- function(mod){
#   # I2
#   # sigma2_v = typical sampling error variance
#   sigma2_v <- sum(1 / mod$vi) * (mod$k - 1) /
#     (sum(1 / mod$vi)^2 - sum((1 / mod$vi)^2))
#   # s^2_t = total variance
#   I2_total <- 100 * (sum(mod$sigma2) / (sum(mod$sigma2) + sigma2_v))
#   I2_each <- 100 * (mod$sigma2 / (sum(mod$sigma2) + sigma2_v))
#   I2s_Shinichi <- c(I2_total, I2_each)
#   
#   # matrix version  
#   W <- solve(mod$V)
#   X <- model.matrix(mod)
#   P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
#   I2_total2 <- 100* (sum(mod$sigma2) / (sum(mod$sigma2) + (mod$k - mod$p) / sum(diag(P))))
#   I2_each2 <- 100* (mod$sigma2 / (sum(mod$sigma2) + (mod$k - mod$p) / sum(diag(P))))
#   I2s_Wolfgang <- c(I2_total2, I2_each2)
#   
#   # CVB
#   CV_total <- (sqrt(sum(mod$sigma2)) / abs(mod$beta[1]))
#   CV_each <- (sqrt(mod$sigma2) / abs(mod$beta[1]))
#   CVs <- c(CV_total, CV_each)
#   
#   # M1
#   M1_total <- (sum(sqrt(mod$sigma2)) / (sum(sqrt(mod$sigma2)) + abs(mod$beta[1])))
#   M1_each <- sqrt(mod$sigma2) / (sum(sqrt(mod$sigma2)) + abs(mod$beta[1]))
#   Ms <- c(M1_total, M1_each)
#   
#   hs <- data.frame(I2s_Shinichi,CVs,Ms)
#   rownames(hs) <- c("Total", mod$s.names)
#   return(hs)
#   
# }

# function to calculate heterogeneity - squared version, including CVH2 and M2
h.calc2 <- function(mod){
  # I2
  # sigma2_v = typical sampling error variance
  sigma2_v <- sum(1 / mod$vi) * (mod$k - 1) /
    (sum(1 / mod$vi)^2 - sum((1 / mod$vi)^2))
  # s^2_t = total variance
  I2_total <- 100 * (sum(mod$sigma2) / (sum(mod$sigma2) + sigma2_v))
  I2_each <- 100 * (mod$sigma2 / (sum(mod$sigma2) + sigma2_v))
  #names(I2_each) <- paste0("I2_", model$s.names)
  #names(I2_total) <- "I2_Total"
  I2s_Shinichi <- c(I2_total, I2_each)
  
  # matrix version  
  W <- solve(mod$V)
  X <- model.matrix(mod)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2_total2 <- 100* (sum(mod$sigma2) / (sum(mod$sigma2) + (mod$k - mod$p) / sum(diag(P))))
  I2_each2 <- 100* (mod$sigma2 / (sum(mod$sigma2) + (mod$k - mod$p) / sum(diag(P))))
  #names(I2_each2) <- paste0("I2_", model$s.names)
  #names(I2_total2) <- "I2_Total2"
  I2s_Wolfgang <- c(I2_total2, I2_each2)
  
  
  # CVH2
  CV_total <- (sum(mod$sigma2) / (mod$beta[1])^2)
  CV_each <- (mod$sigma2 / (mod$beta[1])^2)
  
  #names(CVB_each) <- paste0("CVB_", mod$s.names)
  #names(CVB_total) <- "CVB_total"
  CVHs <- c(CV_total, CV_each)
  
  # M2
  M_total <- (sum(mod$sigma2) / (sum(mod$sigma2) + (mod$beta[1])^2))
  M_each <- (mod$sigma2) / (sum(mod$sigma2) + (mod$beta[1])^2)
  Ms <- c(M_total, M_each)
  
  hs <- data.frame(I2s_Shinichi,CVHs,Ms)
  rownames(hs) <- c("Total", mod$s.names)
  return(hs)
}

# function to estimate typical sampling error variance
sigma2_v <- function(mod){
  sigma2_v <- sum(1 / mod$vi) * (mod$k - 1) /
    (sum(1 / mod$vi)^2 - sum((1 / mod$vi)^2))
  return(sigma2_v)
}

################################################################################
# Importing the raw data 
################################################################################

meta.final_ok_ok <- read.csv("data/processed_data/meta-analysis_egg_hormones_processed_data_full.csv")
meta.final_ok_ok_off <- read.csv("data/processed_data/meta-analysis_egg_hormones_processed_data_offspring.csv")
meta.final_ok_ok_mom <- read.csv("data/processed_data/meta-analysis_egg_hormones_processed_data_mom.csv")

################################################################################
# 4. General information on Statistical Analysis:
################################################################################

# We tested each hypothesis following what we wrote in our pre-registration.

# To run the statistical models, we first created a phylogenetic tree and 
# a covariance matrix for each data base, we then run the statistical models, and
# we finally got a table with results and a figure. Below we provide the reader
# with detailed information on each of these steps:

# PHYLOGENETIC TREE: we followed the methods used by Sanchez-Tojar et 
# al. 2020 - 'The jury is still out regarding the generality of adaptive 
# "transgenerational" effects'. 

# SAMPLING VARIANCE-COVARIANCE MATRIX: We specified sampling variance as a 
# variance-covariance matrix, with the sampling variance for each effect size on 
# the diagonal, and the covariance between these measures as off-diagonal elements. 
# The model assumed a 0.5 correlation between the effect size sample variances 
# with the same StudyID. For a similar approach, see O'Dea et al., 2019
# (https://doi.org/10.1111/faf.12394) or Kim et al. 2022 
# (https://doi.org/10.1111/1365-2656.13554), from which this code was taken from


# It is important to note that each statistical model required us to use 
# different data bases. Whenever this was the case, and only if necessary, we 
# created a new phylogenetic tree and/or covariance matrix, and made a note in 
# the script.

# STATISTICAL MODELS: Group ID and Study ID overlap completely. Following
# what we wrote in our pre-registration, we therefore did not include Group ID 
# in the models. Same with LaboratoryID and Lab_PI (I only kept LaboratoryID).

# TABLE WITH RESULTS AND FIGURE: we used the orchaRd package from Nakawaga et al. 
# (2023); "orchaRd 2.0: An R package for visualising meta-analyses with orchard 
# plots"; doi:  10.1111/2041-210X.14152).

################################################################################
# 5. Meta-analytic model (i.e. intercept-only model)

# This model indicates the overall effect that egg hormones have on fitness

# Data base: general data base (i.e., includes maternal and offspring fitness)
################################################################################

# PHYLOGENETIC TREE - Species names:

# We first searched for species names in the in the Open Tree Taxonomy [@rees2017]
# using the R package `rotl` v.3.0.5 [@michonneau2016]. This allowed us to confirm
# that all species names were correct and that synonyms or typos were not present
# in the database.
# resolved_names <- tnrs_match_names(as.character(unique(meta.final_ok_ok$Species)))
# #
# # Saving the taxonomic data created on the 5th Feb 2025 to speed the
# # process in the future and allow full reproducibility
# save(resolved_names, file = "data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life.RData")

# Loading the taxonomic data created on the 5th Feb 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life.RData") #resolved_names

# # extracting phylogenetic information
# my_tree <- tol_induced_subtree(ott_ids =
#                                resolved_names[,"ott_id"],
#                               label_format = "name")
# # Quick tree plotting
# # plot(my_tree, no.margin = TRUE)
# #
# # We need to check for the existence of polytomies
# is.binary(my_tree)
# # Yes, meaning there are no polytomies. Let's go on.
# 
# # To confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# my_tree$tip.label <- gsub("_", " ", my_tree$tip.label)
# #
# intersect(as.character(my_tree$tip.label),
#            as.character(meta.final_ok_ok$Species))
# 
# # Listed in our database but not in the tree
# setdiff(as.character(meta.final_ok_ok$Species),
#        as.character(my_tree$tip.label))
# #
# # Listed in the tree but not in our database
# setdiff(as.character(my_tree$tip.label),
#          as.character(meta.final_ok_ok$Species))
# # No error or inconsistencies found.
# #
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# # 5th Feb 2025
# save(my_tree, file = "data/outputs/phylogenetic_files/tree.Rdata")

# We can now load the saved tree
load("data/outputs/phylogenetic_files/tree.Rdata") #my_tree

# # Compute branch lengths of tree
# phylo_branch <- compute.brlen(my_tree, method = "Grafen", power = 1)
# 
# # Check if tree is ultrametric
# is.ultrametric(phylo_branch)
# # TRUE
# 
# # Matrix to be included in the models
# phylo_cor <- vcv(phylo_branch, cor = T)
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(phylo_cor, file = "data/outputs/phylogenetic_files/phylo_cor.Rdata")

# we can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor.Rdata") #phylo_cor

# visual exploration of the phylogenetic correlation matrix
ggcorrplot::ggcorrplot(phylo_cor, sig.level = 0.05, lab_size = 1,
                       p.mat = NULL,insig = c("pch", "blank"),
                       pch = 1, pch.col = "black", pch.cex = 1, tl.cex = 1.5) +
  theme(axis.text.x = element_text(size = 7, margin = margin(-2, 0, 0, 0)),
        axis.text.y = element_text(size = 7, margin = margin(0, -2, 0, 0)),
        panel.grid.minor = element_line(size = 3)) +
  geom_tile(fill = "white") +
  geom_tile(height = 0.8, width = 0.8) +
  scale_fill_gradient2(low = "#E69F00",mid = "white", high = "#56B4E9",
                       midpoint = 0.5, breaks = c(0, 1),
                       limit = c(0,1)) + labs(fill = "Correlation")

# Creating a duplicate species variable for the phylogenetic analysis
meta.final_ok_ok$Species_phylo <- meta.final_ok_ok$Species
meta.final_ok_ok$Species_phylo

# VARIANCE-COVARIANCE MATRIX: 
# Creating a var-covar matrix assuming a 0.5 correlation between effect sizes 
# from the same StudyID. covariance = (0.5 * sqrt(vi.1) * sqrt(vi.2))
# Creates a matrix (called 'VCV_ESVar') with the dimensions =  
# n(effect_sizes) x n(effect_sizes)

VCV_ESVar <- matrix(0, nrow = nrow(meta.final_ok_ok), 
                    ncol = nrow(meta.final_ok_ok))

# Names rows and columns for each obsID
rownames(VCV_ESVar) <- meta.final_ok_ok[, "EffectID"]
colnames(VCV_ESVar) <- meta.final_ok_ok[, "EffectID"]

# Finds effect sizes that come from the same study
shared_coord <- which(meta.final_ok_ok[, "StudyID"] %in% 
                        meta.final_ok_ok[duplicated(meta.final_ok_ok[, "StudyID"]), 
                                         "StudyID"] == TRUE)

combinations <- do.call("rbind", tapply(shared_coord, 
                                        meta.final_ok_ok[shared_coord, "StudyID"], 
                                        function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations)[1]) {
  p1 <- combinations[i, 1]
  p2 <- combinations[i, 2]
  p1_p2_cov <- 0.5 * sqrt(meta.final_ok_ok[p1, "cor_var"]) * 
    sqrt(meta.final_ok_ok[p2, "cor_var"])
  VCV_ESVar[p1, p2] <- p1_p2_cov
  VCV_ESVar[p2, p1] <- p1_p2_cov
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(VCV_ESVar) <- meta.final_ok_ok[, "cor_var"]

# # In case you want to visually double check the matrix outside of R
# write.csv(VCV_ESVar, 'data/outputs/variance-covariance_matrices/VCV_ESVar.csv')


################################################################################
# Main effect model

# STATISTICAL MODEL:
meta.model <- rma.mv(cor,
                     VCV_ESVar,
                     mods = ~ 1,
                     random = list(~ 1 | StudyID,
                                   ~ 1 | LaboratoryID,
                                   ~ 1 | PopulationID,
                                   ~ 1 | Species,
                                   ~ 1 | Species_phylo,
                                   ~ 1 | EffectID),
                     method = "REML",
                     R = list(Species_phylo = phylo_cor),
                     test = "t",
                     data = meta.final_ok_ok)

#save(meta.model, file = "data/outputs/statistical_models/meta_model.Rdata")
load(file = "data/outputs/statistical_models/meta_model.Rdata") #meta.model

# Printing the summary results of the model
print(meta.model, digits = 3)
# There is an overall negative effect (-0.072) not statistically significant.

# Printing the results again, but adding the credibility/prediction interval, 
# which uses the heterogeneity to generate an interval that should contain 95%
# of the effect sizes of any future or unknown studies with similar features
# as those included in the current database
predict(meta.model, digits = 3)

# Model funnel plots
par(mfrow = c(1, 1))
funnel(meta.model)

################################################################################
# Heterogeneity quantification and stratification
################################################################################

# Typical sampling variance (which captures the statistical noise of the data, 
# but is rarely reported in the current meta-analytic practice):
sigma2_v(meta.model)


# Total unstandardized raw variance (i.e. total heterogeneity)
round(sum(meta.model$sigma2), 3)


# # I2, CV and M
# round(h.calc(meta.model),2)
round(i2_ml(meta.model), 3)
round(orchaRd::cvh2_ml(meta.model), 3)
round(orchaRd::m2_ml(meta.model), 3)


# Visualize heterogeneity
## make dataframe
h_status <- h.calc2(meta.model)

# adding sigmas
h_status$sigma2s <- c(sum(meta.model$sigma2),
                      meta.model$sigma2[1],
                      meta.model$sigma2[2],
                      meta.model$sigma2[3],
                      meta.model$sigma2[4],
                      meta.model$sigma2[5],
                      meta.model$sigma2[6])
round(h_status$sigma2s, 3)

h_status$levels <- rownames(h_status)  
h_status$levels <- dplyr::recode(h_status$levels, 
                                 "Total" = "Total",  
                                 "EffectID" = "Within",  
                                 "StudyID" = "Among", 
                                 "Species" = "Spp", 
                                 "Species_phylo" = "Phylo",
                                 "LaboratoryID" = "Lab",
                                 "PopulationID" = "Pop")

h_status$levels <- as.factor(h_status$levels)

h_status$levels <- factor(h_status$levels, levels = c("Total", "Phylo", 
                                                      "Spp",
                                                      "Pop",
                                                      "Lab",
                                                      "Among", 
                                                      "Within"))

# plot

# sigma
p.sigma <- ggplot(h_status, aes(levels, sigma2s)) +
  geom_col(alpha = 1,
           color = wes_palette('GrandBudapest1', 4, type = 'discrete')[1],
           fill = wes_palette('GrandBudapest1', 4, type = 'discrete')[1]) +
  labs(y = expression("Variance"), x = "" ,
       title = "Unstandardised heterogeneity metrics") +
  theme_bw() +
  theme(legend.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black")
  )

# I2
p.I2 <- ggplot(h_status, aes(levels, I2s_Shinichi/100)) +
  geom_col(alpha = 1, color = wes_palette('GrandBudapest1', 4, 
                                          type = 'discrete')[2], 
           fill = wes_palette('GrandBudapest1', 4, type = 'discrete')[2]) +
  scale_y_continuous(labels = scales::percent_format()
  ) +
  labs(y = expression(paste(italic(I)[]^2)), x = "" , 
       title = "Source of heterogeneity") + 
  theme_bw() +
  theme(legend.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black")
  ) 

# CV
p.CV <- ggplot(h_status, aes(levels, CVHs)) +
  geom_col(alpha = 1, color = wes_palette('GrandBudapest1', 4, 
                                          type = 'discrete')[3], 
           fill = wes_palette('GrandBudapest1', 4, type = 'discrete')[3]) +
  labs(y = expression(paste(italic(CV)[])), x = "" , 
       title = "Magnitude of heterogeneity") + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + # https://stackoverflow.com/questions/38722202/how-do-i-change-the-number-of-decimal-places-on-axis-labels-in-ggplot2
  theme_bw() +
  theme(legend.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black")
  )   


# M
p.M <- ggplot(h_status, aes(x = levels, y = Ms)) +
  geom_col(alpha = 1, color = wes_palette('GrandBudapest1', 4, 
                                          type = 'discrete')[4], fill = wes_palette('GrandBudapest1', 4, type = 'discrete')[4]) +
  labs(y = expression(paste(italic(M)[])), x = "" , 
       title = "Magnitude of heterogeneity") + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + # https://stackoverflow.com/questions/38722202/how-do-i-change-the-number-of-decimal-places-on-axis-labels-in-ggplot2
  theme_bw() +
  theme(legend.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black")
  )  


p.sigma + p.I2 + p.CV + p.M + plot_layout(ncol = 2, nrow = 2) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold"))


################################################################################
# Main model results

# TABLE
results_meta.model <- mod_results(meta.model, 
                                  mod = "1", 
                                  at = NULL, group = "StudyID")


# FIGURE: overall effect of yolk hormones on fitness
figure_meta.model <- orchaRd::orchard_plot(meta.model, mod = "1", 
                                           group = "Fitness_trait", 
                                           xlab = "Effect size",
                                           transfm = "none",
                                           trunk.size = 2,
                                           branch.size = 2,
                                           colour = TRUE,  
                                           twig.size = 1)

# Important: in the figure there are four colors. This is because we have 
# 4 levels at the moment: mother, offspring, fem and male offspring when 
# recaptured as breeding adults. 
# For the msc, I will modify the colors of the figure so that we have only 
# to two categories (mother and offspring). Like this we can easily see that, irrespective
# of whether the fitness trait was of maternal or offspring fitness, the overall
# null effect represents quite well both categories.

# Overall, egg hormones have a negative effect on fitness traits (-0.073). This 
# effect is not statistically significant.


# I also create a figure to inspect effect sizes for each species. 
# figure_meta.model_spp <- orchaRd::orchard_plot(meta.model, mod = "1",
#                                                group = "Species_phylo",
#                                                xlab = "Effect size",
#                                                transfm = "none",
#                                                trunk.size = 8,
#                                                branch.size = 2,
#                                                colour = TRUE,
#                                                twig.size = 1)


# I will plot the phylogenetic tree and the effect sizes found for each species.
# We opted for doing this, because the results concerning heterogeneity of 
# random effects suggest that 'phylogeny' explains a good amount of differences
# across effect sizes. 
# By plotting the phylogenetic tree together with effect sizes, we can visually
# inspect, and better interpret, the results obtained.

taxa <- resolved_names
ott_in_tree <- ott_id(taxa)[is_in_tree(ott_id(taxa))]
length(ott_id(taxa)) - length(is.na(ott_in_tree)) # all good

# make phylo tree
tree <- suppressWarnings(tol_induced_subtree(ott_ids = ott_id(taxa)))

tree$tip.label <- strip_ott_ids(tree$tip.label, remove_underscores = TRUE)

tree <- compute.brlen(tree, method = "Grafen", power = 1)

tree_matrix <- vcv.phylo(tree, model = "Brownian", corr = T)

ggcorrplot::ggcorrplot(tree_matrix, sig.level = 0.05, lab_size = 4.5, p.mat = NULL,
                       insig = c("pch", "blank"), pch = 1, pch.col = "black", pch.cex = 1, tl.cex = 14) +
  theme(axis.text.x = element_text(size = 10, margin = margin(-2, 0, 0, 0)), axis.text.y = element_text(size = 10,
                                                                                                        margin = margin(0, -2, 0, 0)), panel.grid.minor = element_line(size = 10)) +
  geom_tile(fill = "white") + geom_tile(height = 0.8, width = 0.8) + scale_fill_gradient2(low = "#E69F00",
                                                                                          mid = "white", high = "#56B4E9", midpoint = 0.5, breaks = c(0, 1), limit = c(0,
                                                                                                                                                                       1)) + labs(fill = "Correlation")

meta.final_ok_ok_checking <- meta.final_ok_ok

meta.final_ok_ok_checking <- meta.final_ok_ok_checking %>%
  mutate(search_string = tolower(Species_phylo))

meta.final_ok_ok_checking <- left_join(meta.final_ok_ok_checking, 
                                       dplyr::select(taxa, search_string, unique_name,
                                                     ott_id), by = "search_string")

meta.final_ok_ok_checking <- meta.final_ok_ok_checking %>%
  mutate(spp = unique_name, phylo = unique_name)


db_effect <- data.frame(Species = meta.final_ok_ok_checking$Species_phylo,
                        Cor = meta.final_ok_ok_checking$cor,
                        Cor_var = meta.final_ok_ok_checking$cor_var)

tip.label <- data.frame(Species_phylo = tree$tip.label)

spp <- dplyr::distinct(meta.final_ok_ok_checking, Species_phylo, .keep_all = TRUE) %>%
  dplyr::select(Species_phylo)
spp.2 <- left_join(tip.label, spp, by = "Species_phylo")

# Phylogenetic tree
tree.p1 <- ggtree(tree, layout = "rectangular", cex = 0.4)

# Phylogenetic tree and species name
tree.p2 <- tree.p1 %<+% spp.2 + geom_tiplab(aes(color = "Species_phylo"), size = 3, fontface = "italic",
                                            align = T, offset = 0.05) + geom_tippoint(aes(color = "Species_phylo")) + guides(color = "none") +
  xlim_expand(xlim = c(0, 1.8), panel = "Tree") + scale_color_viridis_d()

# Phylogenetic tree and effect sizes
tree.p3 <- tree.p2 + geom_facet(panel = "Effect size", data = db_effect, 
                                geom = ggstance::geom_pointrangeh,
                                mapping = aes(x = Cor, xmin = Cor_var, xmax = Cor_var, color = "Species_phylo")) +
  theme_tree2() + theme(strip.background = element_rect(fill = "white")) + guides(fill = "none",
                                                                                  color = "none")

#################################################################################
# Main effect model - OFFSPRING FITNESS TRAITS ONLY

# Note that this is non-pre-registered intercept only model that we used to 
# confirm that the overall relationship seen for both mother and offspring fitness
# traits, hold when looking solely at offspring traits. 

################################################################################
# First, I create the phylogenetic tree and matrix.

# # PHYLOGENETIC TREE:
# resolved_names_off <- tnrs_match_names(as.character(unique(meta.final_ok_ok_off$Species)))
# 
# # Saving the taxonomic data created on the 5th Feb 2025 to speed the
# # process in the future and allow full reproducibility
# save(resolved_names_off, file = "data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_off.RData")

# Loading the taxonomic data created on the 5th Feb 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_off.RData") #resolved_names_off

# # extracting phylogenetic information
# my_tree_off <- tol_induced_subtree(ott_ids =
#                                      resolved_names_off[,"ott_id"],
#                                    label_format = "name")
# 
# # # Quick tree plotting
# # plot(my_tree_off, no.margin = TRUE)
# 
# # We need to check for the existence of polytomies
# is.binary(my_tree_off)
# # Yes, meaning there are no polytomies. Let's go on.
# 
# # To confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# my_tree_off$tip.label <- gsub("_", " ", my_tree_off$tip.label)
# 
# intersect(as.character(my_tree_off$tip.label),
#           as.character(meta.final_ok_ok_off$Species))
# 
# # Listed in our database but not in the tree
# setdiff(as.character(meta.final_ok_ok_off$Species),
#         as.character(my_tree_off$tip.label))
# 
# # Listed in the tree but not in our database
# setdiff(as.character(my_tree_off$tip.label),
#         as.character(meta.final_ok_ok_off$Species))
# # No error or inconsistencies found.
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# # 5th Feb 2025
# save(my_tree_off, file = "data/outputs/phylogenetic_files/tree_off.Rdata")

# We can now load the saved tree
load("data/outputs/phylogenetic_files/tree_off.Rdata") #my_tree_off

# # Compute branch lengths of tree
# phylo_branch_off <- compute.brlen(my_tree_off, method = "Grafen", power = 1)
# 
# # Check if tree is ultrametric
# is.ultrametric(phylo_branch_off)
# # TRUE
# 
# # Matrix to be included in the models
# phylo_cor_off <- vcv(phylo_branch_off, cor = T)
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(phylo_cor_off, file = "data/outputs/phylogenetic_files/phylo_cor_off.Rdata")

# we can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor_off.Rdata") #phylo_branch_off

# Creating a duplicate species variable for the phylogenetic analysis
meta.final_ok_ok_off$Species_phylo_off <- meta.final_ok_ok_off$Species



# VARIANCE-COVARIANCE MATRIX:
VCV_ESVar_off <- matrix(0, nrow = nrow(meta.final_ok_ok_off), 
                        ncol = nrow(meta.final_ok_ok_off))

# Names rows and columns for each obsID
rownames(VCV_ESVar_off) <- meta.final_ok_ok_off[, "EffectID"]
colnames(VCV_ESVar_off) <- meta.final_ok_ok_off[, "EffectID"]

# Finds effect sizes that come from the same study
shared_coord_off <- which(meta.final_ok_ok_off[, "StudyID"] %in% 
                            meta.final_ok_ok_off[duplicated(meta.final_ok_ok_off[, "StudyID"]), 
                                                 "StudyID"] == TRUE)

combinations_off <- do.call("rbind", tapply(shared_coord_off, 
                                            meta.final_ok_ok_off[shared_coord_off, "StudyID"], 
                                            function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations_off)[1]) {
  p1_off <- combinations_off[i, 1]
  p2_off <- combinations_off[i, 2]
  p1_p2_cov_off <- 0.5 * sqrt(meta.final_ok_ok_off[p1_off, "cor_var"]) * 
    sqrt(meta.final_ok_ok_off[p2_off, "cor_var"])
  VCV_ESVar_off[p1_off, p2_off] <- p1_p2_cov_off
  VCV_ESVar_off[p2_off, p1_off] <- p1_p2_cov_off
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(VCV_ESVar_off) <- meta.final_ok_ok_off[, "cor_var"]

# # In case you want to visually double check the matrix outside of R
# write.csv(VCV_ESVar_off, 'data/outputs/variance-covariance_matrices/VCV_ESVar_off.csv')


# STATISTICAL MODEL:
meta.model.off <- rma.mv(cor,
                         VCV_ESVar_off,
                         mods = ~ 1,
                         random = list(~ 1 | StudyID,
                                       ~ 1 | LaboratoryID,
                                       ~ 1 | PopulationID,
                                       ~ 1 | Species,
                                       ~ 1 | Species_phylo_off,
                                       ~ 1 | EffectID),
                         method = "REML",
                         R = list(Species_phylo_off = phylo_cor_off),
                         test = "t",
                         data = meta.final_ok_ok_off)

#save(meta.model.off, file = "data/outputs/statistical_models/meta_model.off.Rdata")
load(file = "data/outputs/statistical_models/meta_model.off.Rdata") #meta.model.off

# Printing the summary results of the model
print(meta.model.off, digits = 3)
# There is an overall negative effect (-0.063) not statistically significant.

# Printing the results again, but adding the credibility/prediction interval, 
# which uses the heterogeneity to generate an interval that should contain 95%
# of the effect sizes of any future or unknown studies with similar features
# as those included in the current database
predict(meta.model.off, digits = 3)

# Model funnel plots
par(mfrow = c(1, 1))
funnel(meta.model.off)


################################################################################
# Heterogeneity quantification and stratification - OFFSPRING FITNESS
################################################################################

# Typical sampling variance (which captures the statistical noise of the data, 
# but is rarely reported in the current meta-analytic practice):
sigma2_v(meta.model.off)


# Total unstandardized raw variance (i.e. total heterogeneity)
round(sum(meta.model.off$sigma2), 3)


# # I2, CV and M
# round(h.calc(meta.model.off),2)
round(i2_ml(meta.model.off), 3)
round(cvh2_ml(meta.model.off), 3)
round(m2_ml(meta.model.off), 3)


# Visualize heterogeneity
## make dataframe
h_status <- h.calc2(meta.model.off)

# adding sigmas
h_status$sigma2s <- c(sum(meta.model.off$sigma2),
                      meta.model.off$sigma2[1],
                      meta.model.off$sigma2[2],
                      meta.model.off$sigma2[3],
                      meta.model.off$sigma2[4],
                      meta.model.off$sigma2[5],
                      meta.model.off$sigma2[6])
round(h_status$sigma2s, 3)

h_status$levels <- rownames(h_status)  
h_status$levels <- dplyr::recode(h_status$levels, 
                                 "Total" = "Total",  
                                 "EffectID" = "Within",  
                                 "StudyID" = "Among", 
                                 "Species" = "Spp", 
                                 "Species_phylo" = "Phylo",
                                 "LaboratoryID" = "Lab",
                                 "PopulationID" = "Pop")

h_status$levels <- as.factor(h_status$levels)

h_status$levels <- factor(h_status$levels, levels = c("Total", "Phylo", 
                                                      "Spp",
                                                      "Pop",
                                                      "Lab",
                                                      "Among", 
                                                      "Within"))

# plot

# sigma
p.sigma <- ggplot(h_status, aes(levels, sigma2s)) +
  geom_col(alpha = 1,
           color = wes_palette('GrandBudapest1', 4, type = 'discrete')[1],
           fill = wes_palette('GrandBudapest1', 4, type = 'discrete')[1]) +
  labs(y = expression("Variance"), x = "" ,
       title = "Unstandardised heterogeneity metrics") +
  theme_bw() +
  theme(legend.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black")
  )

# I2
p.I2 <- ggplot(h_status, aes(levels, I2s_Shinichi/100)) +
  geom_col(alpha = 1, color = wes_palette('GrandBudapest1', 4, 
                                          type = 'discrete')[2], 
           fill = wes_palette('GrandBudapest1', 4, type = 'discrete')[2]) +
  scale_y_continuous(labels = scales::percent_format()
  ) +
  labs(y = expression(paste(italic(I)[]^2)), x = "" , 
       title = "Source of heterogeneity") + 
  theme_bw() +
  theme(legend.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black")
  ) 

# CV
p.CV <- ggplot(h_status, aes(levels, CVHs)) +
  geom_col(alpha = 1, color = wes_palette('GrandBudapest1', 4, 
                                          type = 'discrete')[3], 
           fill = wes_palette('GrandBudapest1', 4, type = 'discrete')[3]) +
  labs(y = expression(paste(italic(CV)[])), x = "" , 
       title = "Magnitude of heterogeneity") + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + # https://stackoverflow.com/questions/38722202/how-do-i-change-the-number-of-decimal-places-on-axis-labels-in-ggplot2
  theme_bw() +
  theme(legend.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black")
  )   


# M
p.M <- ggplot(h_status, aes(x = levels, y = Ms)) +
  geom_col(alpha = 1, color = wes_palette('GrandBudapest1', 4, 
                                          type = 'discrete')[4], fill = wes_palette('GrandBudapest1', 4, type = 'discrete')[4]) +
  labs(y = expression(paste(italic(M)[])), x = "" , 
       title = "Magnitude of heterogeneity") + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + # https://stackoverflow.com/questions/38722202/how-do-i-change-the-number-of-decimal-places-on-axis-labels-in-ggplot2
  theme_bw() +
  theme(legend.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black")
  )  


p.sigma + p.I2 + p.CV + p.M + plot_layout(ncol = 2, nrow = 2) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold"))


################################################################################
# Main model results

# TABLE
results_meta.model.off <- mod_results(meta.model.off, 
                                      mod = "1", 
                                      at = NULL, group = "StudyID")


# FIGURE: overall effect of yolk hormones on fitness
figure_meta.model.off <- orchaRd::orchard_plot(meta.model.off, mod = "1", 
                                               group = "StudyID", 
                                               xlab = "Effect size",
                                               transfm = "none",
                                               trunk.size = 2,
                                               branch.size = 2,
                                               twig.size = 1)


################################################################################
# 5. BIOLOGICAL HYPOTHESIS

# BH1 - Meta-regression - Biological Hypothesis 1 (BH.1)

# These 5 predictions (BH1.1 - BH1.5) are expected to occur at the offspring level.
# Hence I will run the models using only the data base for offspring fitness. 
# Following the pre-registration, I will 1) subset the data base for each hormone, 
# create a new phylogenetic tree and covariance for each data base (because the
# number of species changes across data bases), and then I will run each of the 
# models. 

################################################################################
# BH.1.1: Androgens and fitness 

# Prediction: High concentrations of egg androgens (i.e., androstenedione,
# testosterone and 5-α dihydrotestosterone) are positively linked with fitness.

# Data base: rows containing information on androgen hormones and offspring fitness.
################################################################################

# PHYLOGENETIC TREE - Species names:

# I subset the data base to have a new data base that contains information for
# offspring fitness and androgens.
meta.final_ok_ok_off_androgens <- droplevels(subset(meta.final_ok_ok_off, 
                                                    meta.final_ok_ok_off$Hormone_measured_general == "androgens"))

# # Create the phylogenetic tree
# resolved_names_off_androgens <- tnrs_match_names(as.character(unique(meta.final_ok_ok_off_androgens$Species)))
#
# # Saving the taxonomic data created on the 5th Feb 2025 to speed the
# # process in the future and allow full reproducibility
# save(resolved_names_off_androgens, file = "data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_off_androgens.RData")

# Loading the taxonomic data created on the 5th Feb 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_off_androgens.RData") #resolved_names_off_androgens

# # extracting phylogenetic information
# my_tree_off_androgens <- tol_induced_subtree(ott_ids =
#                                                resolved_names_off_androgens[,"ott_id"],
#                                              label_format = "name")
# # # Quick tree plotting
# # plot(my_tree_off_androgens, no.margin = TRUE)
# #
# # We need to check for the existence of polytomies
# is.binary(my_tree_off_androgens)
# # Yes, meaning there are no polytomies. Let's go on.
# #
# # To confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# my_tree_off_androgens$tip.label <- gsub("_", " ", my_tree_off_androgens$tip.label)
# 
# intersect(as.character(my_tree_off_androgens$tip.label),
#           as.character(meta.final_ok_ok_off_androgens$Species))
# 
# # Listed in our database but not in the tree
# setdiff(as.character(meta.final_ok_ok_off_androgens$Species),
#         as.character(my_tree_off_androgens$tip.label))
# 
# # Listed in the tree but not in our database
# setdiff(as.character(my_tree_off_androgens$tip.label),
#         as.character(meta.final_ok_ok_off_androgens$Species))
# # No error or inconsistencies found.
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# # 5th February 2025
# save(my_tree_off_androgens, file = "data/outputs/phylogenetic_files/tree_off_androgens.Rdata")

# # We can now load the saved tree
load("data/outputs/phylogenetic_files/tree_off_androgens.Rdata") #my_tree_off_androgens

# # Compute branch lengths of tree
# phylo_branch_off_androgens <- compute.brlen(my_tree_off_androgens,
#                                             method = "Grafen", power = 1)
# 
# # Check if tree is ultrametric
# is.ultrametric(phylo_branch_off_androgens)
# # TRUE
# 
# # Matrix to be included in the models
# phylo_cor_off_androgens <- vcv(phylo_branch_off_androgens, cor = T)
# #
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(phylo_cor_off_androgens, file = "data/outputs/phylogenetic_files/phylo_cor_off_androgens.Rdata")

# we can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor_off_androgens.Rdata") #phylo_cor

# Creating a duplicate species variable for the phylogenetic analysis
meta.final_ok_ok_off_androgens$Species_phylo_off_androgens <- 
  meta.final_ok_ok_off_androgens$Species



# # VARIANCE-COVARIANCE MATRIX:
VCV_ESVar_off_androgens <- matrix(0, nrow = nrow(meta.final_ok_ok_off_androgens), 
                                  ncol = nrow(meta.final_ok_ok_off_androgens))

# Names rows and columns for each obsID
rownames(VCV_ESVar_off_androgens) <- meta.final_ok_ok_off_androgens[, "EffectID"]
colnames(VCV_ESVar_off_androgens) <- meta.final_ok_ok_off_androgens[, "EffectID"]

# Finds effect sizes that come from the same study
shared_coord_off_androgens <- which(meta.final_ok_ok_off_androgens[, "StudyID"] %in% 
                                      meta.final_ok_ok_off_androgens[duplicated(meta.final_ok_ok_off_androgens[, "StudyID"]), 
                                                                     "StudyID"] == TRUE)

combinations_off_androgens <- do.call("rbind", tapply(shared_coord_off_androgens, 
                                                      meta.final_ok_ok_off_androgens[shared_coord_off_androgens, "StudyID"], 
                                                      function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations_off_androgens)[1]) {
  p1_off_androgens <- combinations_off_androgens[i, 1]
  p2_off_androgens <- combinations_off_androgens[i, 2]
  p1_p2_cov_off_androgens <- 0.5 * sqrt(meta.final_ok_ok_off_androgens[p1_off_androgens, "cor_var"]) * 
    sqrt(meta.final_ok_ok_off_androgens[p2_off_androgens, "cor_var"])
  VCV_ESVar_off_androgens[p1_off_androgens, p2_off_androgens] <- p1_p2_cov_off_androgens
  VCV_ESVar_off_androgens[p2_off_androgens, p1_off_androgens] <- p1_p2_cov_off_androgens
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(VCV_ESVar_off_androgens) <- meta.final_ok_ok_off_androgens[, "cor_var"]

# # In case you want to visually double check the matrix outside of R
# write.csv(VCV_ESVar_off_androgens, 'data/outputs/variance-covariance_matrices/VCV_ESVar_off_androgens.csv')



# STATISTICAL MODEL (Table S5):
meta.regression.bh1.1 <- rma.mv(cor,
                                VCV_ESVar_off_androgens,
                                mods = ~ 1,
                                random = list(~ 1 | StudyID,
                                              ~ 1 | LaboratoryID,
                                              ~ 1 | PopulationID,
                                              ~ 1 | Species,
                                              ~ 1 | Species_phylo_off_androgens,
                                              ~ 1 | EffectID),
                                method = "REML",
                                R = list(Species_phylo_off_androgens = phylo_cor_off_androgens),
                                test = "t",
                                data = meta.final_ok_ok_off_androgens)

#save(meta.regression.bh1.1, file = "data/outputs/statistical_models/meta_regression_bh1_1.RData")
load(file = "data/outputs/statistical_models/meta_regression_bh1_1.RData") #meta.regression.bh1.1

# Printing the summary results of the model
print(meta.regression.bh1.1, digits = 3)
predict(meta.regression.bh1.1, digits = 3)


################################################################################
# Heterogeneity
# Typical sampling variance
sigma2_v(meta.regression.bh1.1)
# Total unstandardized raw variance (i.e. total heterogeneity)
round(sum(meta.regression.bh1.1$sigma2),3)

# adding sigmas
h_status$sigma2s <- c(sum(meta.regression.bh1.1$sigma2),
                      meta.regression.bh1.1$sigma2[1],
                      meta.regression.bh1.1$sigma2[2],
                      meta.regression.bh1.1$sigma2[3],
                      meta.regression.bh1.1$sigma2[4],
                      meta.regression.bh1.1$sigma2[5],
                      meta.regression.bh1.1$sigma2[6])
round(h_status$sigma2s, 2)


# # I2, CV and M
# round(h.calc(meta.regression.bh1.1),2)
round(i2_ml(meta.regression.bh1.1),2)
round(cvh2_ml(meta.regression.bh1.1),2)
round(m2_ml(meta.regression.bh1.1),2)


# TABLE WITH RESULTS:
results_bh1.1 <- mod_results(meta.regression.bh1.1, 
                             mod = "1", 
                             at = NULL, group = "StudyID")

round(as.data.frame(results_bh1.1[[1]])[, c(2:6)], 3)

# FIGURE: overall effect of yolk androgens on fitness
figure_bh1.1 <- orchaRd::orchard_plot(meta.regression.bh1.1, 
                                      mod = "1", 
                                      group = "StudyID", 
                                      xlab = "Effect size",
                                      transfm = "none",
                                      trunk.size = 2,
                                      branch.size = 2,
                                      twig.size = 1)

# The effect size goes into an opposite direction (i.e., negative) than predicted (i.e., positive).
# This effect is ~0.05 and statistically not significant. 

################################################################################
# BH.1.2: corticosterone and fitness

# Prediction: High concentrations of glucocorticoids are negatively linked with 
# fitness.


# Data base: rows containing information on corticosterone and offspring fitness.

################################################################################
# I subset the data base
meta.final_ok_ok_off_cort <- droplevels(subset(meta.final_ok_ok_off, 
                                               meta.final_ok_ok_off$Hormone_measured_general == "corticosterone"))

# PHYLOGENETIC TREE:
# resolved_names_off_cort <- tnrs_match_names(as.character(unique(meta.final_ok_ok_off_cort$Species)))
#
# # Saving the taxonomic data created on the 5th February 2025 to speed the
# # process in the future and allow full reproducibility
# save(resolved_names_off_cort, file = "data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_off_cort.RData")

# Loading the taxonomic data created on the 5th February 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_off_cort.RData") #resolved_names_off_cort

# # extracting phylogenetic information
# my_tree_off_cort <- tol_induced_subtree(ott_ids =
#                                          resolved_names_off_cort[,"ott_id"],
#                                        label_format = "name")
# # # Quick tree plotting
# # plot(my_tree_off_cort, no.margin = TRUE)
# 
# # We need to check for the existence of polytomies
# is.binary(my_tree_off_cort)
# # Yes, meaning there are no polytomies. Let's go on.
# 
# # To confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# my_tree_off_cort$tip.label <- gsub("_", " ", my_tree_off_cort$tip.label)
# 
# intersect(as.character(my_tree_off_cort$tip.label),
#         as.character(meta.final_ok_ok_off_cort$Species))
# 
# # Listed in our database but not in the tree
# setdiff(as.character(meta.final_ok_ok_off_cort$Species),
#        as.character(my_tree_off_cort$tip.label))
# 
# # Listed in the tree but not in our database
# setdiff(as.character(my_tree_off_cort$tip.label),
#       as.character(meta.final_ok_ok_off_cort$Species))
# # No error or inconsistencies found.
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# # 5th February 2025
# save(my_tree_off_cort, file = "data/outputs/phylogenetic_files/tree_off_cort.Rdata")

# We can now load the saved tree
load("data/outputs/phylogenetic_files/tree_off_cort.Rdata") #my_tree_off_cort

# # Compute branch lengths of tree
# phylo_branch_off_cort <- compute.brlen(my_tree_off_cort, method = "Grafen", power = 1)
# 
# # Check if tree is ultrametric
# is.ultrametric(phylo_branch_off_cort)
# # TRUE
# 
# # Matrix to be included in the models
# phylo_cor_off_cort <- vcv(phylo_branch_off_cort, cor = T)
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(phylo_cor_off_cort, file = "data/outputs/phylogenetic_files/phylo_cor_off_cort.Rdata")

# We can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor_off_cort.Rdata") #phylo_cor_off_cort

# Creating a duplicate species variable for the phylogenetic analysis
meta.final_ok_ok_off_cort$Species_phylo_off_cort <- meta.final_ok_ok_off_cort$Species



# VARIANCE-COVARIANCE MATRIX:
VCV_ESVar_off_cort <- matrix(0, nrow = nrow(meta.final_ok_ok_off_cort), 
                             ncol = nrow(meta.final_ok_ok_off_cort))

# Names rows and columns for each obsID
rownames(VCV_ESVar_off_cort) <- meta.final_ok_ok_off_cort[, "EffectID"]
colnames(VCV_ESVar_off_cort) <- meta.final_ok_ok_off_cort[, "EffectID"]

# Finds effect sizes that come from the same study
shared_coord_off_cort <- which(meta.final_ok_ok_off_cort[, "StudyID"] %in% 
                                 meta.final_ok_ok_off_cort[duplicated(meta.final_ok_ok_off_cort[, "StudyID"]), 
                                                           "StudyID"] == TRUE)

combinations_off_cort <- do.call("rbind", tapply(shared_coord_off_cort, 
                                                 meta.final_ok_ok_off_cort[shared_coord_off_cort, "StudyID"], 
                                                 function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations_off_cort)[1]) {
  p1_off_cort <- combinations_off_cort[i, 1]
  p2_off_cort <- combinations_off_cort[i, 2]
  p1_p2_cov_off_cort <- 0.5 * sqrt(meta.final_ok_ok_off_cort[p1_off_cort, "cor_var"]) * 
    sqrt(meta.final_ok_ok_off_cort[p2_off_cort, "cor_var"])
  VCV_ESVar_off_cort[p1_off_cort, p2_off_cort] <- p1_p2_cov_off_cort
  VCV_ESVar_off_cort[p2_off_cort, p1_off_cort] <- p1_p2_cov_off_cort
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(VCV_ESVar_off_cort) <- meta.final_ok_ok_off_cort[, "cor_var"]

# # In case you want to visually double check the matrix outside of R
# write.csv(VCV_ESVar_off_cort, 'data/outputs/variance-covariance_matrices/VCV_ESVar_off_cort.csv')


# STATISTICAL MODEL  (Table S5):
meta.regression.bh1.2 <- rma.mv(cor,
                                VCV_ESVar_off_cort,
                                mods = ~ 1,
                                random = list(~ 1 | StudyID,
                                              ~ 1 | LaboratoryID,
                                              ~ 1 | PopulationID,
                                              ~ 1 | Species,
                                              ~ 1 | Species_phylo_off_cort,
                                              ~ 1 | EffectID),
                                method = "REML",
                                R = list(Species_phylo_off_cort = phylo_cor_off_cort),
                                test = "t",
                                data = meta.final_ok_ok_off_cort)

#save(meta.regression.bh1.2, file = "data/outputs/statistical_models/meta_regression_bh1_2.RData")
load(file = "data/outputs/statistical_models/meta_regression_bh1_2.RData") #meta.regression.bh1.2

# Printing the summary results of the model
print(meta.regression.bh1.2, digits = 3)
predict(meta.regression.bh1.2, digits = 3)


################################################################################
# Heterogeneity
# Typical sampling variance
sigma2_v(meta.regression.bh1.2)
# Total unstandardized raw variance (i.e. total heterogeneity)
round(sum(meta.regression.bh1.2$sigma2),3)


# adding sigmas
h_status$sigma2s <- c(sum(meta.regression.bh1.2$sigma2),
                      meta.regression.bh1.2$sigma2[1],
                      meta.regression.bh1.2$sigma2[2],
                      meta.regression.bh1.2$sigma2[3],
                      meta.regression.bh1.2$sigma2[4],
                      meta.regression.bh1.2$sigma2[5],
                      meta.regression.bh1.2$sigma2[6])
round(h_status$sigma2s, 2)



# # I2, CV and M
# round(h.calc(meta.regression.bh1.2),2)
round(i2_ml(meta.regression.bh1.2),2)
round(cvh2_ml(meta.regression.bh1.2),2)
round(m2_ml(meta.regression.bh1.2),2)


# TABLE WITH RESULTS:
results_bh1.2 <- mod_results(meta.regression.bh1.2, 
                             mod = "1", 
                             at = NULL, group = "StudyID")


# FIGURE: overall effect of yolk cort on fitness
figure_bh1.2 <- orchaRd::orchard_plot(meta.regression.bh1.2, 
                                      mod = "1", 
                                      group = "StudyID", 
                                      xlab = "Effect size",
                                      transfm = "none",
                                      trunk.size = 2,
                                      branch.size = 2,
                                      twig.size = 1)

# The effect size goes into the opposite direction than the one we predicted.
# This effect is quite small (0.076) and statistically not significant. 

################################################################################
# BH.1.3: THs and fitness

# BH.1.3: Thyroid hormones are positively linked with fitness.

# Data base: rows containing information on THs hormones and offspring fitness.

################################################################################
# I subset the data base
meta.final_ok_ok_off_ths <- droplevels(subset(meta.final_ok_ok_off, 
                                              meta.final_ok_ok_off$Hormone_measured_general == "TH3 and TH4"))


# PYLOGENETIC TREE:
# resolved_names_off_ths <- tnrs_match_names(as.character(unique(meta.final_ok_ok_off_ths$Species)))
# #
# # # Saving the taxonomic data created on the 5th February 2025 to speed the
# # # process in the future and allow full reproducibility
# save(resolved_names_off_ths, file = "data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_off_ths.RData")

# Loading the taxonomic data created on the 5th February 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_off_ths.RData") #resolved_names_off_ths

# # extracting phylogenetic information
# my_tree_off_ths <- tol_induced_subtree(ott_ids =
#                                          resolved_names_off_ths[,"ott_id"],
#                                        label_format = "name")
# # # Quick tree plotting
# # plot(my_tree_off_ths, no.margin = TRUE)
# 
# # We need to check for the existence of polytomies
# is.binary(my_tree_off_ths)
# # Yes, meaning there are no polytomies. Let's go on.
# #
# # To confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# my_tree_off_ths$tip.label <- gsub("_", " ", my_tree_off_ths$tip.label)
# 
# 
# # Listed in our database but not in the tree
# setdiff(as.character(meta.final_ok_ok_off_ths$Species),
#         as.character(my_tree_off_ths$tip.label))
# 
# # Listed in the tree but not in our database
# setdiff(as.character(my_tree_off_ths$tip.label),
#         as.character(meta.final_ok_ok_off_ths$Species))
# # No error or inconsistencies found.
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# # 5th February 2025
# save(my_tree_off_ths, file = "data/outputs/phylogenetic_files/tree_off_ths.Rdata")

# We can now load the saved tree
load("data/outputs/phylogenetic_files/tree_off_ths.Rdata") #my_tree_off_ths

# # Compute branch lengths of tree
# phylo_branch_off_ths <- compute.brlen(my_tree_off_ths, method = "Grafen", power = 1)
# 
# # Check if tree is ultrametric
# is.ultrametric(phylo_branch_off_ths)
# # TRUE
# 
# # Matrix to be included in the models
# phylo_cor_off_ths <- vcv(phylo_branch_off_ths, cor = T)
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(phylo_cor_off_ths, file = "data/outputs/phylogenetic_files/phylo_cor_off_ths.Rdata")

# we can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor_off_ths.Rdata") #phylo_branch_off_ths

# Creating a duplicate species variable for the phylogenetic analysis
meta.final_ok_ok_off_ths$Species_phylo_off_ths <- meta.final_ok_ok_off_ths$Species

# Matrix for OFFSPRING
VCV_ESVar_off_ths <- matrix(0, nrow = nrow(meta.final_ok_ok_off_ths), 
                            ncol = nrow(meta.final_ok_ok_off_ths))

# Names rows and columns for each obsID
rownames(VCV_ESVar_off_ths) <- meta.final_ok_ok_off_ths[, "EffectID"]
colnames(VCV_ESVar_off_ths) <- meta.final_ok_ok_off_ths[, "EffectID"]

# Finds effect sizes that come from the same study
shared_coord_off_ths <- which(meta.final_ok_ok_off_ths[, "StudyID"] %in% 
                                meta.final_ok_ok_off_ths[duplicated(meta.final_ok_ok_off_ths[, "StudyID"]), 
                                                         "StudyID"] == TRUE)

combinations_off_ths <- do.call("rbind", tapply(shared_coord_off_ths, 
                                                meta.final_ok_ok_off_ths[shared_coord_off_ths, "StudyID"], 
                                                function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations_off_ths)[1]) {
  p1_off_ths <- combinations_off_ths[i, 1]
  p2_off_ths <- combinations_off_ths[i, 2]
  p1_p2_cov_off_ths <- 0.5 * sqrt(meta.final_ok_ok_off_ths[p1_off_ths, "cor_var"]) * 
    sqrt(meta.final_ok_ok_off_ths[p2_off_ths, "cor_var"])
  VCV_ESVar_off_ths[p1_off_ths, p2_off_ths] <- p1_p2_cov_off_ths
  VCV_ESVar_off_ths[p2_off_ths, p1_off_ths] <- p1_p2_cov_off_ths
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(VCV_ESVar_off_ths) <- meta.final_ok_ok_off_ths[, "cor_var"]

# # In case you want to visually double check the matrix outside of R
# write.csv(VCV_ESVar_off_ths, 'data/outputs/variance-covariance_matrices/VCV_ESVar_off_ths.csv')



# STATISTICAL MODEL (Table S5):
meta.regression.bh1.3 <- rma.mv(cor,
                                VCV_ESVar_off_ths,
                                mods = ~ 1,
                                random = list(~ 1 | StudyID,
                                              ~ 1 | LaboratoryID,
                                              ~ 1 | PopulationID,
                                              ~ 1 | Species,
                                              ~ 1 | Species_phylo_off_ths,
                                              ~ 1 | EffectID),
                                method = "REML",
                                R = list(Species_phylo_off_ths = phylo_cor_off_ths),
                                test = "t",
                                data = meta.final_ok_ok_off_ths)

#save(meta.regression.bh1.3, file = "data/outputs/statistical_models/meta_regression_bh1_3.RData")
load(file = "data/outputs/statistical_models/meta_regression_bh1_3.RData") #meta.regression.bh1.3

# Printing the summary results of the model
print(meta.regression.bh1.3, digits = 3)
predict(meta.regression.bh1.3, digits = 3)


################################################################################
# Heterogeneity
# Typical sampling variance
sigma2_v(meta.regression.bh1.3)
# Total unstandardized raw variance (i.e. total heterogeneity)
round(sum(meta.regression.bh1.3$sigma2),3)


# adding sigmas
h_status$sigma2s <- c(sum(meta.regression.bh1.3$sigma2),
                      meta.regression.bh1.3$sigma2[1],
                      meta.regression.bh1.3$sigma2[2],
                      meta.regression.bh1.3$sigma2[3],
                      meta.regression.bh1.3$sigma2[4],
                      meta.regression.bh1.3$sigma2[5],
                      meta.regression.bh1.3$sigma2[6])
round(h_status$sigma2s, 2)


# # I2, CV and M
# round(h.calc(meta.regression.bh1.3),2)
round(i2_ml(meta.regression.bh1.3),2)
round(cvh2_ml(meta.regression.bh1.3),2)
round(m2_ml(meta.regression.bh1.3),2)


# TABLE WITH RESULTS:
results_bh1.3 <- mod_results(meta.regression.bh1.3, 
                             mod = "1", 
                             at = NULL, group = "StudyID")



# FIGURE: overall effect of yolk ths on fitness
figure_bh1.3 <- orchaRd::orchard_plot(meta.regression.bh1.3, 
                                      mod = "1", 
                                      group = "StudyID", 
                                      xlab = "Effect size",
                                      transfm = "none",
                                      trunk.size = 2,
                                      branch.size = 2,
                                      twig.size = 1)


# The effect size goes into the same direction as the one we predicted. We found 
# a small positive effect (0.02) that is statistically not significant. 


################################################################################
# BH.1.4 & BH.1.5: These hypothesis were in relation to progesterone and estradiol.

# We planned to explore their association with fitness without any prediction on 
# the directionality due to the few studies that have measured these hormones.
# However, in our final data base, we did not have any studies looking at the 
# relationship between progesterone and estradiol with fitness. Hence, we cannot
# test this.
################################################################################


################################################################################
# BH2 - Sex-specific fitness consequences

# We predicted males to have an average larger effect size than females. This 
# effect is predicted to be independent from the hormone and trait measured. The
# direction of the effect is predicted to be in the same direction as stated for 
# BH.1.

# Data base: Predictions were done for offspring fitness. Hence, I will use the
# data base that includes information only for the offspring, and within this
# data base, I will keep those rows that contain information on the sex of the
# offspring.

# Before going into the statistical model, I need to first create a new 
# phylogenetic tree and a covariance matrix. This is because from 400 rows, 
# we only have 49 rows with data about the sex of the offspring. 

################################################################################

# I subset the data base so I work only with those rows for which I do have info
# on the sex of the offspring.
meta.final_ok_ok_off_sex <- droplevels(subset(meta.final_ok_ok_off, 
                                              meta.final_ok_ok_off$Off_sex != "NA"))

# # PHYLOGENETIC TREE:
# resolved_names_off_sex <- tnrs_match_names(as.character(unique(meta.final_ok_ok_off_sex$Species)))
# 
# # Saving the taxonomic data created on the 5th February 2025 to speed the
# # process in the future and allow full reproducibility
# save(resolved_names_off_sex, file = "data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_off_sex.RData")

# Loading the taxonomic data created on the 5th February 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_off_sex.RData") #resolved_names_off_sex

# # extracting phylogenetic information
# my_tree_off_sex <- tol_induced_subtree(ott_ids =
#                                          resolved_names_off_sex[,"ott_id"],
#                                        label_format = "name")
# # # Quick tree plotting
# # plot(my_tree_off_sex, no.margin = TRUE)
# 
# # We need to check for the existence of polytomies
# is.binary(my_tree_off_sex)
# # Yes, meaning there are no polytomies. Let's go on.
# 
# # To confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# my_tree_off_sex$tip.label <- gsub("_", " ", my_tree_off_sex$tip.label)
# 
# intersect(as.character(my_tree_off_sex$tip.label),
#           as.character(meta.final_ok_ok_off_sex$Species))
# 
# # Listed in our database but not in the tree
# setdiff(as.character(meta.final_ok_ok_off_sex$Species),
#         as.character(my_tree_off_sex$tip.label))
# 
# # Listed in the tree but not in our database
# setdiff(as.character(my_tree_off_sex$tip.label),
#         as.character(meta.final_ok_ok_off_sex$Species))
# # No error or inconsistencies found.
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# # 5th February 2025
# save(my_tree_off_sex, file = "data/outputs/phylogenetic_files/tree_off_sex.Rdata")

# We can now load the saved tree
load("data/outputs/phylogenetic_files/tree_off_sex.Rdata") #my_tree_off_sex

# # Compute branch lengths of tree
# phylo_branch_off_sex <- compute.brlen(my_tree_off_sex, method = "Grafen", power = 1)
# 
# # Check if tree is ultrametric
# is.ultrametric(phylo_branch_off_sex)
# # TRUE
# 
# # Matrix to be included in the models
# phylo_cor_off_sex <- vcv(phylo_branch_off_sex, cor = T)
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(phylo_cor_off_sex, file = "data/outputs/phylogenetic_files/phylo_cor_off_sex.Rdata")

# we can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor_off_sex.Rdata") #phylo_cor_off_sex

# Creating a duplicate species variable for the phylogenetic analysis
meta.final_ok_ok_off_sex$Species_phylo_off_sex <- meta.final_ok_ok_off_sex$Species



# VARIANCE-COVARIANCE MATRIX:
VCV_ESVar_off_sex <- matrix(0, nrow = nrow(meta.final_ok_ok_off_sex), 
                            ncol = nrow(meta.final_ok_ok_off_sex))

# Names rows and columns for each obsID
rownames(VCV_ESVar_off_sex) <- meta.final_ok_ok_off_sex[, "EffectID"]
colnames(VCV_ESVar_off_sex) <- meta.final_ok_ok_off_sex[, "EffectID"]

# Finds effect sizes that come from the same study
shared_coord_off_sex <- which(meta.final_ok_ok_off_sex[, "StudyID"] %in% 
                                meta.final_ok_ok_off_sex[duplicated(meta.final_ok_ok_off_sex[, "StudyID"]), 
                                                         "StudyID"] == TRUE)

combinations_off_sex <- do.call("rbind", tapply(shared_coord_off_sex, 
                                                meta.final_ok_ok_off_sex[shared_coord_off_sex, "StudyID"], 
                                                function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations_off_sex)[1]) {
  p1_off_sex <- combinations_off_sex[i, 1]
  p2_off_sex <- combinations_off_sex[i, 2]
  p1_p2_cov_off_sex <- 0.5 * sqrt(meta.final_ok_ok_off_sex[p1_off_sex, "cor_var"]) * 
    sqrt(meta.final_ok_ok_off_sex[p2_off_sex, "cor_var"])
  VCV_ESVar_off_sex[p1_off_sex, p2_off_sex] <- p1_p2_cov_off_sex
  VCV_ESVar_off_sex[p2_off_sex, p1_off_sex] <- p1_p2_cov_off_sex
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(VCV_ESVar_off_sex) <- meta.final_ok_ok_off_sex[, "cor_var"]

# # In case you want to visually double check the matrix outside of R
# write.csv(VCV_ESVar_off_sex, 'data/outputs/variance-covariance_matrices/VCV_ESVar_off_sex.csv')

nrow(meta.final_ok_ok_off_sex)


# STATISTICAL MODEL:
# The new data set contains 45 rows that have hormones that are predicted to have
# a positive effect on offspring fitness traits (i.e., there are no papers
# testing the effect of corticosterone). Hence the variable "Hormone_measured_general" 
# has only one level and interactions (which are needed to test our pre-registered
# prediction) cannot be estimated. We can only test if the effect size of 
# different hormones on fitness traits differ for males and females.

# The model that we can run based on our prediction and the final sample size
# The "-1" nest to the moderator removes the intercept so the results are easier
# to interpret. 

meta.regression.bh2 <- rma.mv(cor,
                              VCV_ESVar_off_sex,
                              mods = ~ (Hormone_measured_general + Off_sex +
                                          Hormone_measured_general * Off_sex) - 1,
                              random = list(~ 1 | StudyID,
                                            ~ 1 | LaboratoryID,
                                            ~ 1 | PopulationID,
                                            ~ 1 | Species,
                                            ~ 1 | Species_phylo_off_sex,
                                            ~ 1 | EffectID),
                              method = "REML",
                              R = list(Species_phylo_off_sex = phylo_cor_off_sex),
                              test = "t",
                              data = meta.final_ok_ok_off_sex)

#save(meta.regression.bh2, file = "data/outputs/statistical_models/meta_regression_bh2.RData")
load(file = "data/outputs/statistical_models/meta_regression_bh2.RData") #meta.regression.bh2

# We had a Convergence Problem with the rma.mv() Function. 
# By using the command "verbose = TRUE" we obtain the information that starting 
# around the iteration 297 there is no further change in the log likelihood when 
# rounded to 4 decimal placed. 
# However, because the default settings of the model had a very small threshold, 
# the model could not for determined when convergence occurred.
# We therefore modified this threshold with the commands "control = list(...)".
# Results remain the same and we no longer have a error message.

# Printing the summary results of the model
print(meta.regression.bh2, digits = 3)
# predict(meta.regression.bh2, digits = 3)


# Calculate marginal R2 with r2_ml
R2.method_bh2 <- r2_ml(meta.regression.bh2)
round(R2.method_bh2 * 100, 1)


# In order to better understand and plot the results, I will create an artificial
# variable that includes this interaction.
meta.final_ok_ok_off_sex$hormone_sex <- paste(meta.final_ok_ok_off_sex$Hormone_measured_general,
                                              meta.final_ok_ok_off_sex$Off_sex,
                                              sep = "_")


meta.regression.bh2_artificial <- rma.mv(cor,
                                         VCV_ESVar_off_sex,
                                         mods = ~ hormone_sex - 1,
                                         random = list(~ 1 | StudyID,
                                                       ~ 1 | LaboratoryID,
                                                       ~ 1 | PopulationID,
                                                       ~ 1 | Species,
                                                       ~ 1 | Species_phylo_off_sex,
                                                       ~ 1 | EffectID),
                                         method = "REML",
                                         #verbose = TRUE,
                                         control = list(rel.tol = 1e-8),
                                         R = list(Species_phylo_off_sex = phylo_cor_off_sex),
                                         test = "t",
                                         data = meta.final_ok_ok_off_sex)

#save(meta.regression.bh2_artificial, file = "data/outputs/statistical_models/meta_regression_bh2_artificial.RData")
load(file = "data/outputs/statistical_models/meta_regression_bh2_artificial.RData") #meta.regression.bh2_artificial

# Printing the summary results of the model
print(meta.regression.bh2_artificial, digits = 3)

# Calculate marginal R2 with r2_ml
R2.method_bh2_artificial <- r2_ml(meta.regression.bh2_artificial)
round(R2.method_bh2_artificial * 100, 1)

# TABLE WITH RESULTS:
results_bh2 <- orchaRd::mod_results(meta.regression.bh2_artificial, 
                                    mod = "hormone_sex", 
                                    group = "StudyID", 
                                    subset = TRUE)

round(as.data.frame(results_bh2[[1]])[, c(2:6)], 3)


# FIGURE: effect of each group of hormones on fitness traits.
fig_sex_differences <- orchaRd::orchard_plot(meta.regression.bh2_artificial, 
                                             mod = "hormone_sex", 
                                             group = "StudyID", 
                                             xlab = "Effect size",
                                             trunk.size = 2,
                                             branch.size = 2,
                                             twig.size = 1,
                                             colour = FALSE)


# For both male and female offspring, egg hormones have a negative effect on their
# fitness traits. Overall, this effect tends to be bigger for androgens than thyroid hormones.
# The effect of androgens on male offspring seems to be bigger than for females, 
# whereas the effect of THs on male offspring seems to be smaller than for females. 


# Pair-wise comparisons including post-hoc Wald tests
bh2_artificial_pc <- rma.mv(cor,
                            VCV_ESVar_off_sex,
                            mods = ~ hormone_sex,
                            random = list(~ 1 | StudyID,
                                          ~ 1 | LaboratoryID,
                                          ~ 1 | PopulationID,
                                          ~ 1 | Species,
                                          ~ 1 | Species_phylo_off_sex,
                                          ~ 1 | EffectID),
                            method = "REML",
                            R = list(Species_phylo_off_sex = phylo_cor_off_sex),
                            test = "t",
                            data = meta.final_ok_ok_off_sex)

#save(bh2_artificial_pc, file = "data/outputs/statistical_models/bh2_artificial_pc.RData")
load(file = "data/outputs/statistical_models/bh2_artificial_pc.RData") #bh2_artificial_pc

print(bh2_artificial_pc, digits = 3)


# I now want to test for statistical differences between groups. Male and fem
# information for androgens can be obtained from the outcome of the table. For thyroid
# hormones I need to do Post-hoc Wald tests: 
car::linearHypothesis(bh2_artificial_pc, rbind(c(0,0,1,-1)))


################################################################################
# BH3 - Relationship between maternal egg hormones and offspring fitness proxies
# vary depending on age of the offspring.
################################################################################

# BH3.1. - We predicted that the relationships between egg hormones and offspring
# fitness proxies become weaker from hatching day to before (e.g. nesting period 
# other than hatching date) and after offspring independence (e.g. adulthood)

# Data base: predictions were done for offspring fitness. Hence, I will use the 
# data base that includes information only for the offspring.

################################################################################

# I first check if we have enough data points for each hormone and life history
# stage to be used in the statistical model:
table(meta.final_ok_ok_off$Off_life_stage,
      meta.final_ok_ok_off$Hormone_measured_general)

# Levels of offspring life history stage: hatching, before independence, and 
# after independence. 

# For corticosterone and THs we do not have enough data points to test 
# this prediction. In particular, for corticosterone we have less than 5 point at 
# hatching and after independence (2 and 1 data points, respectively). For THs
# we do not have any information for the effect of these hormones at hatching and
# after independence.

# Therefore, we can only test the effect of androgens on fitness over different
# life history stages. Therefore, I will only work with the data base that includes
# information for androgens.


# I will re-order the levels so they are ordered chronologically.
meta.final_ok_ok_off_androgens$Off_life_stage <- factor(meta.final_ok_ok_off_androgens$Off_life_stage,
                                                        levels = c("hatching", 
                                                                   "before independence", 
                                                                   "after independence"))

# STATISTICAL MODEL:
meta.regression.bh3.1 <- rma.mv(cor,
                                VCV_ESVar_off_androgens,
                                mods = ~ Off_life_stage - 1,
                                random = list(~ 1 | StudyID,
                                              ~ 1 | LaboratoryID,
                                              ~ 1 | PopulationID,
                                              ~ 1 | Species,
                                              ~ 1 | Species_phylo_off_androgens,
                                              ~ 1 | EffectID),
                                method = "REML",
                                R = list(Species_phylo_off_androgens = phylo_cor_off_androgens),
                                test = "t",
                                data = meta.final_ok_ok_off_androgens)

#save(meta.regression.bh3.1, file = "data/outputs/statistical_models/meta_regression_bh3_1.RData")
load(file = "data/outputs/statistical_models/meta_regression_bh3_1.RData") #meta.regression.bh3.1

# Printing the summary results of the model
print(meta.regression.bh3.1, digits = 3)
# predict(meta.regression.bh3.1, digits = 3)

# Calculate marginal R2 with r2_ml
R2.method_bh3.1 <- r2_ml(meta.regression.bh3.1)
round(R2.method_bh3.1 * 100, 1)


# TABLE WITH RESULTS:
results_bh3.1 <- orchaRd::mod_results(meta.regression.bh3.1, 
                                      mod = "Off_life_stage", 
                                      group = "StudyID", 
                                      subset = TRUE)

round(as.data.frame(results_bh3.1[[1]])[, c(2:6)], 3)


# FIGURE: effect of each group of hormones on fitness traits.
fig_bh3.1 <- orchaRd::orchard_plot(meta.regression.bh3.1, 
                                   mod = "Off_life_stage", 
                                   group = "StudyID", 
                                   xlab = "Effect size",
                                   trunk.size = 2,
                                   branch.size = 2,
                                   twig.size = 1)

# Androgen hormones have a positive effect on fitness traits only at hatching. 
# Later on, the effect is negative. Effect sizes increase over time 
# (hatching < before independence < after independence).
# None of the effects is statistically significant.



# Pair-wise comparisons including post-hoc Wald tests
meta.regression.bh3.1_pc <- rma.mv(cor,
                                   VCV_ESVar_off_androgens,
                                   mods = ~ Off_life_stage,
                                   random = list(~ 1 | StudyID,
                                                 ~ 1 | LaboratoryID,
                                                 ~ 1 | PopulationID,
                                                 ~ 1 | Species,
                                                 ~ 1 | Species_phylo_off_androgens,
                                                 ~ 1 | EffectID),
                                   method = "REML",
                                   R = list(Species_phylo_off_androgens = phylo_cor_off_androgens),
                                   test = "t",
                                   data = meta.final_ok_ok_off_androgens)

#save(meta.regression.bh3.1_pc, file = "data/outputs/statistical_models/meta_regression_bh3_1_pc.RData")
load(file = "data/outputs/statistical_models/meta_regression_bh3_1_pc.RData") #meta.regression.bh3.1_pc

print(meta.regression.bh3.1_pc, digits = 3)


# # I now want to test for statistical differences between groups. Information 
# between hatching and before independence and hatching and after independence
# can be obtained from the outcome of the table. For comparisons between 
# before and after hatching I need to run a Post-hoc Wald test: 
car::linearHypothesis(meta.regression.bh3.1_pc, rbind(c(0,1,-1)))


###############################################################################
# BH.3.2 - For altricial species, we predict that the relationships between egg 
# hormones and offspring fitness proxies become weaker throughout the nesting phase.

# Data base: I have to work with the data base that only includes offspring fitness.
# And from this data base subset those rows that contain information on 
# nestling relative age (i.e., offspring age/average nesting time)
################################################################################

meta.final_ok_ok_relative_age <- droplevels((subset(meta.final_ok_ok_off,
                                                    meta.final_ok_ok_off$Off_relative_age != "NA")))

nrow(meta.final_ok_ok_relative_age)

# I check if we have enough data points for each hormone:
table(meta.final_ok_ok_relative_age$Hormone_measured_general)
# yes

table(meta.final_ok_ok_relative_age$Hormone_measured_general, meta.final_ok_ok_relative_age$StudyID)

# # PHYLOGENETIC TREE:
# resolved_names_relative_age <- tnrs_match_names(as.character(unique(meta.final_ok_ok_relative_age$Species)))
# 
# # Saving the taxonomic data created on the 5th February 2025 to speed the
# # process in the future and allow full reproducibility
# save(resolved_names_relative_age, file = "data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_relative_age.RData")

# Loading the taxonomic data created on the 5th February 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_relative_age.RData") #resolved_names_relative_age

# # extracting phylogenetic information
# my_tree_relative_age <- tol_induced_subtree(ott_ids =
#                                               resolved_names_relative_age[,"ott_id"],
#                                             label_format = "name")
# # # Quick tree plotting
# # plot(my_tree_relative_age, no.margin = TRUE)
# 
# # We need to check for the existence of polytomies
# is.binary(my_tree_relative_age)
# # Yes, meaning there are no polytomies. Let's go on.
# 
# # To confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# my_tree_relative_age$tip.label <- gsub("_", " ", my_tree_relative_age$tip.label)
# 
# intersect(as.character(my_tree_relative_age$tip.label),
#           as.character(meta.final_ok_ok_relative_age$Species))
# 
# # Listed in our database but not in the tree
# setdiff(as.character(meta.final_ok_ok_relative_age$Species),
#         as.character(my_tree_relative_age$tip.label))
# 
# # Listed in the tree but not in our database
# setdiff(as.character(my_tree_relative_age$tip.label),
#         as.character(meta.final_ok_ok_relative_age$Species))
# # No error or inconsistencies found.
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# # 5th February 2025
# save(my_tree_relative_age, file = "data/outputs/phylogenetic_files/tree_relative_age.Rdata")

# We can now load the saved tree
load("data/outputs/phylogenetic_files/tree_relative_age.Rdata") #my_tree_relative_age

# # Compute branch lengths of tree
# phylo_branch_relative_age <- compute.brlen(my_tree_relative_age, method = "Grafen", power = 1)
# #
# # Check if tree is ultrametric
# is.ultrametric(phylo_branch_relative_age)
# # TRUE
# 
# # Matrix to be included in the models
# phylo_cor_relative_age <- vcv(phylo_branch_relative_age, cor = T)
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(phylo_cor_relative_age, file = "data/outputs/phylogenetic_files/phylo_cor_relative_age.Rdata")

# we can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor_relative_age.Rdata") #phylo_cor_relative_age

# Creating a duplicate species variable for the phylogenetic analysis
meta.final_ok_ok_relative_age$Species_phylo_relative_age <- meta.final_ok_ok_relative_age$Species



# VARIANCE-COVARIANCE MATRIX:
VCV_ESVar_relative_age <- matrix(0, nrow = nrow(meta.final_ok_ok_relative_age), 
                                 ncol = nrow(meta.final_ok_ok_relative_age))

# Names rows and columns for each obsID
rownames(VCV_ESVar_relative_age) <- meta.final_ok_ok_relative_age[, "EffectID"]
colnames(VCV_ESVar_relative_age) <- meta.final_ok_ok_relative_age[, "EffectID"]

# Finds effect sizes that come from the same study
shared_coord_relative_age <- which(meta.final_ok_ok_relative_age[, "StudyID"] %in% 
                                     meta.final_ok_ok_relative_age[duplicated(meta.final_ok_ok_relative_age[, "StudyID"]), 
                                                                   "StudyID"] == TRUE)

combinations_relative_age <- do.call("rbind", tapply(shared_coord_relative_age, 
                                                     meta.final_ok_ok_relative_age[shared_coord_relative_age, "StudyID"], 
                                                     function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations_relative_age)[1]) {
  p1_relative_age <- combinations_relative_age[i, 1]
  p2_relative_age <- combinations_relative_age[i, 2]
  p1_p2_cov_relative_age <- 0.5 * sqrt(meta.final_ok_ok_relative_age[p1_relative_age, "cor_var"]) * 
    sqrt(meta.final_ok_ok_relative_age[p2_relative_age, "cor_var"])
  VCV_ESVar_relative_age[p1_relative_age, p2_relative_age] <- p1_p2_cov_relative_age
  VCV_ESVar_relative_age[p2_relative_age, p1_relative_age] <- p1_p2_cov_relative_age
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(VCV_ESVar_relative_age) <- meta.final_ok_ok_relative_age[, "cor_var"]

# In case you want to visually double check the matrix outside of R
# write.csv(VCV_ESVar_relative_age, 'data/outputs/variance-covariance_matrices/VCV_ESVar_relative_age.csv')


# STATISTICAL ANALYSIS:
meta.regression.bh3.2 <- rma.mv(cor,
                                VCV_ESVar_relative_age,
                                mods = ~ Hormone_measured_general * Off_relative_age,
                                random = list(~ 1 | StudyID,
                                              ~ 1 | LaboratoryID,
                                              ~ 1 | PopulationID,
                                              ~ 1 | Species,
                                              ~ 1 | Species_phylo_relative_age,
                                              ~ 1 | EffectID),
                                method = "REML",
                                R = list(Species_phylo_relative_age = phylo_cor_relative_age),
                                test = "t",
                                data = meta.final_ok_ok_relative_age)

#save(meta.regression.bh3.2, file = "data/outputs/statistical_models/meta_regression_bh3_2.RData")
load(file = "data/outputs/statistical_models/meta_regression_bh3_2.RData") #meta.regression.bh3.2

# Printing the summary results of the model
print(meta.regression.bh3.2, digits = 3)

# Calculate marginal R2 with r2_ml
R2.method_bh3.2 <- r2_ml(meta.regression.bh3.2)
round(R2.method_bh3.2 * 100, 1)

# TABLE WITH RESULTS FOR A CONTINUOUS VARIABLE:
results_bh3.2 <- orchaRd::mod_results(meta.regression.bh3.2,
                                      mod = "Off_relative_age", 
                                      group = "StudyID",
                                      weights = "prop",
                                      by = "Hormone_measured_general")

# FIGURE WITH A CONTINOUS VARIABLE:
figure_bh3.2 <- orchaRd::bubble_plot(results_bh3.2,
                                     group = "StudyID",
                                     mod = "Off_relative_age",
                                     xlab = "Off_relative_age",
                                     legend.pos = "bottom.right")

# Androgens and THs effect sizes decrease over offspring relative age, whereas
# glucocorticoid effect size increases. 

# # exploring results further
# results_bh3.2.alternative <- orchaRd::mod_results(meta.regression.bh3.2,
#                                                   mod = "Hormone_measured_general",
#                                                   group = "StudyID",
#                                                   weights = "prop",
#                                                   by = "Off_relative_age")
# 
# # I believe this is the mean effect for each hormone separately after accounting
# # for (perhaps better said conditional on) Off_relative_age
# figure_bh3.2.alternative <- orchaRd::orchard_plot(results_bh3.2.alternative,
#                                                   mod = "Off_relative_age",
#                                                   group = "StudyID",
#                                                   xlab = "Effect size",
#                                                   trunk.size = 6,
#                                                   branch.size = 2,
#                                                   twig.size = 1)


################################################################################
# 7 - METHODOLOGICAL HYPOTHESIS

# MH.1 - The relationship between maternal egg hormones and (maternal and offspring)
# fitness is different for studies that are correlational vs experimental.

################################################################################

# MH.1.1 - The effect should be stronger for experimental vs correlational 
# studies.

# Data base: the general data base (i.e., the one that includes maternal and 
# offspring fitness data)

################################################################################
table(meta.final_ok_ok$Study_type,
      meta.final_ok_ok$Hormone_measured_general)
# There are no correlational studies for TH3 and TH4. Given that we cannot compare
# the effect of thyroid hormones for correlational and experimental studies, 
# we will only study this question for androgens and glucocorticoids.


#Therefore, I need to first subset from the general data base, the rows that
# contain information only for androgens and glucocorticoids.
meta.final_ok_ok_androgens_glucocorticoids <- droplevels(subset(meta.final_ok_ok, 
                                                                meta.final_ok_ok$Hormone_measured_general != "TH3 and TH4"))


# # PHYLOGENETIC TREE:
# resolved_names_and_gcs <- tnrs_match_names(as.character(unique(meta.final_ok_ok_androgens_glucocorticoids$Species)))
# 
# # Saving the taxonomic data created on the 5th February 2025 to speed the
# # process in the future and allow full reproducibility
# save(resolved_names_and_gcs, file = "data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_and_gcs.RData")

# Loading the taxonomic data created on the 5th February 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_and_gcs.RData") 

# # Extracting phylogenetic information
# my_tree_and_gcs <- tol_induced_subtree(ott_ids =
#                                          resolved_names_and_gcs[,"ott_id"],
#                                        label_format = "name")
# # # Quick tree plotting
# # plot(my_tree_and_gcs, no.margin = TRUE)
# # 
# # We need to check for the existence of polytomies
# is.binary(my_tree_and_gcs)
# # Yes, meaning there are no polytomies. Let's go on.
# #
# # To confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# my_tree_and_gcs$tip.label <- gsub("_", " ", my_tree_and_gcs$tip.label)
# 
# intersect(as.character(my_tree_and_gcs$tip.label),
#           as.character(meta.final_ok_ok_androgens_glucocorticoids$Species))
# 
# # Listed in our database but not in the tree
# setdiff(as.character(meta.final_ok_ok_androgens_glucocorticoids$Species),
#         as.character(my_tree_and_gcs$tip.label))
# 
# # Listed in the tree but not in our database
# setdiff(as.character(my_tree_and_gcs$tip.label),
#         as.character(meta.final_ok_ok_androgens_glucocorticoids$Species))
# # No error or inconsistencies found.
# # 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# # 5th February 2025
# save(my_tree_and_gcs, file = "data/outputs/phylogenetic_files/tree_and_gcs.Rdata")

# We can now load the saved tree
load("data/outputs/phylogenetic_files/tree_and_gcs.Rdata") #my_tree_and_gcs

# # Compute branch lengths of tree
# phylo_branch_and_gcs <- compute.brlen(my_tree_and_gcs, method = "Grafen", power = 1)
# #
# # Check if tree is ultrametric
# is.ultrametric(phylo_branch_and_gcs)
# # TRUE
# 
# # Matrix to be included in the models
# phylo_cor_and_gcs <- vcv(phylo_branch_and_gcs, cor = T)
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(phylo_cor_and_gcs, file = "data/outputs/phylogenetic_files/phylo_cor_and_gcs.Rdata")

# we can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor_and_gcs.Rdata") #phylo_cor_and_gcs

# Creating a duplicate species variable for the phylogenetic analysis
meta.final_ok_ok_androgens_glucocorticoids$Species_phylo_and_gcs <- meta.final_ok_ok_androgens_glucocorticoids$Species

# VARIANCE-COVARIANCE MATRIX:
VCV_ESVar_and_gcs <- matrix(0, nrow = nrow(meta.final_ok_ok_androgens_glucocorticoids), 
                            ncol = nrow(meta.final_ok_ok_androgens_glucocorticoids))

# Names rows and columns for each obsID
rownames(VCV_ESVar_and_gcs) <- meta.final_ok_ok_androgens_glucocorticoids[, "EffectID"]
colnames(VCV_ESVar_and_gcs) <- meta.final_ok_ok_androgens_glucocorticoids[, "EffectID"]

# Finds effect sizes that come from the same study
shared_coord_and_gcs <- which(meta.final_ok_ok_androgens_glucocorticoids[, "StudyID"] %in% 
                                meta.final_ok_ok_androgens_glucocorticoids[duplicated(meta.final_ok_ok_androgens_glucocorticoids[, "StudyID"]), 
                                                                           "StudyID"] == TRUE)

combinations_and_gcs <- do.call("rbind", tapply(shared_coord_and_gcs, 
                                                meta.final_ok_ok_androgens_glucocorticoids[shared_coord_and_gcs, "StudyID"], 
                                                function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations_and_gcs)[1]) {
  p1_and_gcs <- combinations_and_gcs[i, 1]
  p2_and_gcs <- combinations_and_gcs[i, 2]
  p1_p2_cov_and_gcs <- 0.5 * sqrt(meta.final_ok_ok_androgens_glucocorticoids[p1_and_gcs, "cor_var"]) * 
    sqrt(meta.final_ok_ok_androgens_glucocorticoids[p2_and_gcs, "cor_var"])
  VCV_ESVar_and_gcs[p1_and_gcs, p2_and_gcs] <- p1_p2_cov_and_gcs
  VCV_ESVar_and_gcs[p2_and_gcs, p1_and_gcs] <- p1_p2_cov_and_gcs
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(VCV_ESVar_and_gcs) <- meta.final_ok_ok_androgens_glucocorticoids[, "cor_var"]

# # In case you want to visually double check the matrix outside of R
# write.csv(VCV_ESVar_and_gcs, 'data/outputs/variance-covariance_matrices/VCV_ESVar_and_gcs.csv')



# STATISTICAL MODEL:
meta.regression.mh1.1 <- rma.mv(cor,
                                VCV_ESVar_and_gcs,
                                mods = ~ Hormone_measured_general * Study_type,
                                random = list(~ 1 | StudyID,
                                              ~ 1 | LaboratoryID,
                                              ~ 1 | PopulationID,
                                              ~ 1 | Species,
                                              ~ 1 | Species_phylo_and_gcs,
                                              ~ 1 | EffectID),
                                method = "REML",
                                R = list(Species_phylo_and_gcs = phylo_cor_and_gcs),
                                test = "t",
                                data = meta.final_ok_ok_androgens_glucocorticoids)

#save(meta.regression.mh1.1, file = "data/outputs/statistical_models/meta_regression_mh1_1.RData")
load(file = "data/outputs/statistical_models/meta_regression_mh1_1.RData") #meta.regression.mh1.1

# Printing the summary results of the model
print(meta.regression.mh1.1, digits = 3)

# Calculate marginal R2 with r2_ml
R2.method_mh1.1 <- r2_ml(meta.regression.mh1.1)
round(R2.method_mh1.1 * 100, 1)


# In order to better understand and plot the results, I will create an artificial
# variable that includes this interaction.
meta.final_ok_ok_androgens_glucocorticoids$hormone_studytype <- paste(meta.final_ok_ok_androgens_glucocorticoids$Hormone_measured_general,
                                                                      meta.final_ok_ok_androgens_glucocorticoids$Study_type,
                                                                      sep = "_")

# I will now run a new model with this artificial variable.
meta.regression.mh1.1_artificial <- rma.mv(cor,
                                           VCV_ESVar_and_gcs,
                                           mods = ~ hormone_studytype - 1,
                                           random = list(~ 1 | StudyID,
                                                         ~ 1 | LaboratoryID,
                                                         ~ 1 | PopulationID,
                                                         ~ 1 | Species,
                                                         ~ 1 | Species_phylo_and_gcs,
                                                         ~ 1 | EffectID),
                                           method = "REML",
                                           R = list(Species_phylo_and_gcs = phylo_cor_and_gcs),
                                           test = "t",
                                           data = meta.final_ok_ok_androgens_glucocorticoids)

#save(meta.regression.mh1.1_artificial, file = "data/outputs/statistical_models/meta_regression_mh1_1_artificial.RData")
load(file = "data/outputs/statistical_models/meta_regression_mh1_1_artificial.RData") #meta.regression.mh1.1_artificial

# Printing the summary results of the model
print(meta.regression.mh1.1_artificial, digits = 3)

# Calculate marginal R2 with r2_ml
R2.method_mh1.1_artificial <- r2_ml(meta.regression.mh1.1_artificial)
round(R2.method_mh1.1_artificial * 100, 1)

# TABLE WITH RESULTS:
results_mh1.1_artificial <- orchaRd::mod_results(meta.regression.mh1.1_artificial, 
                                                 mod = "hormone_studytype", 
                                                 group = "StudyID", 
                                                 subset = TRUE)

round(as.data.frame(results_mh1.1_artificial[[1]])[,c(2:6)], 3)


# FIGURE: effect of each group of hormones on fitness traits.
fig_mh1.1_artificial <- orchaRd::orchard_plot(meta.regression.mh1.1_artificial, 
                                              mod = "hormone_studytype", 
                                              group = "StudyID", 
                                              xlab = "Effect size",
                                              trunk.size = 2,
                                              branch.size = 2,
                                              twig.size = 1)

# Experimental studies have a negative effect for both androgens and glucocticoids.
# Correlational studies show oppostive trends for each hormone  (and: neg; gcs: pos).
# None of these effects are statistically significant.


# Pair-wise comparisons including post-hoc Wald tests
meta.regression.mh1.1_artificial_pc <- rma.mv(cor,
                                              VCV_ESVar_and_gcs,
                                              mods = ~ hormone_studytype,
                                              random = list(~ 1 | StudyID,
                                                            ~ 1 | LaboratoryID,
                                                            ~ 1 | PopulationID,
                                                            ~ 1 | Species,
                                                            ~ 1 | Species_phylo_and_gcs,
                                                            ~ 1 | EffectID),
                                              method = "REML",
                                              R = list(Species_phylo_and_gcs = phylo_cor_and_gcs),
                                              test = "t",
                                              data = meta.final_ok_ok_androgens_glucocorticoids)

#save(meta.regression.mh1.1_artificial_pc, file = "data/outputs/statistical_models/meta_regression_mh1_1_artificial_pc.RData")
load(file = "data/outputs/statistical_models/meta_regression_mh1_1_artificial_pc.RData") #meta.regression.mh1.1_artificial_pc

print(meta.regression.mh1.1_artificial_pc, digits = 3)


# I now compare statistical differences between groups.
# For androgens, the difference between cor and exp studies can be obtained from
# the output of the model. For glucocorticoids, I do a Post-hoc Wald test:
# Glucocorticoids - correlational vs experimental: 
car::linearHypothesis(meta.regression.mh1.1_artificial_pc, rbind(c(0,0,1,-1)))

################################################################################

# MH.1.2: The type of experimental manipulation can have differential fitness
# consequences (stronger effect for direct manipulations to the egg compared to
# the female)

# Data base: I work with the data base that is only for the offspring.

################################################################################

# DATA BASE: I need to subset the data base to work only with experimental studies.
meta.final_ok_ok_off_exp <- droplevels(subset(meta.final_ok_ok_off, 
                                              meta.final_ok_ok_off$Study_type == "experimental"))

nrow(meta.final_ok_ok_off_exp)

# IMPORTANT: The column name 'Exp_mother_egg" has two levels: mother and egg; they
# provide information on whether the experiment was done on the mother or on the
# egg. Now, in this column there are also NAs. These 'NAs' are because in these 
# studies the experimental manipulation was not done directly on the mother, but
# rather indirectly (e.g., by manipulating the number of males around). However, 
# since the aim was to manipulate the state of the mother to see how this affects
# egg hormone deposition, I will replace these 'NAs' by mother.
nrow_na_exp <- which(is.na(meta.final_ok_ok_off_exp$Exp_mother_egg))
meta.final_ok_ok_off_exp$Exp_mother_egg[nrow_na_exp] <- "mother"


table(meta.final_ok_ok_off_exp$Exp_mother_egg,
      meta.final_ok_ok_off_exp$Hormone_measured_general)
# There are no levels for manipulations done on the mother for THs. Hence, we 
# cannot include thyroid hormones in our analysis. 

# I will subset the database:
meta.final_ok_ok_off_exp_and_gcs <- droplevels(subset(meta.final_ok_ok_off_exp, 
                                                      meta.final_ok_ok_off_exp$Hormone_measured_general != "TH3 and TH4"))

# # PHYLOGENETIC TREE:
# resolved_names_exp_and_gcs <- tnrs_match_names(as.character(unique(meta.final_ok_ok_off_exp_and_gcs$Species)))
# 
# # Saving the taxonomic data created on the 5th February 2025 to speed the
# # process in the future and allow full reproducibility
# save(resolved_names_exp_and_gcs, file = "data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_exp_and_gcs.RData")

# Loading the taxonomic data created on the 5th February 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_exp_and_gcs.RData") #resolved_names_exp

# # extracting phylogenetic information
# my_tree_exp_and_gcs <- tol_induced_subtree(ott_ids =
#                                           resolved_names_exp_and_gcs[,"ott_id"],
#                                           label_format = "name")
# # # Quick tree plotting
# # plot(my_tree_exp_and_gcs, no.margin = TRUE)
# 
# # We need to check for the existence of polytomies
# is.binary(my_tree_exp_and_gcs)
# # Yes, meaning there are no polytomies. Let's go on.
# 
# # To confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# my_tree_exp_and_gcs$tip.label <- gsub("_", " ", my_tree_exp_and_gcs$tip.label)
# 
# intersect(as.character(my_tree_exp_and_gcs$tip.label),
#          as.character(meta.final_ok_ok_off_exp_and_gcs$Species))
# 
# # Listed in our database but not in the tree
# setdiff(as.character(meta.final_ok_ok_off_exp_and_gcs$Species),
#       as.character(my_tree_exp_and_gcs$tip.label))
# 
# # Listed in the tree but not in our database
# setdiff(as.character(my_tree_exp_and_gcs$tip.label),
#       as.character(meta.final_ok_ok_off_exp_and_gcs$Species))
# # No error or inconsistencies found.
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(my_tree_exp_and_gcs, file = "data/outputs/phylogenetic_files/tree_exp_and_gcs.Rdata")

# We can now load the saved tree
load("data/outputs/phylogenetic_files/tree_exp_and_gcs.Rdata") #my_tree_exp

# # Compute branch lengths of tree
# phylo_branch_exp_and_gcs <- compute.brlen(my_tree_exp_and_gcs, method = "Grafen", power = 1)
# 
# # Check if tree is ultrametric
# is.ultrametric(phylo_branch_exp_and_gcs)
# # TRUE
# 
# # Matrix to be included in the models
# phylo_cor_exp_and_gcs <- vcv(phylo_branch_exp_and_gcs, cor = T)
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(phylo_cor_exp_and_gcs, file = "data/outputs/phylogenetic_files/phylo_cor_exp_and_gcs.Rdata")

# we can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor_exp_and_gcs.Rdata") #phylo_cor_exp

# Creating a duplicate species variable for the phylogenetic analysis
meta.final_ok_ok_off_exp_and_gcs$Species_phylo_exp_and_gcs <- meta.final_ok_ok_off_exp_and_gcs$Species


# VARIANCE-COVARIANCE MATRIX:
VCV_ESVar_exp_and_gcs <- matrix(0, nrow = nrow(meta.final_ok_ok_off_exp_and_gcs), 
                                ncol = nrow(meta.final_ok_ok_off_exp_and_gcs))

# Names rows and columns for each obsID
rownames(VCV_ESVar_exp_and_gcs) <- meta.final_ok_ok_off_exp_and_gcs[, "EffectID"]
colnames(VCV_ESVar_exp_and_gcs) <- meta.final_ok_ok_off_exp_and_gcs[, "EffectID"]

# Finds effect sizes that come from the same study
shared_coord_exp_and_gcs <- which(meta.final_ok_ok_off_exp_and_gcs[, "StudyID"] %in% 
                                    meta.final_ok_ok_off_exp_and_gcs[duplicated(meta.final_ok_ok_off_exp_and_gcs[, "StudyID"]), 
                                                                     "StudyID"] == TRUE)

combinations_exp_and_gcs <- do.call("rbind", tapply(shared_coord_exp_and_gcs, 
                                                    meta.final_ok_ok_off_exp_and_gcs[shared_coord_exp_and_gcs, "StudyID"], 
                                                    function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations_exp_and_gcs)[1]) {
  p1_exp_and_gcs <- combinations_exp_and_gcs[i, 1]
  p2_exp_and_gcs <- combinations_exp_and_gcs[i, 2]
  p1_p2_cov_exp_and_gcs <- 0.5 * sqrt(meta.final_ok_ok_off_exp_and_gcs[p1_exp_and_gcs, "cor_var"]) * 
    sqrt(meta.final_ok_ok_off_exp_and_gcs[p2_exp_and_gcs, "cor_var"])
  VCV_ESVar_exp_and_gcs[p1_exp_and_gcs, p2_exp_and_gcs] <- p1_p2_cov_exp_and_gcs
  VCV_ESVar_exp_and_gcs[p2_exp_and_gcs, p1_exp_and_gcs] <- p1_p2_cov_exp_and_gcs
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(VCV_ESVar_exp_and_gcs) <- meta.final_ok_ok_off_exp_and_gcs[, "cor_var"]

# # In case you want to visually double check the matrix outside of R
# write.csv(VCV_ESVar_exp_and_gcs, 'data/outputs/variance-covariance_matrices/VCV_ESVar_exp_and_gcs.csv')


# STATISTICAL ANALYSIS:
meta.regression.mh1.2 <- rma.mv(cor,
                                VCV_ESVar_exp_and_gcs,
                                mods = ~ Hormone_measured_general * Exp_mother_egg,
                                random = list(~ 1 | StudyID,
                                              ~ 1 | LaboratoryID,
                                              ~ 1 | PopulationID,
                                              ~ 1 | Species,
                                              ~ 1 | Species_phylo_exp_and_gcs,
                                              ~ 1 | EffectID),
                                method = "REML",
                                R = list(Species_phylo_exp_and_gcs = phylo_cor_exp_and_gcs),
                                test = "t",
                                data = meta.final_ok_ok_off_exp_and_gcs)

#save(meta.regression.mh1.2, file = "data/outputs/statistical_models/meta_regression_mh1_2.RData")
load(file = "data/outputs/statistical_models/meta_regression_mh1_2.RData") #meta.regression.mh1.2

# Printing the summary results of the model
print(meta.regression.mh1.2, digits = 3)

# Calculate marginal R2 with r2_ml
R2.method_mh1.2 <- r2_ml(meta.regression.mh1.2)
round(R2.method_mh1.2 * 100, 1)


# Like I did for the previous hypothesis, I will create an artificial variable
# that includes the interaction between the hormone meassured and the time when
# the manipulation was done (i.e., mother or egg)
meta.final_ok_ok_off_exp_and_gcs$hormone_exp_stage <- paste(meta.final_ok_ok_off_exp_and_gcs$Hormone_measured_general,
                                                            meta.final_ok_ok_off_exp_and_gcs$Exp_mother_egg,
                                                            sep = "_")


# STATISTICAL MODEL:
meta.regression.mh1.2_artificial <- rma.mv(cor,
                                           VCV_ESVar_exp_and_gcs,
                                           mods = ~ hormone_exp_stage - 1,
                                           random = list(~ 1 | StudyID,
                                                         ~ 1 | LaboratoryID,
                                                         ~ 1 | PopulationID,
                                                         ~ 1 | Species,
                                                         ~ 1 | Species_phylo_exp_and_gcs,
                                                         ~ 1 | EffectID),
                                           method = "REML",
                                           R = list(Species_phylo_exp_and_gcs = phylo_cor_exp_and_gcs),
                                           test = "t",
                                           data = meta.final_ok_ok_off_exp_and_gcs)

#save(meta.regression.mh1.2_artificial, file = "data/outputs/statistical_models/meta_regression_mh1_2_artificial.RData")
load(file = "data/outputs/statistical_models/meta_regression_mh1_2_artificial.RData") #meta.regression.mh1.2_artificial

# Printing the summary results of the model
print(meta.regression.mh1.2_artificial, digits = 3)

# Calculate marginal R2 with r2_ml
R2.method_mh1.2_artificial <- r2_ml(meta.regression.mh1.2_artificial)
round(R2.method_mh1.2_artificial * 100, 1)


# TABLE WITH RESULTS:
results_mh1.2_artificial <- orchaRd::mod_results(meta.regression.mh1.2_artificial, 
                                                 mod = "hormone_exp_stage", 
                                                 group = "StudyID", 
                                                 subset = TRUE)

round(as.data.frame(results_mh1.2_artificial[[1]])[,c(2:6)], 3)


# FIGURE: effect of each group of hormones on fitness traits.
fig_mh1.2_artificial <- orchaRd::orchard_plot(meta.regression.mh1.2_artificial, 
                                              mod = "hormone_exp_stage", 
                                              group = "StudyID", 
                                              xlab = "Effect size",
                                              trunk.size = 2,
                                              branch.size = 2,
                                              twig.size = 1)

# Both hormones have a negative effect independently of the type of experimental 
# manipulation. 
# None of these effects are statistically significant.


# Pair-wise comparisons including post-hoc Wald tests
meta.regression.mh1.2_artificial_pc <- rma.mv(cor,
                                              VCV_ESVar_exp_and_gcs,
                                              mods = ~ hormone_exp_stage,
                                              random = list(~ 1 | StudyID,
                                                            ~ 1 | LaboratoryID,
                                                            ~ 1 | PopulationID,
                                                            ~ 1 | Species,
                                                            ~ 1 | Species_phylo_exp_and_gcs,
                                                            ~ 1 | EffectID),
                                              method = "REML",
                                              R = list(Species_phylo_exp_and_gcs = phylo_cor_exp_and_gcs),
                                              test = "t",
                                              data = meta.final_ok_ok_off_exp_and_gcs)

#save(meta.regression.mh1.2_artificial_pc, file = "data/outputs/statistical_models/meta_regression_mh1_2_artificial_pc.RData")
load(file = "data/outputs/statistical_models/meta_regression_mh1_2_artificial_pc.RData") #meta.regression.mh1.2_artificial_pc


print(meta.regression.mh1.2_artificial_pc, digits = 3)


# I now test if there are statistical differences between groups. For androgens, 
# I can obtain this information from he output of the model. For glucocorticoids, 
# I need to run a Post-hoc Wald test:
car::linearHypothesis(meta.regression.mh1.2_artificial_pc, rbind(c(0, 0, 1, -1)))



################################################################################
# MH.1.3 - The dose can also influence the results. Studies done with hormone
# modifications within the natural range of the hormone should have a stronger
# effect than studies using supra-physiological values.

# Data base: I work with the data base that is only for the offspring.

################################################################################

# This information is divided for experimental studies done in the eggs vs the mom.
# I therefore need to create a new variable containing this information.

meta.final_ok_ok_off_exp$Exp_dose <- "NA"

meta.final_ok_ok_off_exp$Exp_dose <- as.character(meta.final_ok_ok_off_exp$Exp_dose)
meta.final_ok_ok_off_exp$Horm_mother_cateogrical_dose <- as.character(meta.final_ok_ok_off_exp$Horm_mother_cateogrical_dose)
meta.final_ok_ok_off_exp$Horm_egg_cateogrical_dose <- as.character(meta.final_ok_ok_off_exp$Horm_egg_cateogrical_dose)

meta.final_ok_ok_off_exp$Exp_dose <- meta.final_ok_ok_off_exp$Horm_mother_cateogrical_dose

nrow_na_egg_dose <- which(is.na(meta.final_ok_ok_off_exp$Exp_dose))
meta.final_ok_ok_off_exp$Exp_dose[nrow_na_egg_dose] <- meta.final_ok_ok_off_exp$Horm_egg_cateogrical_dose[nrow_na_egg_dose]

meta.final_ok_ok_off_exp$Exp_dose <- as.character(meta.final_ok_ok_off_exp$Exp_dose)

############
# There are some cells for which we do not have information (i.e., 'NAs'). I will
# remove them from the final data base
meta.final_ok_ok_off_exp_dose <- droplevels(subset(meta.final_ok_ok_off_exp, 
                                                   meta.final_ok_ok_off_exp$Exp_dose != "NA"))

# In the final column, we have 3 levels: "supraphysiological", "within natural
# range" and "within natural range (but probably not for all)". I will replace
# the last category for "within natural range" because the authors claimed in their
# papers that the manipulations were done within that range.

meta.final_ok_ok_off_exp_dose$Exp_dose[meta.final_ok_ok_off_exp_dose$Exp_dose == 
                                         "within natural range (but probably not for all)"] <- "within natural range"

meta.final_ok_ok_off_exp_dose$Exp_dose <- as.factor(meta.final_ok_ok_off_exp_dose$Exp_dose)
nrow(meta.final_ok_ok_off_exp_dose)
length(unique(meta.final_ok_ok_off_exp_dose$Species))


# I now check if I have enough levels for testing the interactions in the model.
table(meta.final_ok_ok_off_exp_dose$Hormone_measured_general,
      meta.final_ok_ok_off_exp_dose$Exp_dose)

# I do not have data points for testing how supraphysiological concentrations of 
# corticosterone (i.e., the hormone that is expected to have a negative effect 
# on fitness) and THs influences fitness. The only group of hormones for which we
# have enough data to be compared are androgens. Therefore, I will only work with
# the columns that have information for androgens and the type of manipulation.

meta.final_ok_ok_off_exp_dose_androgens <- droplevels(subset(meta.final_ok_ok_off_exp_dose, 
                                                             meta.final_ok_ok_off_exp_dose$Hormone_measured_general == "androgens"))
nrow(meta.final_ok_ok_off_exp_dose_androgens)
length(unique(meta.final_ok_ok_off_exp_dose_androgens$Species))


# # PHYLOGENETIC TREE:
# # I need to create a new phylogenetic tree and covariance matrix because I have
# # a smaller data set.
# resolved_names_exp_dose_androgens <- tnrs_match_names(as.character(unique(meta.final_ok_ok_off_exp_dose_androgens$Species)))
# 
# # Saving the taxonomic data created on the 5th February 2025 to speed the
# # process in the future and allow full reproducibility
# save(resolved_names_exp_dose_androgens, file = "data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_exp_dose_androgens.RData")

# Loading the taxonomic data created on the 5th February 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_exp_dose_androgens.RData") #resolved_names_exp_dose_androgens

# # extracting phylogenetic information
# my_tree_exp_dose_androgens <- tol_induced_subtree(ott_ids =
#                                                     resolved_names_exp_dose_androgens[,"ott_id"],
#                                                   label_format = "name")
# 
# # # Quick tree plotting
# # plot(my_tree_exp_dose_androgens, no.margin = TRUE)
# 
# # We need to check for the existence of polytomies
# is.binary(my_tree_exp_dose_androgens)
# # Yes, meaning there are no polytomies. Let's go on.
# 
# # To confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# my_tree_exp_dose_androgens$tip.label <- gsub("_", " ",
#                                              my_tree_exp_dose_androgens$tip.label)
# 
# intersect(as.character(my_tree_exp_dose_androgens$tip.label),
#           as.character(meta.final_ok_ok_off_exp_dose_androgens$Species))
# 
# # Listed in our database but not in the tree
# setdiff(as.character(meta.final_ok_ok_off_exp_dose_androgens$Species),
#         as.character(my_tree_exp_dose_androgens$tip.label))
# 
# # Listed in the tree but not in our database
# setdiff(as.character(my_tree_exp_dose_androgens$tip.label),
#         as.character(meta.final_ok_ok_off_exp_dose_androgens$Species))
# # No error or inconsistencies found.
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(my_tree_exp_dose_androgens, file = "data/outputs/phylogenetic_files/tree_exp_dose_androgens.Rdata")

# We can now load the saved tree
load("data/outputs/phylogenetic_files/tree_exp_dose_androgens.Rdata") #my_tree_exp_dose_androgens

# # Compute branch lengths of tree
# phylo_branch_exp_dose_androgens <- compute.brlen(my_tree_exp_dose_androgens,
#                                                  method = "Grafen", power = 1)
# 
# # Check if tree is ultrametric
# is.ultrametric(phylo_branch_exp_dose_androgens)
# # TRUE
# 
# # Matrix to be included in the models
# phylo_cor_exp_dose_androgens <- vcv(phylo_branch_exp_dose_androgens, cor = T)
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(phylo_cor_exp_dose_androgens, file = "data/outputs/phylogenetic_files/phylo_cor_exp_dose_androgens.Rdata")

# we can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor_exp_dose_androgens.Rdata") #phylo_cor_exp_dose_androgens

# Creating a duplicate species variable for the phylogenetic analysis
meta.final_ok_ok_off_exp_dose_androgens$Species_phylo_exp_dose_androgens <- meta.final_ok_ok_off_exp_dose_androgens$Species


# VARIANCE-COVARIANCE MATRIX:
VCV_ESVar_exp_dose_androgens <- matrix(0, nrow = nrow(meta.final_ok_ok_off_exp_dose_androgens), 
                                       ncol = nrow(meta.final_ok_ok_off_exp_dose_androgens))

# Names rows and columns for each obsID
rownames(VCV_ESVar_exp_dose_androgens) <- meta.final_ok_ok_off_exp_dose_androgens[, "EffectID"]
colnames(VCV_ESVar_exp_dose_androgens) <- meta.final_ok_ok_off_exp_dose_androgens[, "EffectID"]

# Finds effect sizes that come from the same study
shared_coord_exp_dose_androgens <- which(meta.final_ok_ok_off_exp_dose_androgens[, "StudyID"] %in% 
                                           meta.final_ok_ok_off_exp_dose_androgens[duplicated(meta.final_ok_ok_off_exp_dose_androgens[, "StudyID"]), 
                                                                                   "StudyID"] == TRUE)

combinations_exp_dose_androgens <- do.call("rbind", tapply(shared_coord_exp_dose_androgens, 
                                                           meta.final_ok_ok_off_exp_dose_androgens[shared_coord_exp_dose_androgens, "StudyID"], 
                                                           function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations_exp_dose_androgens)[1]) {
  p1_exp_dose_androgens <- combinations_exp_dose_androgens[i, 1]
  p2_exp_dose_androgens <- combinations_exp_dose_androgens[i, 2]
  p1_p2_cov_exp_dose_androgens <- 0.5 * sqrt(meta.final_ok_ok_off_exp_dose_androgens[p1_exp_dose_androgens, "cor_var"]) * 
    sqrt(meta.final_ok_ok_off_exp_dose_androgens[p2_exp_dose_androgens, "cor_var"])
  VCV_ESVar_exp_dose_androgens[p1_exp_dose_androgens, p2_exp_dose_androgens] <- p1_p2_cov_exp_dose_androgens
  VCV_ESVar_exp_dose_androgens[p2_exp_dose_androgens, p1_exp_dose_androgens] <- p1_p2_cov_exp_dose_androgens
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(VCV_ESVar_exp_dose_androgens) <- meta.final_ok_ok_off_exp_dose_androgens[, "cor_var"]

# # In case you want to visually double check the matrix outside of R
# write.csv(VCV_ESVar_exp_dose_androgens, 'data/outputs/variance-covariance_matrices/VCV_ESVar_exp_dose_androgens.csv')


# STATISTICAL ANALYSIS:
meta.regression.mh1.3 <- rma.mv(cor,
                                VCV_ESVar_exp_dose_androgens,
                                mods = ~ Exp_dose - 1,
                                random = list(~ 1 | StudyID,
                                              ~ 1 | LaboratoryID,
                                              ~ 1 | PopulationID,
                                              ~ 1 | Species,
                                              ~ 1 | Species_phylo_exp_dose_androgens,
                                              ~ 1 | EffectID),
                                method = "REML",
                                R = list(Species_phylo_exp_dose_androgens = phylo_cor_exp_dose_androgens),
                                test = "t",
                                data = meta.final_ok_ok_off_exp_dose_androgens)

#save(meta.regression.mh1.3, file = "data/outputs/statistical_models/meta_regression_mh1_3.RData")
load(file = "data/outputs/statistical_models/meta_regression_mh1_3.RData") #meta.regression.mh1.3


# Printing the summary results of the model
print(meta.regression.mh1.3, digits = 3)

# Calculate marginal R2 with r2_ml
R2.method_mh1.3 <- r2_ml(meta.regression.mh1.3)
round(R2.method_mh1.3 * 100, 2)


# TABLE WITH RESULTS:
results_mh1.3 <- orchaRd::mod_results(meta.regression.mh1.3, 
                                      mod = "Exp_dose", 
                                      group = "StudyID", 
                                      subset = TRUE)

round(as.data.frame(results_mh1.3[[1]])[,c(2:6)], 3)


# FIGURE: effect of each group of hormones on fitness traits.
fig_mh1.3 <- orchaRd::orchard_plot(meta.regression.mh1.3, 
                                   mod = "Exp_dose", 
                                   group = "StudyID", 
                                   xlab = "Effect size",
                                   trunk.size = 2,
                                   branch.size = 2,
                                   twig.size = 1)

# Androgen hormones have positive effects on offpring fitness only under 
# supraphysiological conditions. If not, the effect is negative.
# No statistical significance.


# STATISTICAL ANALYSIS to check for statistical differences between the two groups:
meta.regression.mh1.3_pc <- rma.mv(cor,
                                   VCV_ESVar_exp_dose_androgens,
                                   mods = ~ Exp_dose,
                                   random = list(~ 1 | StudyID,
                                                 ~ 1 | LaboratoryID,
                                                 ~ 1 | PopulationID,
                                                 ~ 1 | Species,
                                                 ~ 1 | Species_phylo_exp_dose_androgens,
                                                 ~ 1 | EffectID),
                                   method = "REML",
                                   R = list(Species_phylo_exp_dose_androgens = phylo_cor_exp_dose_androgens),
                                   test = "t",
                                   data = meta.final_ok_ok_off_exp_dose_androgens)

#save(meta.regression.mh1.3_pc, file = "data/outputs/statistical_models/meta_regression_mh1_3_pc.RData")
load(file = "data/outputs/statistical_models/meta_regression_mh1_3_pc.RData") #meta.regression.mh1.3_pc


# Printing the summary results of the model
print(meta.regression.mh1.3_pc, digits = 3)



################################################################################
# 8) SENSITIVITY ANALYSES

# We will perform two sensitivity analyses. One that was not pre-registered but
# we consider important to do base on the transformation of the effect sizes we
# decided to do, and the other one following the pre-registration.

# 8.1 -  In the first one, we will check how strong is the effect size of choice 
# (Pearson or Biserial correlations) affecting the results. For this, we will
# run three intercept-only meta-analytic models: one a) using only the data from 
# studies were we transform the data into a Pearson correlation (193 rows), 
# a second one b) using the papers where we did Biserial correlations (277 rows), and
# a third one c) using log-response ratio coefficients as a response variable instead
# of Biserial correlation. 

# 8.2 - In the second one, that follow the pre-registration, we will test the robustness
# of the results for predictions BH.1.1-BH.1.5 by fitting a phylogenetic 
# multilevel meta-regression with the same random effects structure as used above
# and hormone type as a moderator. The data base we use to run this analysis is 
# the one that only includes offspring data.
################################################################################

# 8.1. Sensitivity analysis to check for the influence of the chosen effect sizes.

# I first need to separate the data base: the studies from where we transformed 
# the data into Pearson correlations from the ones where we used Biserial cor.

# a) Data base that contains Pearson correlations: 
meta.final_ok_sensitivity_pearson <- droplevels(subset(meta.final_ok_ok, 
                                                       meta.final_ok_ok$Effect_size_type != "Mean"))

nrow(meta.final_ok_sensitivity_pearson)


# I will now run an intercept-only meta-analytic model. But before doing this
# I need to check the species that are present in this subset. There are 2 species
# less than in the general data base (18 vs 20), and also less entries. I will 
# therefore calculate a new phylogenetic tree and variance-covariance matrix.

# # PHYLOGENETIC TREE:
# resolved_names_pearson <- tnrs_match_names(as.character(unique(meta.final_ok_sensitivity_pearson$Species)))
# 
# # Saving the taxonomic data created on the 5th February 2025 to speed the
# # process in the future and allow full reproducibility
# save(resolved_names_pearson, file = "data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_pearson.RData")

# Loading the taxonomic data created on the 5th February 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_pearson.RData") #resolved_names_pearson

# # extracting phylogenetic information
# my_tree_pearson <- tol_induced_subtree(ott_ids =
#                                          resolved_names_pearson[,"ott_id"],
#                                        label_format = "name")
# 
# # # Quick tree plotting
# # plot(my_tree_pearson, no.margin = TRUE)
# 
# # We need to check for the existence of polytomies
# is.binary(my_tree_pearson)
# # Yes, meaning there are no polytomies. Let's go on.
# 
# # To confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# my_tree_pearson$tip.label <- gsub("_", " ", my_tree_pearson$tip.label)
# 
# intersect(as.character(my_tree_pearson$tip.label),
#           as.character(meta.final_ok_sensitivity_pearson$Species))
# 
# # Listed in our database but not in the tree
# setdiff(as.character(meta.final_ok_sensitivity_pearson$Species),
#         as.character(my_tree_pearson$tip.label))
# 
# # Listed in the tree but not in our database
# setdiff(as.character(my_tree_pearson$tip.label),
#         as.character(meta.final_ok_sensitivity_pearson$Species))
# # No error or inconsistencies found.
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# # 5th February 2025
# save(my_tree_pearson, file = "data/outputs/phylogenetic_files/tree_sens.pearson.Rdata")

# We can now load the saved tree
load("data/outputs/phylogenetic_files/tree_sens.pearson.Rdata") #my_tree_pearson

# # Compute branch lengths of tree
# phylo_branch_pearson <- compute.brlen(my_tree_pearson,
#                                       method = "Grafen", power = 1)
# 
# # Check if tree is ultrametric
# is.ultrametric(phylo_branch_pearson)
# # TRUE
# 
# # Matrix to be included in the models
# phylo_cor_pearson <- vcv(phylo_branch_pearson, cor = T)
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(phylo_cor_pearson, file = "data/outputs/phylogenetic_files/phylo_cor_pearson.Rdata")

# we can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor_pearson.Rdata") #phylo_cor_pearson

# Creating a duplicate species variable for the phylogenetic analysis
meta.final_ok_sensitivity_pearson$Species_phylo_pearson <- 
  meta.final_ok_sensitivity_pearson$Species



# VARIANCE-COVARIANCE MATRIX:
VCV_ESVar_pearson <- matrix(0, nrow = nrow(meta.final_ok_sensitivity_pearson), 
                            ncol = nrow(meta.final_ok_sensitivity_pearson))

# Names rows and columns for each obsID
rownames(VCV_ESVar_pearson) <- meta.final_ok_sensitivity_pearson[, "EffectID"]
colnames(VCV_ESVar_pearson) <- meta.final_ok_sensitivity_pearson[, "EffectID"]

# Finds effect sizes that come from the same study
shared_coord_pearson <- which(meta.final_ok_sensitivity_pearson[, "StudyID"] %in% 
                                meta.final_ok_sensitivity_pearson[duplicated(meta.final_ok_sensitivity_pearson[, "StudyID"]), 
                                                                  "StudyID"] == TRUE)

combinations_pearson <- do.call("rbind", tapply(shared_coord_pearson, 
                                                meta.final_ok_sensitivity_pearson[shared_coord_pearson, "StudyID"], 
                                                function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations_pearson)[1]) {
  p1_pearson <- combinations_pearson[i, 1]
  p2_pearson <- combinations_pearson[i, 2]
  p1_p2_cov_pearson <- 0.5 * sqrt(meta.final_ok_sensitivity_pearson[p1_pearson, "cor_var"]) * 
    sqrt(meta.final_ok_sensitivity_pearson[p2_pearson, "cor_var"])
  VCV_ESVar_pearson[p1_pearson, p2_pearson] <- p1_p2_cov_pearson
  VCV_ESVar_pearson[p2_pearson, p1_pearson] <- p1_p2_cov_pearson
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(VCV_ESVar_pearson) <- meta.final_ok_sensitivity_pearson[, "cor_var"]

# In case you want to visually double check the matrix outside of R
#write.csv(VCV_ESVar_pearson, 'data/outputs/variance-covariance_matrices/VCV_ESVar_pearson.csv')


# STATISTICAL ANALYSIS: 
sensitivy.model.intercept.pearson <- rma.mv(cor,
                                            VCV_ESVar_pearson,
                                            mods = ~ 1,
                                            random = list(~ 1 | StudyID,
                                                          ~ 1 | LaboratoryID,
                                                          ~ 1 | PopulationID,
                                                          ~ 1 | Species,
                                                          ~ 1 | Species_phylo_pearson,
                                                          ~ 1 | EffectID),
                                            method = "REML",
                                            R = list(Species_phylo_pearson = phylo_cor_pearson),
                                            test = "t",
                                            data = meta.final_ok_sensitivity_pearson)
#save(sensitivy.model.intercept.pearson, file = "data/outputs/statistical_models/sensitivy_model_intercept_pearson.RData")
load(file = "data/outputs/statistical_models/sensitivy_model_intercept_pearson.RData") #sensitivy.model.intercept.pearson


# Printing the summary results of the model
print(sensitivy.model.intercept.pearson, digits = 3)

# Model funnel plots
par(mfrow = c(1, 1))
funnel(sensitivy.model.intercept.pearson)

# Estimating heterogeneity as I2 (Nakagawa and Santos 2012)
I2.model_pearson <- i2_ml(sensitivy.model.intercept.pearson)
round(I2.model_pearson, 1)


# TABLE:
results_sensitivy.model.intercept.pearson <- orchaRd::mod_results(sensitivy.model.intercept.pearson, 
                                                                  group = "StudyID", 
                                                                  subset = TRUE)

round(as.data.frame(results_sensitivy.model.intercept.pearson[[1]])[,c(2:6)], 3)


# FIGURE: effect of each group of hormones on fitness traits.
fig_hormones_fitness_intercept.pearson <- orchaRd::orchard_plot(sensitivy.model.intercept.pearson, 
                                                                group = "StudyID", 
                                                                xlab = "Effect size",
                                                                trunk.size = 2,
                                                                branch.size = 2,
                                                                twig.size = 1)

# When considering only those studies that have effect sizes as Pearson correlations
# we obtain a different, yet comparable,result than when pooling both Pearson and 
# pearson correlations. By considering only studies doing Pearson correlation we 
# obtain a small, positive effect size (0.04), whereas by pooling both models
# we obtain a small negative effect size (-0.07). In both models effects are
# statistically not significant.


##########################
# b) Data base that contains the data base with Biserial correlations

meta.final_ok_sensitivity_biserial <- droplevels(subset(meta.final_ok_ok, 
                                                        meta.final_ok_ok$Effect_size_type == "Mean"))

nrow(meta.final_ok_sensitivity_biserial)

length(unique(meta.final_ok_sensitivity_biserial$StudyID))


# I will now run an intercept-only meta-analytic model. But before doing this
# I need to check the species that are present in this subset. There are 4 species
# less than in the general data base (16 vs 20), and also less entries. I will 
# therefore calculate a new phylogenetic tree and variance-covariance matrix.

# # PHYLOGENETIC TREE:
# resolved_names_biserial <- tnrs_match_names(as.character(unique(meta.final_ok_sensitivity_biserial$Species)))
# 
# # Saving the taxonomic data created on the 5th February 2025 to speed the
# # process in the future and allow full reproducibility
# save(resolved_names_biserial, file = "data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_biserial.RData")

# Loading the taxonomic data created on the 5th February 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_biserial.RData") #resolved_names_biserial

# # extracting phylogenetic information
# my_tree_biserial <- tol_induced_subtree(ott_ids =
#                                           resolved_names_biserial[,"ott_id"],
#                                         label_format = "name")
# 
# # # Quick tree plotting
# # plot(my_tree_biserial, no.margin = TRUE)
# 
# # We need to check for the existence of polytomies
# is.binary(my_tree_biserial)
# # Yes, meaning there are no polytomies. Let's go on.
# 
# # To confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# my_tree_biserial$tip.label <- gsub("_", " ", my_tree_biserial$tip.label)
# 
# intersect(as.character(my_tree_biserial$tip.label),
#           as.character(meta.final_ok_sensitivity_biserial$Species))
# 
# # Listed in our database but not in the tree
# setdiff(as.character(meta.final_ok_sensitivity_biserial$Species),
#         as.character(my_tree_biserial$tip.label))
# #
# # Listed in the tree but not in our database
# setdiff(as.character(my_tree_biserial$tip.label),
#         as.character(meta.final_ok_sensitivity_biserial$Species))
# # No error or inconsistencies found.
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# # 5th February 2025
# save(my_tree_biserial, file = "data/outputs/phylogenetic_files/tree_sens.pearson.Rdata")

# We can now load the saved tree
load("data/outputs/phylogenetic_files/tree_sens.pearson.Rdata") #my_tree_biserial

# # Compute branch lengths of tree
# phylo_branch_biserial <- compute.brlen(my_tree_biserial,
#                                        method = "Grafen", power = 1)
# 
# # Check if tree is ultrametric
# is.ultrametric(phylo_branch_biserial)
# # TRUE
# 
# # Matrix to be included in the models
# phylo_cor_biserial <- vcv(phylo_branch_biserial, cor = T)
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(phylo_cor_biserial, file = "data/outputs/phylogenetic_files/phylo_cor_biserial.Rdata")

# we can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor_biserial.Rdata") #phylo_cor

# Creating a duplicate species variable for the phylogenetic analysis
meta.final_ok_sensitivity_biserial$Species_phylo_biserial <- 
  meta.final_ok_sensitivity_biserial$Species


# VARIANCE-COVARIANCE MATRIX:
VCV_ESVar_biserial <- matrix(0, nrow = nrow(meta.final_ok_sensitivity_biserial), 
                             ncol = nrow(meta.final_ok_sensitivity_biserial))

# Names rows and columns for each obsID
rownames(VCV_ESVar_biserial) <- meta.final_ok_sensitivity_biserial[, "EffectID"]
colnames(VCV_ESVar_biserial) <- meta.final_ok_sensitivity_biserial[, "EffectID"]

# Finds effect sizes that come fbiserial the same study
shared_coord_biserial <- which(meta.final_ok_sensitivity_biserial[, "StudyID"] %in% 
                                 meta.final_ok_sensitivity_biserial[duplicated(meta.final_ok_sensitivity_biserial[, "StudyID"]), 
                                                                    "StudyID"] == TRUE)

combinations_biserial <- do.call("rbind", tapply(shared_coord_biserial, 
                                                 meta.final_ok_sensitivity_biserial[shared_coord_biserial, "StudyID"], 
                                                 function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations_biserial)[1]) {
  p1_biserial <- combinations_biserial[i, 1]
  p2_biserial <- combinations_biserial[i, 2]
  p1_p2_cov_biserial <- 0.5 * sqrt(meta.final_ok_sensitivity_biserial[p1_biserial, "cor_var"]) * 
    sqrt(meta.final_ok_sensitivity_biserial[p2_biserial, "cor_var"])
  VCV_ESVar_biserial[p1_biserial, p2_biserial] <- p1_p2_cov_biserial
  VCV_ESVar_biserial[p2_biserial, p1_biserial] <- p1_p2_cov_biserial
} 


# Enters previously calculated effect size sampling variances into diagonals 
diag(VCV_ESVar_biserial) <- meta.final_ok_sensitivity_biserial[, "cor_var"]

# In case you want to visually double check the matrix outside of R
# write.csv(VCV_ESVar_biserial, 'data/outputs/variance-covariance_matrices/VCV_ESVar_biserial.csv')


# STATISTICAL ANALYSIS: 
sensitivy.model.intercept.biserial <- rma.mv(cor,
                                             VCV_ESVar_biserial,
                                             mods = ~ 1,
                                             random = list(~ 1 | StudyID,
                                                           ~ 1 | LaboratoryID,
                                                           ~ 1 | PopulationID,
                                                           ~ 1 | Species,
                                                           ~ 1 | Species_phylo_biserial,
                                                           ~ 1 | EffectID),
                                             method = "REML",
                                             R = list(Species_phylo_biserial = phylo_cor_biserial),
                                             test = "t",
                                             data = meta.final_ok_sensitivity_biserial)

#save(sensitivy.model.intercept.biserial, file = "data/outputs/statistical_models/sensitivy_model_intercept_biserial.RData")
load(file = "data/outputs/statistical_models/sensitivy_model_intercept_biserial.RData") #sensitivy.model.intercept.biserial

# Printing the summary results of the model
print(sensitivy.model.intercept.biserial, digits = 3)

# Model funnel plots
par(mfrow = c(1, 1))
funnel(sensitivy.model.intercept.biserial)

# Estimating heterogeneity as I2 (Nakagawa and Santos 2012)
I2.model_biserial <- i2_ml(sensitivy.model.intercept.biserial)
round(I2.model_biserial, 1)


# TABLE:
results_sensitivy.model.intercept.biserial <- orchaRd::mod_results(sensitivy.model.intercept.biserial, 
                                                                   group = "StudyID", 
                                                                   subset = TRUE)

round(as.data.frame(results_sensitivy.model.intercept.biserial[[1]])[,c(2:6)], 3)


# FIGURE: effect of each group of hormones on fitness traits.
fig_hormones_fitness_intercept.biserial <- orchaRd::orchard_plot(sensitivy.model.intercept.biserial, 
                                                                 group = "StudyID", 
                                                                 xlab = "Effect size",
                                                                 trunk.size = 2,
                                                                 branch.size = 2,
                                                                 twig.size = 1)

# When considering only those studies that were biserial transformed, we get 
# an overall negative effect (same as when we analyze all data together)


##########################
# c) Data base that contains the data base with log transformation

# Calculating the ROM (log transform ratio of means) correlation coefficients 
# and their variance

meta.final_ok_sensitivity_rom <- droplevels(subset(meta.final_ok_ok, 
                                                   meta.final_ok_ok$Effect_size_type == "Mean"))

nrow(meta.final_ok_sensitivity_rom)


db.rom <- escalc(measure = "ROM", 
                 n2i = as.numeric(meta.final_ok_sensitivity_rom$control_n_total),
                 n1i = as.numeric(meta.final_ok_sensitivity_rom$exp_n_total),
                 m2i = as.numeric(meta.final_ok_sensitivity_rom$Control_Mean),
                 m1i = as.numeric(meta.final_ok_sensitivity_rom$Exp_Mean),
                 sd2i = as.numeric(meta.final_ok_sensitivity_rom$sd_control),
                 sd1i = as.numeric(meta.final_ok_sensitivity_rom$sd_exp))

# yi is the correlation, vi the variance
# i will therefore add this information to the general db.

meta.final_ok_sensitivity_rom$cor <- db.rom$yi
meta.final_ok_sensitivity_rom$cor_var <- db.rom$vi

# we had to exclude three values for which the sampling variances were too large
# specifically, 281.16431 and 66.6377, and 2.185199, so that the models could run
# for the latter, we suspect that SD rather than SEM were reported
meta.final_ok_sensitivity_rom <- meta.final_ok_sensitivity_rom[meta.final_ok_sensitivity_rom$cor_var < 2, ]


nrow(meta.final_ok_sensitivity_rom)
length(unique(meta.final_ok_sensitivity_rom$StudyID))


# # data simulation for exploring those large variances
# hist(rnorm(mean=0.019,sd=0.936,n=74),breaks=100)
# hist(rnorm(mean=0.019,sd=0.936,n=37),breaks=100)
# hist(rnorm(mean=0.015,sd=1.033,n=22),breaks=100)
# hist(rnorm(mean=0.198,sd=1.091,n=24),breaks=100)
# hist(rnorm(mean=2.600,sd=11.62755,n=20),breaks=100)
# hist(rnorm(mean=2.700,sd=7.833901,n=17),breaks=100)
# hist(rnorm(mean=0.019,sd=0.936,n=1000),breaks=100)
# hist(rnorm(mean=0.015,sd=1.033,n=1000),breaks=100)
# hist(rnorm(mean=0.198,sd=1.091,n=1000),breaks=100)
# hist(rnorm(mean=2.600,sd=11.62755,n=1000),breaks=100)
# hist(rnorm(mean=2.700,sd=7.833901,n=1000),breaks=100)

# I will now run an intercept-only meta-analytic model. But before doing this
# I need to check the species that are present in this subset. There are 4 species
# less than in the general data base (16 vs 20), and also less entries. I will 
# therefore calculate a new phylogenetic tree and variance-covariance matrix.

# # PHYLOGENETIC TREE:
# resolved_names_rom <- tnrs_match_names(as.character(unique(meta.final_ok_sensitivity_rom$Species)))
# 
# # Saving the taxonomic data created on the 5th February 2025 to speed the
# # process in the future and allow full reproducibility
# save(resolved_names_rom, file = "data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_rom.RData")

# Loading the taxonomic data created on the 5th February 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_rom.RData") #resolved_names_rom

# # extracting phylogenetic information
# my_tree_rom <- tol_induced_subtree(ott_ids =
#                                      resolved_names_rom[,"ott_id"],
#                                    label_format = "name")
# 
# # # Quick tree plotting
# # plot(my_tree_rom, no.margin = TRUE)
# 
# # We need to check for the existence of polytomies
# is.binary(my_tree_rom)
# # Yes, meaning there are no polytomies. Let's go on.
# 
# # To confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# my_tree_rom$tip.label <- gsub("_", " ", my_tree_rom$tip.label)
# 
# intersect(as.character(my_tree_rom$tip.label),
#           as.character(meta.final_ok_sensitivity_rom$Species))
# 
# # Listed in our database but not in the tree
# setdiff(as.character(meta.final_ok_sensitivity_rom$Species),
#         as.character(my_tree_rom$tip.label))
# 
# # Listed in the tree but not in our database
# setdiff(as.character(my_tree_rom$tip.label),
#         as.character(meta.final_ok_sensitivity_rom$Species))
# # No error or inconsistencies found.
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# # 5th February 2025
# save(my_tree_rom, file = "data/outputs/phylogenetic_files/tree_sens.pearson.Rdata")

# We can now load the saved tree
load("data/outputs/phylogenetic_files/tree_sens.pearson.Rdata") #my_tree_rom

# # Compute branch lengths of tree
# phylo_branch_rom <- compute.brlen(my_tree_rom,
#                                   method = "Grafen", power = 1)
# 
# # Check if tree is ultrametric
# is.ultrametric(phylo_branch_rom)
# # TRUE
# 
# # Matrix to be included in the models
# phylo_cor_rom <- vcv(phylo_branch_rom, cor = T)
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(phylo_cor_rom, file = "data/outputs/phylogenetic_files/phylo_cor_rom.Rdata")

# we can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor_rom.Rdata") #phylo_cor_rom

# Creating a duplicate species variable for the phylogenetic analysis
meta.final_ok_sensitivity_rom$Species_phylo_rom <- 
  meta.final_ok_sensitivity_rom$Species


# VARIANCE-COVARIANCE MATRIX:
VCV_ESVar_rom <- matrix(0, nrow = nrow(meta.final_ok_sensitivity_rom), 
                        ncol = nrow(meta.final_ok_sensitivity_rom))

# Names rows and columns for each obsID
rownames(VCV_ESVar_rom) <- meta.final_ok_sensitivity_rom[, "EffectID"]
colnames(VCV_ESVar_rom) <- meta.final_ok_sensitivity_rom[, "EffectID"]

# Finds effect sizes that come from the same study
shared_coord_rom <- which(meta.final_ok_sensitivity_rom[, "StudyID"] %in% 
                            meta.final_ok_sensitivity_rom[duplicated(meta.final_ok_sensitivity_rom[, "StudyID"]), 
                                                          "StudyID"] == TRUE)

combinations_rom <- do.call("rbind", tapply(shared_coord_rom, 
                                            meta.final_ok_sensitivity_rom[shared_coord_rom, "StudyID"], 
                                            function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations_rom)[1]) {
  p1_rom <- combinations_rom[i, 1]
  p2_rom <- combinations_rom[i, 2]
  p1_p2_cov_rom <- 0.5 * sqrt(meta.final_ok_sensitivity_rom[p1_rom, "cor_var"]) * 
    sqrt(meta.final_ok_sensitivity_rom[p2_rom, "cor_var"])
  VCV_ESVar_rom[p1_rom, p2_rom] <- p1_p2_cov_rom
  VCV_ESVar_rom[p2_rom, p1_rom] <- p1_p2_cov_rom
} 


# Enters previously calculated effect size sampling variances into diagonals 
diag(VCV_ESVar_rom) <- meta.final_ok_sensitivity_rom[, "cor_var"]

# In case you want to visually double check the matrix outside of R
# write.csv(VCV_ESVar_rom, 'data/outputs/variance-covariance_matrices/VCV_ESVar_rom.csv')


# STATISTICAL ANALYSIS: 
sensitivy.model.intercept.rom <- rma.mv(cor,
                                        VCV_ESVar_rom,
                                        mods = ~ 1,
                                        random = list(~ 1 | StudyID,
                                                      ~ 1 | LaboratoryID,
                                                      ~ 1 | PopulationID,
                                                      ~ 1 | Species,
                                                      ~ 1 | Species_phylo_rom,
                                                      ~ 1 | EffectID),
                                        method = "REML",
                                        R = list(Species_phylo_rom = phylo_cor_rom),
                                        test = "t",
                                        data = meta.final_ok_sensitivity_rom)

#save(sensitivy.model.intercept.rom, file = "data/outputs/statistical_models/sensitivy_model_intercept_rom.RData")
load(file = "data/outputs/statistical_models/sensitivy_model_intercept_rom.RData") #sensitivy.model.intercept.rom


# Printing the summary results of the model
print(sensitivy.model.intercept.rom, digits = 3)

# Model funnel plots
par(mfrow = c(1, 1))
funnel(sensitivy.model.intercept.rom)

# Estimating heterogeneity as I2 (Nakagawa and Santos 2012)
I2.model_rom <- i2_ml(sensitivy.model.intercept.rom)
round(I2.model_rom, 1)


# TABLE:
results_sensitivy.model.intercept.rom <- orchaRd::mod_results(sensitivy.model.intercept.rom, 
                                                              group = "StudyID", 
                                                              subset = TRUE)

round(as.data.frame(results_sensitivy.model.intercept.rom[[1]])[,c(2:6)], 3)



# FIGURE: effect of each group of hormones on fitness traits.
fig_hormones_fitness_intercept.rom <- orchaRd::orchard_plot(sensitivy.model.intercept.rom, 
                                                            group = "StudyID", 
                                                            xlab = "Effect size",
                                                            trunk.size = 2,
                                                            branch.size = 2,
                                                            twig.size = 1)


###############################################################################
# 8.2. - Sensitivity analysis

# As stated above, I will test the robustness of the results for predictions 
# BH.1.1-BH.1.5 by fitting a phylogenetic multilevel meta-regression with the 
# same random effects structure as used above and hormone type as a moderator. 

# The data base we use to run this analysis is the one that only includes
# offspring data.

# We reported these results in the main text: section 3.2.1.

##############################################################################

# STATISTICAL ANALYSIS:
sensitivy.model.bh1 <- rma.mv(cor,
                              VCV_ESVar_off,
                              mods = ~ Hormone_measured_general - 1,
                              random = list(~ 1 | StudyID,
                                            ~ 1 | LaboratoryID,
                                            ~ 1 | PopulationID,
                                            ~ 1 | Species,
                                            ~ 1 | Species_phylo_off,
                                            ~ 1 | EffectID),
                              method = "REML",
                              R = list(Species_phylo_off = phylo_cor_off),
                              test = "t",
                              data = meta.final_ok_ok_off)

#save(sensitivy.model.bh1, file = "data/outputs/statistical_models/sensitivy_model_bh1.RData")
load(file = "data/outputs/statistical_models/sensitivy_model_bh1.RData") #sensitivy.model.bh1


# Printing the summary results of the model
print(sensitivy.model.bh1, digits = 3)
# predict(sensitivy.model.bh1, digits = 3)

# All hormones have a negative effect on fitness traits, and are not statistically
# significant. 


# Calculate marginal R2 with r2_ml
R2.sensitivy.model.bh1 <- r2_ml(sensitivy.model.bh1)
round(R2.sensitivy.model.bh1 * 100, 1)


# TABLE:
results_sensitivy.model.bh1 <- orchaRd::mod_results(sensitivy.model.bh1, 
                                                    mod = "Hormone_measured_general", 
                                                    group = "StudyID", 
                                                    subset = TRUE)

round(as.data.frame(results_sensitivy.model.bh1[[1]])[,c(2:6)], 3)

# FIGURE: effect of each group of hormones on fitness traits.
fig_hormones_fitness <- orchaRd::orchard_plot(sensitivy.model.bh1, 
                                              mod = "Hormone_measured_general", 
                                              group = "StudyID", 
                                              xlab = "Effect size",
                                              trunk.size = 2,
                                              branch.size = 2,
                                              twig.size = 1)

# When all hormones are fitted in the same model, each group of hormones has a 
# negative effect on fitness traits. These effects are not
# statistically significant.


# Pair-wise comparisons including post-hoc Wald tests
sensitivy.model.bh1.pc <- rma.mv(cor,
                                 VCV_ESVar_off,
                                 mods = ~ Hormone_measured_general,
                                 random = list(~ 1 | StudyID,
                                               ~ 1 | LaboratoryID,
                                               ~ 1 | PopulationID,
                                               ~ 1 | Species,
                                               ~ 1 | Species_phylo_off,
                                               ~ 1 | EffectID),
                                 method = "REML",
                                 R = list(Species_phylo_off = phylo_cor_off),
                                 test = "t",
                                 data = meta.final_ok_ok_off)

#save(sensitivy.model.bh1.pc, file = "data/outputs/statistical_models/sensitivy_model_bh1_pc.RData")
load(file = "data/outputs/statistical_models/sensitivy_model_bh1_pc.RData") #sensitivy.model.bh1.pc

print(sensitivy.model.bh1.pc, digits = 3)

# Post-hoc Wald tests to test for statistically significant differences between 
# all levels

# Information on group comparison between androgens and glucocorticoids and 
# androgens and THs are obtained from the main table. I need to further compare the group

# glucocorticoids vs thyroids:
car::linearHypothesis(sensitivy.model.bh1.pc, rbind(c(0,1,-1))) 

###############################################################################
# 8.3. - Sensitivity analysis

# Using a vector of sampling variances rather than the variance-covariance matrix 
# for i) the intercept only model for both the 
# i) entire data base and ii) only offspring fitness traits, and the iii) 
# the hormone type meta regression for the offspring database.

# Note that this sensitivity analysis was not pre-registered, but done afterwards
# to further confirm the robustness of our results.

##############################################################################

# 8.3.i -  Intercept only model for the entire database
meta.model_corvar <- rma.mv(cor,
                            cor_var,
                            mods = ~ 1,
                            random = list(~ 1 | StudyID,
                                          ~ 1 | LaboratoryID,
                                          ~ 1 | PopulationID,
                                          ~ 1 | Species,
                                          ~ 1 | Species_phylo,
                                          ~ 1 | EffectID),
                            method = "REML",
                            R = list(Species_phylo = phylo_cor),
                            test = "t",
                            data = meta.final_ok_ok)

#save(meta.model, file = "data/outputs/statistical_models/meta.model_corvar.Rdata")
load(file = "data/outputs/statistical_models/meta.model_corvar.Rdata") #meta.model

# Printing the summary results of the model
print(meta.model_corvar, digits = 3)

# Printing the results again, but adding the credibility/prediction interval, 
# which uses the heterogeneity to generate an interval that should contain 95%
# of the effect sizes of any future or unknown studies with similar features
# as those included in the current database
predict(meta.model_corvar, digits = 3)

# Model funnel plots
par(mfrow = c(1, 1))
funnel(meta.model_corvar)

# FIGURE: overall effect of yolk hormones on fitness
figure_meta.model_corvar <- orchaRd::orchard_plot(meta.model_corvar, mod = "1", 
                                                  group = "StudyID", 
                                                  xlab = "Effect size",
                                                  transfm = "none",
                                                  trunk.size = 2,
                                                  branch.size = 2,
                                                  twig.size = 1)



# 8.3.ii -  Intercept only model for the offspring database
meta.model.off_corvar <- rma.mv(cor,
                                cor_var,
                                mods = ~ 1,
                                random = list(~ 1 | StudyID,
                                              ~ 1 | LaboratoryID,
                                              ~ 1 | PopulationID,
                                              ~ 1 | Species,
                                              ~ 1 | Species_phylo_off,
                                              ~ 1 | EffectID),
                                method = "REML",
                                R = list(Species_phylo_off = phylo_cor_off),
                                test = "t",
                                data = meta.final_ok_ok_off)

#save(meta.model.off_corvar, file = "data/outputs/statistical_models/meta.model.off_corvar.Rdata")
load(file = "data/outputs/statistical_models/meta.model.off_corvar.Rdata") #meta.model.off

# Printing the summary results of the model
print(meta.model.off_corvar, digits = 3)

# Printing the results again, but adding the credibility/prediction interval, 
# which uses the heterogeneity to generate an interval that should contain 95%
# of the effect sizes of any future or unknown studies with similar features
# as those included in the current database
predict(meta.model.off_corvar, digits = 3)

# Model funnel plots
par(mfrow = c(1, 1))
funnel(meta.model.off_corvar)

# FIGURE: overall effect of yolk hormones on fitness
figure_meta.model.off_corvar <- orchaRd::orchard_plot(meta.model.off_corvar, mod = "1", 
                                                      group = "StudyID", 
                                                      xlab = "Effect size",
                                                      transfm = "none",
                                                      trunk.size = 2,
                                                      branch.size = 2,
                                                      twig.size = 1)


# 8.3.iii -  Hormone type meta regression for the offspring database.
meta.model.off_perhorm_corvar <- rma.mv(cor,
                                        cor_var,
                                        mods = ~ Hormone_measured_general - 1,
                                        random = list(~ 1 | StudyID,
                                                      ~ 1 | LaboratoryID,
                                                      ~ 1 | PopulationID,
                                                      ~ 1 | Species,
                                                      ~ 1 | Species_phylo_off,
                                                      ~ 1 | EffectID),
                                        method = "REML",
                                        R = list(Species_phylo_off = phylo_cor_off),
                                        test = "t",
                                        data = meta.final_ok_ok_off)

#save(sensitivy.model.bh1, file = "data/outputs/statistical_models/meta.model.off_perhorm_corvar.RData")
load(file = "data/outputs/statistical_models/meta.model.off_perhorm_corvar.RData") #sensitivy.model.bh1


# Printing the summary results of the model
print(meta.model.off_perhorm_corvar, digits = 3)
#predict(meta.model.off_perhorm_corvar, digits = 3)

# All hormones have a negative effect on fitness traits, and are not statistically
# significant. 


# Calculate marginal R2 with r2_ml
R2.meta.model.off_perhorm_corvar <- r2_ml(meta.model.off_perhorm_corvar)
round(R2.meta.model.off_perhorm_corvar * 100, 1)


# TABLE:
results_meta.model.off_perhorm_corvar <- orchaRd::mod_results(meta.model.off_perhorm_corvar, 
                                                              mod = "Hormone_measured_general", 
                                                              group = "StudyID", 
                                                              subset = TRUE)

round(as.data.frame(results_meta.model.off_perhorm_corvar[[1]])[,c(2:6)], 3)

# FIGURE: effect of each group of hormones on fitness traits.
fig_meta.model.off_perhorm_corvar <- orchaRd::orchard_plot(meta.model.off_perhorm_corvar, 
                                                           mod = "Hormone_measured_general", 
                                                           group = "StudyID", 
                                                           xlab = "Effect size",
                                                           trunk.size = 2,
                                                           branch.size = 2,
                                                           twig.size = 1)



###############################################################################
# 9) BIOLOGICAL EXPLORATORY HYPOTHESIS

# BEH.1 – Associations between each type of maternal egg hormone group and 
# maternal fitness

# Data base: the one that includes information only for the mother.

# We cannot test this exploratory hypothesis for thyroid hormones since we only
# have 3 effect sizes for these hormones (which is smaller than 5, the number
# we establishes in our pre-registration). Hence, we will test the maternal 
# fitness consequences of androgens and glucocorticoids.
###############################################################################
meta.final_ok_ok_mom_and_gcs <- droplevels(subset(meta.final_ok_ok_mom, 
                                                  meta.final_ok_ok_mom$Hormone_measured_general != "TH3 and TH4"))

# # PHYLOGENETIC TREE:
# resolved_names_mom_and_gcs <- tnrs_match_names(as.character(unique(meta.final_ok_ok_mom_and_gcs$Species)))
# 
# # Saving the taxonomic data created on the 5th February 2025 to speed the
# # process in the future and allow full reproducibility
# save(resolved_names_mom_and_gcs, file = "data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_mom_and_gcs.RData")

# Loading the taxonomic data created on the 5th February 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_mom_and_gcs.RData") #resolved_names_mom_and_gcs

# # extracting phylogenetic information
# my_tree_mom_and_gcs <- tol_induced_subtree(ott_ids =
#                                      resolved_names_mom_and_gcs[,"ott_id"],
#                                     label_format = "name")
# # # Quick tree plotting
# # plot(my_tree_mom_and_gcs, no.margin = TRUE)
# 
# # We need to check for the existence of polytomies
# is.binary(my_tree_mom_and_gcs)
# # Yes, meaning there are no polytomies. Let's go on.
# 
# # To confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# my_tree_mom_and_gcs$tip.label <- gsub("_", " ", my_tree_mom_and_gcs$tip.label)
# 
# intersect(as.character(my_tree_mom_and_gcs$tip.label),
#           as.character(meta.final_ok_ok_mom_and_gcs$Species))
# 
# # Listed in our database but not in the tree
# setdiff(as.character(meta.final_ok_ok_mom_and_gcs$Species),
#         as.character(my_tree_mom_and_gcs$tip.label))
# 
# # Listed in the tree but not in our database
# setdiff(as.character(my_tree_mom_and_gcs$tip.label),
#         as.character(meta.final_ok_ok_mom_and_gcs$Species))
# # No error or inconsistencies found.
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# # 5th February 2025
# save(my_tree_mom_and_gcs, file = "data/outputs/phylogenetic_files/tree_mom_and_gcs.Rdata")

# We can now load the saved tree
load("data/outputs/phylogenetic_files/tree_mom_and_gcs.Rdata") #my_tree_mom_and_gcs

# # Compute branch lengths of tree
# phylo_branch_mom_and_gcs <- compute.brlen(my_tree_mom_and_gcs, method = "Grafen", power = 1)
# 
# # Check if tree is ultrametric
# is.ultrametric(phylo_branch_mom_and_gcs)
# # TRUE
# 
# # Matrix to be included in the models
# phylo_cor_mom_and_gcs <- vcv(phylo_branch_mom_and_gcs, cor = T)
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(phylo_cor_mom_and_gcs, file = "data/outputs/phylogenetic_files/phylo_cor_mom_and_gcs.Rdata")

# we can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor_mom_and_gcs.Rdata") #phylo_cor_mom_and_gcs

# Creating a duplicate species variable for the phylogenetic analysis
meta.final_ok_ok_mom_and_gcs$Species_phylo_mom_and_gcs <- meta.final_ok_ok_mom_and_gcs$Species


# VARIANCE-COVARIANCE MATRIX:
VCV_ESVar_mom_and_gcs <- matrix(0, nrow = nrow(meta.final_ok_ok_mom_and_gcs), 
                                ncol = nrow(meta.final_ok_ok_mom_and_gcs))

# Names rows and columns for each obsID
rownames(VCV_ESVar_mom_and_gcs) <- meta.final_ok_ok_mom_and_gcs[, "EffectID"]
colnames(VCV_ESVar_mom_and_gcs) <- meta.final_ok_ok_mom_and_gcs[, "EffectID"]

# Finds effect sizes that come from the same study
shared_coord_mom_and_gcs <- which(meta.final_ok_ok_mom_and_gcs[, "StudyID"] %in% 
                                    meta.final_ok_ok_mom_and_gcs[duplicated(meta.final_ok_ok_mom_and_gcs[, "StudyID"]), 
                                                                 "StudyID"] == TRUE)

combinations_mom_and_gcs <- do.call("rbind", tapply(shared_coord_mom_and_gcs, 
                                                    meta.final_ok_ok_mom_and_gcs[shared_coord_mom_and_gcs, "StudyID"], 
                                                    function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations_mom_and_gcs)[1]) {
  p1_mom_and_gcs <- combinations_mom_and_gcs[i, 1]
  p2_mom_and_gcs <- combinations_mom_and_gcs[i, 2]
  p1_p2_cov_mom_and_gcs <- 0.5 * sqrt(meta.final_ok_ok_mom_and_gcs[p1_mom_and_gcs, "cor_var"]) * 
    sqrt(meta.final_ok_ok_mom_and_gcs[p2_mom_and_gcs, "cor_var"])
  VCV_ESVar_mom_and_gcs[p1_mom_and_gcs, p2_mom_and_gcs] <- p1_p2_cov_mom_and_gcs
  VCV_ESVar_mom_and_gcs[p2_mom_and_gcs, p1_mom_and_gcs] <- p1_p2_cov_mom_and_gcs
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(VCV_ESVar_mom_and_gcs) <- meta.final_ok_ok_mom_and_gcs[, "cor_var"]

# # In case you want to visually double check the matrix outside of R
# write.csv(VCV_ESVar_mom_and_gcs, 'data/outputs/variance-covariance_matrices/VCV_ESVar_mom_and_gcs.csv')



# STATISTICAL ANALYSIS:
meta.regression.beh1 <- rma.mv(cor,
                               VCV_ESVar_mom_and_gcs,
                               mods = ~ Hormone_measured_general - 1,
                               random = list(~ 1 | StudyID,
                                             ~ 1 | LaboratoryID,
                                             ~ 1 | PopulationID,
                                             ~ 1 | Species,
                                             ~ 1 | Species_phylo_mom_and_gcs,
                                             ~ 1 | EffectID),
                               method = "REML",
                               R = list(Species_phylo_mom_and_gcs = phylo_cor_mom_and_gcs),
                               test = "t",
                               data = meta.final_ok_ok_mom_and_gcs)

#save(meta.regression.beh1, file = "data/outputs/statistical_models/meta_regression_beh1.RData")
load(file = "data/outputs/statistical_models/meta_regression_beh1.RData") #meta.regression.beh1


# Printing the summary results of the model
print(meta.regression.beh1, digits = 3)

# Calculate marginal R2 with r2_ml
R2.method_beh1 <- r2_ml(meta.regression.beh1)
round(R2.method_beh1 * 100, 1)



# TABLE:
results_sensitivy.model.beh1 <- orchaRd::mod_results(meta.regression.beh1, 
                                                     mod = "Hormone_measured_general", 
                                                     group = "StudyID", 
                                                     subset = TRUE)

round(as.data.frame(results_sensitivy.model.beh1[[1]])[,c(2:6)], 3)


# FIGURE: effect of each group of hormones on maternal fitness traits.
fig_hormones_beh1 <- orchaRd::orchard_plot(meta.regression.beh1, 
                                           mod = "Hormone_measured_general", 
                                           group = "StudyID", 
                                           xlab = "Effect size",
                                           trunk.size = 2,
                                           branch.size = 2,
                                           twig.size = 1)

# Negative effect androgens, positive effect for cort. All effect sizes
# are quite small and not statistically significant effect.


# Pair-wise comparisons including post-hoc Wald tests
meta.regression.beh1_pc <- rma.mv(cor,
                                  VCV_ESVar_mom_and_gcs,
                                  mods = ~ Hormone_measured_general,
                                  random = list(~ 1 | StudyID,
                                                ~ 1 | LaboratoryID,
                                                ~ 1 | PopulationID,
                                                ~ 1 | Species,
                                                ~ 1 | Species_phylo_mom_and_gcs,
                                                ~ 1 | EffectID),
                                  method = "REML",
                                  R = list(Species_phylo_mom_and_gcs = phylo_cor_mom_and_gcs),
                                  test = "t",
                                  data = meta.final_ok_ok_mom_and_gcs)

#save(meta.regression.beh1_pc, file = "data/outputs/statistical_models/meta_regression_beh1_pc.RData")
load(file = "data/outputs/statistical_models/meta_regression_beh1_pc.RData") #meta.regression.beh1_pc

print(meta.regression.beh1_pc, digits = 3)

################################################################################
# BEH.2 – Associations between each type of maternal hormone group and offspring 
# fitness proxies, regarding whether hormones are studied/manipulated in yolk vs.
# albumen.

# Data base: offspring.

################################################################################

# In the current data base, we have the information on where the hormones
# were measured/manipulated divided in two columns. One for correlational
# studies (variable name 'Site_hormone_measured') and one for experimental
# studies (variable name 'Site_egg_exp'). To test this prediction, I need to 
# combine both columns into a new one.

# I create this variable and include firs the data for correlational studies
meta.final_ok_ok_off$Site_measured <- meta.final_ok_ok_off$Site_horm_measured

# I now complete those 'NAs' with the information from experimental studies
meta.final_ok_ok_off$Site_measured <- as.character(meta.final_ok_ok_off$Site_measured)
meta.final_ok_ok_off$Site_egg_exp <- as.character(meta.final_ok_ok_off$Site_egg_exp)

nrow_na_site_horm <- which(is.na(meta.final_ok_ok_off$Site_measured)) 
meta.final_ok_ok_off$Site_measured[nrow_na_site_horm] <- meta.final_ok_ok_off$Site_egg_exp[nrow_na_site_horm]

meta.final_ok_ok_off$Site_measured <- as.factor(meta.final_ok_ok_off$Site_measured)


# I now check that I have enough levels for running the model as we wrote down
# in the pre-registration.
table(meta.final_ok_ok_off$Hormone_measured_general, 
      meta.final_ok_ok_off$Site_measured)


# I cannot test this hypothesis for androgens and THs because we do not have
# enough levels. For androgens, we only have 2 rows that contain data on the
# effect of these hormones when measured in the albumen. This is a smaller number 
# than the number we stated in our pre-registration as the minimum number of 
# levels required (i.e., 5). For THs, there are no papers measuring these
# hormones in the albumen.

# For corticosterone, we do have information. 19 papers measured this hormone
# in the albumen, and there are 8 papers that injected cort in the albumen but
# indicate that it migrates to the yolk. Hence, following what the authors wrote
# down, it should be considered as a hormone measured in the yolk. I will therefore
# rename 'albumen (but they say that it migrates to the yolk)' as 'yolk'.

meta.final_ok_ok_off$Site_measured[meta.final_ok_ok_off$Site_measured == 
                                     "albumen (but they say that it migrates to the yolk)"] <- "yolk"

# Good. Now, I will subset the data base to work with those rows that only
# contain information on corticosterone (i.e., the hormone for which we have
# enough data for testing our prediction) and fitness of the offspring.

meta.final_ok_ok_off_site_cort <- droplevels(subset(meta.final_ok_ok_off, 
                                                    meta.final_ok_ok_off$Hormone_measured_general == "corticosterone"))

# I now check the levels.
table(meta.final_ok_ok_off_site_cort$Hormone_measured_general,
      meta.final_ok_ok_off_site_cort$Site_measured)
# Good


# Now, there are some rows for which we have NAs. Hence, I will remove them 
# so I can run the phylogenetic tree, the variance-covariance matrix and the 
# statistical model with the appropiate data base.
meta.final_ok_ok_off_site_cort <- droplevels(subset(meta.final_ok_ok_off_site_cort, 
                                                    meta.final_ok_ok_off_site_cort$Site_measured != "NA"))
nrow(meta.final_ok_ok_off_site_cort)


# # PHYLOGENETIC TREE:
# resolved_names_off_site_cort <- tnrs_match_names(as.character(unique(meta.final_ok_ok_off_site_cort$Species)))
# 
# # Saving the taxonomic data created on the 5th February 2025 to speed the
# # process in the future and allow full reproducibility
# save(resolved_names_off_site_cort, file = "data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_off_site_cort.RData")

# Loading the taxonomic data created on the 5th February 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_off_site_cort.RData") #resolved_names_off_site_cort

# # extracting phylogenetic information
# my_tree_off_site_cort <- tol_induced_subtree(ott_ids =
#                                              resolved_names_off_site_cort[,"ott_id"],
#                                              label_format = "name")
# # # Quick tree plotting
# # plot(my_tree_off_site_cort, no.margin = TRUE)
# 
# # We need to check for the existence of polytomies
# is.binary(my_tree_off_site_cort)
# # Yes, meaning there are no polytomies. Let's go on.
# 
# # To confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# my_tree_off_site_cort$tip.label <- gsub("_", " ", my_tree_off_site_cort$tip.label)
# 
# intersect(as.character(my_tree_off_site_cort$tip.label),
#           as.character(meta.final_ok_ok_off_site_cort$Species))
# 
# # Listed in our database but not in the tree
# setdiff(as.character(meta.final_ok_ok_off_site_cort$Species),
#         as.character(my_tree_off_site_cort$tip.label))
# 
# # Listed in the tree but not in our database
# setdiff(as.character(my_tree_off_site_cort$tip.label),
#         as.character(meta.final_ok_ok_off_site_cort$Species))
# # No error or inconsistencies found.
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# # 5th February 2025
# save(my_tree_off_site_cort, file = "data/outputs/phylogenetic_files/tree_off_site_cort.Rdata")

# We can now load the saved tree
load("data/outputs/phylogenetic_files/tree_off_site_cort.Rdata") #my_tree_off_site_cort

# # Compute branch lengths of tree
# phylo_branch_off_site_cort <- compute.brlen(my_tree_off_site_cort, method = "Grafen", power = 1)
# 
# # Check if tree is ultrametric
# is.ultrametric(phylo_branch_off_site_cort)
# # TRUE
# 
# # Matrix to be included in the models
# phylo_cor_off_site_cort <- vcv(phylo_branch_off_site_cort, cor = T)
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(phylo_cor_off_site_cort, file = "data/outputs/phylogenetic_files/phylo_cor_off_site_cort.Rdata")

# we can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor_off_site_cort.Rdata") #phylo_cor_off_site_cort

# Creating a duplicate species variable for the phylogenetic analysis
meta.final_ok_ok_off_site_cort$Species_phylo_off_site_cort <- meta.final_ok_ok_off_site_cort$Species

# VARIANCE-COVARIANCE MATRIX:
VCV_ESVar_off_site_cort <- matrix(0, nrow = nrow(meta.final_ok_ok_off_site_cort), 
                                  ncol = nrow(meta.final_ok_ok_off_site_cort))

# Names rows and columns for each obsID
rownames(VCV_ESVar_off_site_cort) <- meta.final_ok_ok_off_site_cort[, "EffectID"]
colnames(VCV_ESVar_off_site_cort) <- meta.final_ok_ok_off_site_cort[, "EffectID"]

# Finds effect sizes that come from the same study
shared_coord_off_site_cort <- which(meta.final_ok_ok_off_site_cort[, "StudyID"] %in% 
                                      meta.final_ok_ok_off_site_cort[duplicated(meta.final_ok_ok_off_site_cort[, "StudyID"]), 
                                                                     "StudyID"] == TRUE)

combinations_off_site_cort <- do.call("rbind", tapply(shared_coord_off_site_cort, 
                                                      meta.final_ok_ok_off_site_cort[shared_coord_off_site_cort, "StudyID"], 
                                                      function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations_off_site_cort)[1]) {
  p1_off_site_cort <- combinations_off_site_cort[i, 1]
  p2_off_site_cort <- combinations_off_site_cort[i, 2]
  p1_p2_cov_off_site_cort <- 0.5 * sqrt(meta.final_ok_ok_off_site_cort[p1_off_site_cort, "cor_var"]) * 
    sqrt(meta.final_ok_ok_off_site_cort[p2_off_site_cort, "cor_var"])
  VCV_ESVar_off_site_cort[p1_off_site_cort, p2_off_site_cort] <- p1_p2_cov_off_site_cort
  VCV_ESVar_off_site_cort[p2_off_site_cort, p1_off_site_cort] <- p1_p2_cov_off_site_cort
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(VCV_ESVar_off_site_cort) <- meta.final_ok_ok_off_site_cort[, "cor_var"]

# # In case you want to visually double check the matrix outside of R
# write.csv(VCV_ESVar_off_site_cort, 'data/outputs/variance-covariance_matrices/VCV_ESVar_off_site_cort.csv')


# STATISTICAL ANALYSIS:
meta.regression.beh2 <- rma.mv(cor,
                               VCV_ESVar_off_site_cort,
                               mods = ~ Site_measured - 1,
                               random = list(~ 1 | StudyID,
                                             ~ 1 | LaboratoryID,
                                             ~ 1 | PopulationID,
                                             ~ 1 | Species,
                                             ~ 1 | Species_phylo_off_site_cort,
                                             ~ 1 | EffectID),
                               method = "REML",
                               R = list(Species_phylo_off_site_cort = phylo_cor_off_site_cort),
                               test = "t",
                               data = meta.final_ok_ok_off_site_cort)

#save(meta.regression.beh2, file = "data/outputs/statistical_models/meta_regression_beh2.RData")
load(file = "data/outputs/statistical_models/meta_regression_beh2.RData") #meta.regression.beh2

# Printing the summary results of the model
print(meta.regression.beh2, digits = 3)

# Calculate marginal R2 with r2_ml
R2.method_beh2 <- r2_ml(meta.regression.beh2)
round(R2.method_beh2 * 100, 1)



# TABLE:
results_hormones_beh2 <- orchaRd::mod_results(meta.regression.beh2, 
                                              mod = "Site_measured", 
                                              group = "StudyID", 
                                              subset = TRUE)

round(as.data.frame(results_hormones_beh2[[1]])[,c(2:6)], 3)


# FIGURE: effect of each group of hormones on maternal fitness traits.
fig_hormones_beh2 <- orchaRd::orchard_plot(meta.regression.beh2, 
                                           mod = "Site_measured", 
                                           group = "StudyID", 
                                           xlab = "Effect size",
                                           trunk.size = 2,
                                           branch.size = 2,
                                           twig.size = 1)

# Positive effect of corticosterone when measured both in the yolk and
# in the albumen . None of these effects are statistically significant.



# Pair-wise comparisons including post-hoc Wald tests
meta.regression.beh2_pc <- rma.mv(cor,
                                  VCV_ESVar_off_site_cort,
                                  mods = ~ Site_measured,
                                  random = list(~ 1 | StudyID,
                                                ~ 1 | LaboratoryID,
                                                ~ 1 | PopulationID,
                                                ~ 1 | Species,
                                                ~ 1 | Species_phylo_off_site_cort,
                                                ~ 1 | EffectID),
                                  method = "REML",
                                  R = list(Species_phylo_off_site_cort = phylo_cor_off_site_cort),
                                  test = "t",
                                  data = meta.final_ok_ok_off_site_cort)

#save(meta.regression.beh2_pc, file = "data/outputs/statistical_models/meta_regression_beh2_pc.RData")
load(file = "data/outputs/statistical_models/meta_regression_beh2_pc.RData") #meta.regression.beh2_pc

print(meta.regression.beh2_pc, digits = 3)

################################################################################
# 10) METHODOLOGICAL EXPLANATORY HYPOTHESIS 

# MEH.1 - The association between maternal yolk hormones and maternal and offspring 
# fitness proxies vary depending on the sampling technique used (i.e., egg biopsy
# vs. entire egg collected). 

# Data base: general data base (i.e., the one that includes information on both
# maternal and offspring fitness)

################################################################################
# I first check if we have enough levels to test this hypothesis.
table(meta.final_ok_ok$Hormone_measured_general, 
      meta.final_ok_ok$Egg_sampl_method)

# I only have information for studies removing entire eggs or doing biopsies for
# androgens. For the other groups of hormones I do not have this information.
# Therefore, I cannot test this methodological hypothesis. I will therefore 
# work only with the data base for androgens.

meta.final_ok_ok_androgens <- droplevels(subset(meta.final_ok_ok, 
                                                meta.final_ok_ok$Hormone_measured_general == "androgens"))

# I now remove those rows that contain no information on how the egg was sampled.
meta.final_ok_ok_androgens_method <- droplevels(subset(meta.final_ok_ok_androgens,
                                                       meta.final_ok_ok_androgens$Egg_sampl_method != "NA"))

nrow(meta.final_ok_ok_androgens_method)



# # PHYLOGENETIC TREE:
# resolved_names_androgens_method <- tnrs_match_names(as.character(unique(meta.final_ok_ok_androgens_method$Species)))
# 
# # Saving the taxonomic data created on the 5th February 2025 to speed the
# # process in the future and allow full reproducibility
# save(resolved_names_androgens_method, file = "data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_androgens_method.RData")

# Loading the taxonomic data created on the 5th February 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_androgens_method.RData") #resolved_names_androgens_method

# # extracting phylogenetic information
# my_tree_androgens_method <- tol_induced_subtree(ott_ids =
#                                                   resolved_names_androgens_method[,"ott_id"],
#                                                 label_format = "name")
# # # Quick tree plotting
# # plot(my_tree_androgens_method, no.margin = TRUE)
# 
# # We need to check for the existence of polytomies
# is.binary(my_tree_androgens_method)
# # Yes, meaning there are no polytomies. Let's go on.
# 
# # To confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# my_tree_androgens_method$tip.label <- gsub("_", " ", my_tree_androgens_method$tip.label)
# 
# intersect(as.character(my_tree_androgens_method$tip.label),
#           as.character(meta.final_ok_ok_androgens_method$Species))
# 
# # Listed in our database but not in the tree
# setdiff(as.character(meta.final_ok_ok_androgens_method$Species),
#         as.character(my_tree_androgens_method$tip.label))
# 
# # Listed in the tree but not in our database
# setdiff(as.character(my_tree_androgens_method$tip.label),
#         as.character(meta.final_ok_ok_androgens_method$Species))
# # No error or inconsistencies found.
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# # 5th February 2025
# save(my_tree_androgens_method, file = "data/outputs/phylogenetic_files/tree_androgens_method.Rdata")

# We can now load the saved tree
load("data/outputs/phylogenetic_files/tree_androgens_method.Rdata") #my_tree_androgens_method

# # Compute branch lengths of tree
# phylo_branch_androgens_method <- compute.brlen(my_tree_androgens_method, method = "Grafen", power = 1)
# 
# # Check if tree is ultrametric
# is.ultrametric(phylo_branch_androgens_method)
# # TRUE
# 
# # Matrix to be included in the models
# phylo_cor_androgens_method <- vcv(phylo_branch_androgens_method, cor = T)
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(phylo_cor_androgens_method, file = "data/outputs/phylogenetic_files/phylo_cor_androgens_method.Rdata")

# we can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor_androgens_method.Rdata") #phylo_cor_androgens_method

# Creating a duplicate species variable for the phylogenetic analysis
meta.final_ok_ok_androgens_method$Species_phylo_androgens_method <- meta.final_ok_ok_androgens_method$Species

# VARIANCE-COVARIANCE MATRIX:
VCV_ESVar_androgens_method <- matrix(0, nrow = nrow(meta.final_ok_ok_androgens_method), 
                                     ncol = nrow(meta.final_ok_ok_androgens_method))

# Names rows and columns for each obsID
rownames(VCV_ESVar_androgens_method) <- meta.final_ok_ok_androgens_method[, "EffectID"]
colnames(VCV_ESVar_androgens_method) <- meta.final_ok_ok_androgens_method[, "EffectID"]

# Finds effect sizes that come from the same study
shared_coord_androgens_method <- which(meta.final_ok_ok_androgens_method[, "StudyID"] %in% 
                                         meta.final_ok_ok_androgens_method[duplicated(meta.final_ok_ok_androgens_method[, "StudyID"]), 
                                                                           "StudyID"] == TRUE)

combinations_androgens_method <- do.call("rbind", tapply(shared_coord_androgens_method, 
                                                         meta.final_ok_ok_androgens_method[shared_coord_androgens_method, "StudyID"], 
                                                         function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations_androgens_method)[1]) {
  p1_androgens_method <- combinations_androgens_method[i, 1]
  p2_androgens_method <- combinations_androgens_method[i, 2]
  p1_p2_cov_androgens_method <- 0.5 * sqrt(meta.final_ok_ok_androgens_method[p1_androgens_method, "cor_var"]) * 
    sqrt(meta.final_ok_ok_androgens_method[p2_androgens_method, "cor_var"])
  VCV_ESVar_androgens_method[p1_androgens_method, p2_androgens_method] <- p1_p2_cov_androgens_method
  VCV_ESVar_androgens_method[p2_androgens_method, p1_androgens_method] <- p1_p2_cov_androgens_method
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(VCV_ESVar_androgens_method) <- meta.final_ok_ok_androgens_method[, "cor_var"]

# # In case you want to visually double check the matrix outside of R
# write.csv(VCV_ESVar_androgens_method, 'data/outputs/variance-covariance_matrices/VCV_ESVar_androgens_method.csv')

# STATISTICAL ANALYSIS:
meta.regression.meh1 <- rma.mv(cor,
                               VCV_ESVar_androgens_method,
                               mods = ~ Egg_sampl_method - 1,
                               random = list(~ 1 | StudyID,
                                             ~ 1 | LaboratoryID,
                                             ~ 1 | PopulationID,
                                             ~ 1 | Species,
                                             ~ 1 | Species_phylo_androgens_method,
                                             ~ 1 | EffectID),
                               method = "REML",
                               R = list(Species_phylo_androgens_method = phylo_cor_androgens_method),
                               test = "t",
                               control=list(rel.tol=1e-8),
                               data = meta.final_ok_ok_androgens_method)

#save(meta.regression.meh1, file = "data/outputs/statistical_models/meta_regression_meh1.RData")
load(file = "data/outputs/statistical_models/meta_regression_meh1.RData") #meta.regression.meh1

# Printing the summary results of the model
print(meta.regression.meh1, digits = 3)

# Calculate marginal R2 with r2_ml
R2.method_meh1 <- r2_ml(meta.regression.meh1) 
round(R2.method_meh1 * 100, 1)



# TABLE:
results_hormones_meh1 <- orchaRd::mod_results(meta.regression.meh1, 
                                              mod = "Egg_sampl_method", 
                                              group = "StudyID", 
                                              subset = TRUE) 


round(as.data.frame(results_hormones_meh1[[1]])[,c(2:6)],3 )


# FIGURE: effect of each group of hormones on maternal fitness traits.
fig_hormones_meh1 <- orchaRd::orchard_plot(meta.regression.meh1, 
                                           mod = "Egg_sampl_method", 
                                           group = "StudyID", 
                                           xlab = "Effect size",
                                           trunk.size = 2,
                                           branch.size = 2,
                                           twig.size = 1)

# Negative effect of androgens on maternal and offspring fitness if researchers
# take a biopsy, but positive effect if they removed the entire egg. 
# Effects are not statistically significant. 



# Pair-wise comparisons including post-hoc Wald tests
meta.regression.meh1_pc <- rma.mv(cor,
                                  VCV_ESVar_androgens_method,
                                  mods = ~ Egg_sampl_method,
                                  random = list(~ 1 | StudyID,
                                                ~ 1 | LaboratoryID,
                                                ~ 1 | PopulationID,
                                                ~ 1 | Species,
                                                ~ 1 | Species_phylo_androgens_method,
                                                ~ 1 | EffectID),
                                  method = "REML",
                                  R = list(Species_phylo_androgens_method = phylo_cor_androgens_method),
                                  test = "t",
                                  data = meta.final_ok_ok_androgens_method)

#save(meta.regression.meh1_pc, file = "data/outputs/statistical_models/meta_regression_meh1_pc.RData")
load(file = "data/outputs/statistical_models/meta_regression_meh1_pc.RData") #meta.regression.meh1_pc

print(meta.regression.meh1_pc, digits = 3)


################################################################################
# 11) PUBLICATION BIAS TESTS

# Following Nakagawa et al. 2022 general recommendations:
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13724

# Data base: general data base (i.e., includes maternal and offspring fitness)
################################################################################

# PHYLOGENETIC TREE - Species names:
# Loading the taxonomic data created on the 5th February 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life.RData") #resolved_names

# We can now load the saved tree
load("data/outputs/phylogenetic_files/tree.Rdata") #my_tree

# we can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor.Rdata") #phylo_cor

# VARIANCE-COVARIANCE MATRIX: 
# Matrix already created above for the main model: VCV_ESVar


################################################################################
# Small-study bias
################################################################################

# We follow the general guidelines of Nakagawa et al. 2022, MEE (but be aware of 
# the published corrigendum!). Instead of calculating effective sample size as 
# (Nr + Nu)/(Nr * Nu), which we could only do for part of the dataset (e.g. not 
# for the correlations), we will simply use the the inverse of the sample size
# and its sqrt. As a sensitivity test, we will use n/4 rather than simply n, but 
# results are expected to be very similar (not the slope though!).


meta.final_ok_ok$inv_esz <-  1/meta.final_ok_ok$final_n
meta.final_ok_ok$inv_esz.4 <-  1/(meta.final_ok_ok$final_n/4)
meta.final_ok_ok$sqrt_inv_esz  <-  with(meta.final_ok_ok, sqrt(inv_esz))
meta.final_ok_ok$sqrt_inv_esz.4 <-  with(meta.final_ok_ok, sqrt(inv_esz.4))

# STATISTICAL ANALYSIS:
meta.regression.small.study.effects <- rma.mv(cor,
                                              VCV_ESVar,
                                              mods = ~ 1 + sqrt_inv_esz,
                                              random = list(~ 1 | StudyID,
                                                            ~ 1 | LaboratoryID,
                                                            ~ 1 | PopulationID,
                                                            ~ 1 | Species,
                                                            ~ 1 | Species_phylo,
                                                            ~ 1 | EffectID),
                                              method = "REML",
                                              R = list(Species_phylo = phylo_cor),
                                              test = "t",
                                              data = meta.final_ok_ok)

#save(meta.regression.small.study.effects, file = "data/outputs/statistical_models/meta_regression_small_study_effects.RData")
load(file = "data/outputs/statistical_models/meta_regression_small_study_effects.RData") #meta.regression.small.study.effects

# Printing the summary results of the model
print(meta.regression.small.study.effects, digits = 3)
#predict(meta.regression.small.study.effects, digits = 3)

# Calculate marginal R2 with r2_ml
R2.meta.regression.small.study.effects <- r2_ml(meta.regression.small.study.effects) 
round(R2.meta.regression.small.study.effects * 100, 1)


# TABLE WITH RESULTS FOR A CONTINUOUS VARIABLE:
results_meta.regression.small.study.effects <- orchaRd::mod_results(meta.regression.small.study.effects,
                                                                    mod = "sqrt_inv_esz", 
                                                                    group = "StudyID",
                                                                    weights = "prop")

# FIGURE WITH A CONTINOUS VARIABLE:
figure_meta.regression.small.study.effects <- orchaRd::bubble_plot(meta.regression.small.study.effects, 
                                                                   mod = "sqrt_inv_esz", 
                                                                   group = "StudyID",
                                                                   xlab = "sqrt_inv_esz",
                                                                   legend.pos = "bottom.right")



# STATISTICAL ANALYSIS: sensitivity test using n/4 rather than simply n
# which confirms the results. Not reported in the manuscript or the supplements
meta.regression.small.study.effects.n4 <- rma.mv(cor,
                                                 VCV_ESVar,
                                                 mods = ~ 1 + sqrt_inv_esz.4,
                                                 random = list(~ 1 | StudyID,
                                                               ~ 1 | LaboratoryID,
                                                               ~ 1 | PopulationID,
                                                               ~ 1 | Species,
                                                               ~ 1 | Species_phylo,
                                                               ~ 1 | EffectID),
                                                 method = "REML",
                                                 R = list(Species_phylo = phylo_cor),
                                                 test = "t",
                                                 data = meta.final_ok_ok)
#save(meta.regression.small.study.effects.n4, file = "data/outputs/statistical_models/meta_regression_small_study_effects_n4.RData")
load(file = "data/outputs/statistical_models/meta_regression_small_study_effects_n4.RData") #meta.regression.small.study.effects.n4

# Printing the summary results of the model
print(meta.regression.small.study.effects.n4, digits = 3)

# Calculate marginal R2 with r2_ml
R2.meta.regression.small.study.effects.n4 <- r2_ml(meta.regression.small.study.effects.n4) 
round(R2.meta.regression.small.study.effects.n4 * 100, 1)


# TABLE WITH RESULTS FOR A CONTINUOUS VARIABLE:
results_meta.regression.small.study.effects.n4 <- orchaRd::mod_results(meta.regression.small.study.effects.n4,
                                                                       mod = "sqrt_inv_esz.4",
                                                                       group = "StudyID",
                                                                       weights = "prop")

# FIGURE WITH A CONTINOUS VARIABLE:
figure_meta.regression.small.study.effects.n4 <- orchaRd::bubble_plot(meta.regression.small.study.effects.n4,
                                                                      mod = "sqrt_inv_esz.4",
                                                                      group = "StudyID",
                                                                      xlab = "sqrt_inv_esz.4",
                                                                      legend.pos = "bottom.right")


# STATISTICAL ANALYSIS: testing the small-study effects per hormone
meta.regression.small.study.effects.by.hormone <- rma.mv(cor,
                                                         VCV_ESVar,
                                                         mods = ~ 1 + Hormone_measured_general * sqrt_inv_esz,
                                                         random = list(~ 1 | StudyID,
                                                                       ~ 1 | LaboratoryID,
                                                                       ~ 1 | PopulationID,
                                                                       ~ 1 | Species,
                                                                       ~ 1 | Species_phylo,
                                                                       ~ 1 | EffectID),
                                                         method = "REML",
                                                         R = list(Species_phylo = phylo_cor),
                                                         test = "t",
                                                         data = meta.final_ok_ok,
                                                         control=list(rel.tol=1e-8))

#save(meta.regression.small.study.effects.by.hormone, file = "data/outputs/statistical_models/meta_regression_small_study_effects_by_hormone.RData")
load(file = "data/outputs/statistical_models/meta_regression_small_study_effects_by_hormone.RData") #meta.regression.small.study.effects.by.hormone


# Printing the summary results of the model
print(meta.regression.small.study.effects.by.hormone, digits = 3)

# Calculate marginal R2 with r2_ml
R2.meta.regression.small.study.effects.by.hormone <- r2_ml(meta.regression.small.study.effects.by.hormone) 
round(R2.meta.regression.small.study.effects.by.hormone * 100, 1)


# # TABLE WITH RESULTS FOR A CONTINUOUS VARIABLE:
results_meta.regression.small.study.effects.by.hormone <- orchaRd::mod_results(meta.regression.small.study.effects.by.hormone,
                                                                               mod = "Hormone_measured_general", 
                                                                               group = "StudyID",
                                                                               weights = "prop",
                                                                               by = "sqrt_inv_esz")

round(as.data.frame(results_meta.regression.small.study.effects.by.hormone[[1]])[,c(2:7)], 3)

# FIGURE WITH A CONTINOUS VARIABLE:
figure_results_meta.regression.small.study.effects.by.hormone.1 <- orchaRd::bubble_plot(meta.regression.small.study.effects.by.hormone,
                                                                                        group = "StudyID",
                                                                                        mod = "sqrt_inv_esz",
                                                                                        xlab = "sqrt_inv_esz",
                                                                                        by = "Hormone_measured_general",
                                                                                        legend.pos = "bottom.right")

# Showing the effect per hormone too
figure_meta.regression.small.study.effects.by.hormone <- orchaRd::orchard_plot(results_meta.regression.small.study.effects.by.hormone,
                                                                               mod = "sqrt_inv_esz",
                                                                               group = "StudyID",
                                                                               xlab = "Effect size",
                                                                               trunk.size = 1,
                                                                               branch.size = 2,
                                                                               twig.size = 1)

# I calculate the number of studies and effect sizes per hormone type to be
# included in the manuscript.

# Androgens
unique(meta.final_ok_ok_androgens$EffectID)
unique(meta.final_ok_ok_androgens$StudyID)

# Glucocorticoids
meta.final_ok_ok_cort <- droplevels(subset(meta.final_ok_ok, meta.final_ok_ok$Hormone_measured_general == 'corticosterone'))
unique(meta.final_ok_ok_cort$EffectID)
unique(meta.final_ok_ok_cort$StudyID)


# Thyroids
meta.final_ok_ok_ths <- droplevels(subset(meta.final_ok_ok, meta.final_ok_ok$Hormone_measured_general == 'TH3 and TH4'))
unique(meta.final_ok_ok_ths$EffectID)
unique(meta.final_ok_ok_ths$StudyID)


# STATISTICAL ANALYSIS: testing the small-study effects for each hormone group
# in separate analysis. 
# We opt to do these tests because sample sizes differed between hormone groups
# and we would like to test for the potential asymmetry more independently.
# We fear that by considering all three groups together, group-specific patterns 
# might get masked. 

# ANDROGENS:
# I need to create a new phylogenetic tree for this data set.
# 
# # PHYLOGENETIC TREE:
# resolved_names_androgens <- tnrs_match_names(as.character(unique(meta.final_ok_ok_androgens$Species)))
# 
# # Saving the taxonomic data created on the 5th February 2025 to speed the
# # process in the future and allow full reproducibility
# save(resolved_names_androgens, file = "data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_androgens.RData")

# Loading the taxonomic data created on the 5th February 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_androgens.RData") #resolved_names_androgens

# # extracting phylogenetic information
# my_tree_androgens <- tol_induced_subtree(ott_ids =
#                                          resolved_names_androgens[,"ott_id"],
#                                          label_format = "name")
# # # Quick tree plotting
# # plot(my_tree_androgens, no.margin = TRUE)
# 
# # We need to check for the existence of polytomies
# is.binary(my_tree_androgens)
# # Yes, meaning there are no polytomies. Let's go on.
# 
# # To confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# my_tree_androgens$tip.label <- gsub("_", " ", my_tree_androgens$tip.label)
# 
# intersect(as.character(my_tree_androgens$tip.label),
#           as.character(meta.final_ok_ok_androgens$Species))
# 
# # Listed in our database but not in the tree
# setdiff(as.character(meta.final_ok_ok_androgens$Species),
#         as.character(my_tree_androgens$tip.label))
# 
# # Listed in the tree but not in our database
# setdiff(as.character(my_tree_androgens$tip.label),
#         as.character(meta.final_ok_ok_androgens$Species))
# # No error or inconsistencies found.
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# # 5th February 2025
# save(my_tree_androgens, file = "data/outputs/phylogenetic_files/tree_androgens.Rdata")

# We can now load the saved tree
load("data/outputs/phylogenetic_files/tree_androgens.Rdata") #my_tree_androgens

# # Compute branch lengths of tree
# phylo_branch_androgens <- compute.brlen(my_tree_androgens, method = "Grafen", power = 1)
# 
# # Check if tree is ultrametric
# is.ultrametric(phylo_branch_androgens)
# # TRUE
# 
# # Matrix to be included in the models
# phylo_cor_androgens <- vcv(phylo_branch_androgens, cor = T)
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(phylo_cor_androgens, file = "data/outputs/phylogenetic_files/phylo_cor_androgens.Rdata")

# we can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor_androgens.Rdata") #phylo_cor_androgens

# Creating a duplicate species variable for the phylogenetic analysis
meta.final_ok_ok_androgens$Species_phylo_androgens <- meta.final_ok_ok_androgens$Species


# VARIANCE-COVARIANCE MATRIX:
VCV_ESVar_androgens <- matrix(0, nrow = nrow(meta.final_ok_ok_androgens), 
                              ncol = nrow(meta.final_ok_ok_androgens))

# Names rows and columns for each obsID
rownames(VCV_ESVar_androgens) <- meta.final_ok_ok_androgens[, "EffectID"]
colnames(VCV_ESVar_androgens) <- meta.final_ok_ok_androgens[, "EffectID"]

# Finds effect sizes that come from the same study
shared_coord_androgens <- which(meta.final_ok_ok_androgens[, "StudyID"] %in% 
                                  meta.final_ok_ok_androgens[duplicated(meta.final_ok_ok_androgens[, "StudyID"]), 
                                                             "StudyID"] == TRUE)

combinations_androgens <- do.call("rbind", tapply(shared_coord_androgens, 
                                                  meta.final_ok_ok_androgens[shared_coord_androgens, "StudyID"], 
                                                  function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations_androgens)[1]) {
  p1_androgens <- combinations_androgens[i, 1]
  p2_androgens <- combinations_androgens[i, 2]
  p1_p2_cov_androgens <- 0.5 * sqrt(meta.final_ok_ok_androgens[p1_androgens, "cor_var"]) * 
    sqrt(meta.final_ok_ok_androgens[p2_androgens, "cor_var"])
  VCV_ESVar_androgens[p1_androgens, p2_androgens] <- p1_p2_cov_androgens
  VCV_ESVar_androgens[p2_androgens, p1_androgens] <- p1_p2_cov_androgens
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(VCV_ESVar_androgens) <- meta.final_ok_ok_androgens[, "cor_var"]

# # In case you want to visually double check the matrix outside of R
# write.csv(VCV_ESVar_androgens, 'data/outputs/variance-covariance_matrices/VCV_ESVar_androgens.csv')


# I can now create the new variables to test for small-study effects.
meta.final_ok_ok_androgens$inv_esz <-  1/meta.final_ok_ok_androgens$final_n
meta.final_ok_ok_androgens$inv_esz.4 <-  1/(meta.final_ok_ok_androgens$final_n/4)
meta.final_ok_ok_androgens$sqrt_inv_esz  <-  with(meta.final_ok_ok_androgens, sqrt(inv_esz))
meta.final_ok_ok_androgens$sqrt_inv_esz.4 <-  with(meta.final_ok_ok_androgens, sqrt(inv_esz.4))


# STATISTICAL ANALYSIS:
meta.regression.small.study.effects_androgens <- rma.mv(cor,
                                                        VCV_ESVar_androgens,
                                                        mods = ~ 1 + sqrt_inv_esz,
                                                        random = list(~ 1 | StudyID,
                                                                      ~ 1 | LaboratoryID,
                                                                      ~ 1 | PopulationID,
                                                                      ~ 1 | Species,
                                                                      ~ 1 | Species_phylo_androgens,
                                                                      ~ 1 | EffectID),
                                                        method = "REML",
                                                        R = list(Species_phylo_androgens = phylo_cor_androgens),
                                                        test = "t",
                                                        data = meta.final_ok_ok_androgens)

#save(meta.regression.small.study.effects_androgens, file = "data/outputs/statistical_models/meta_regression_small_study_effects_androgens.RData")
load(file = "data/outputs/statistical_models/meta_regression_small_study_effects_androgens.RData") #meta.regression.small.study.effects

# Printing the summary results of the model
print(meta.regression.small.study.effects_androgens, digits = 3)
#predict(meta.regression.small.study.effects_androgens, digits = 3)

# Calculate marginal R2 with r2_ml
R2.meta.regression.small.study.effects_androgens <- r2_ml(meta.regression.small.study.effects_androgens) 
round(R2.meta.regression.small.study.effects_androgens * 100, 1)


# TABLE WITH RESULTS FOR A CONTINUOUS VARIABLE:
results_meta.regression.small.study.effects_androgens <- orchaRd::mod_results(meta.regression.small.study.effects_androgens,
                                                                              mod = "sqrt_inv_esz", 
                                                                              group = "StudyID",
                                                                              weights = "prop")

# FIGURE WITH A CONTINOUS VARIABLE:
figure_meta.regression.small.study.effects_androgens <- orchaRd::bubble_plot(meta.regression.small.study.effects_androgens, 
                                                                             mod = "sqrt_inv_esz", 
                                                                             group = "StudyID",
                                                                             xlab = "sqrt_inv_esz",
                                                                             legend.pos = "bottom.right")


# Glucocorticoids
# I need to create a new phylogenetic tree for this data set.

# # PHYLOGENETIC TREE:
# resolved_names_cort <- tnrs_match_names(as.character(unique(meta.final_ok_ok_cort$Species)))
# 
# # Saving the taxonomic data created on the 5th February 2025 to speed the
# # process in the future and allow full reproducibility
# save(resolved_names_cort, file = "data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_cort.RData")

# Loading the taxonomic data created on the 5th February 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_cort.RData") #resolved_names_cort

# # extracting phylogenetic information
# my_tree_cort <- tol_induced_subtree(ott_ids =
#                                          resolved_names_cort[,"ott_id"],
#                                          label_format = "name")
# # # Quick tree plotting
# # plot(my_tree_cort, no.margin = TRUE)
# 
# # We need to check for the existence of polytomies
# is.binary(my_tree_cort)
# # Yes, meaning there are no polytomies. Let's go on.
# 
# # To confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# my_tree_cort$tip.label <- gsub("_", " ", my_tree_cort$tip.label)
# 
# intersect(as.character(my_tree_cort$tip.label),
#           as.character(meta.final_ok_ok_cort$Species))
# 
# # Listed in our database but not in the tree
# setdiff(as.character(meta.final_ok_ok_cort$Species),
#         as.character(my_tree_cort$tip.label))
# 
# # Listed in the tree but not in our database
# setdiff(as.character(my_tree_cort$tip.label),
#         as.character(meta.final_ok_ok_cort$Species))
# # No error or inconsistencies found.
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# # 5th February 2025
# save(my_tree_cort, file = "data/outputs/phylogenetic_files/tree_cort.Rdata")

# We can now load the saved tree
load("data/outputs/phylogenetic_files/tree_cort.Rdata") #my_tree_cort

# # Compute branch lengths of tree
# phylo_branch_cort <- compute.brlen(my_tree_cort, method = "Grafen", power = 1)
# 
# # Check if tree is ultrametric
# is.ultrametric(phylo_branch_cort)
# # TRUE
# 
# # Matrix to be included in the models
# phylo_cor_cort <- vcv(phylo_branch_cort, cor = T)
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(phylo_cor_cort, file = "data/outputs/phylogenetic_files/phylo_cor_cort.Rdata")

# we can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor_cort.Rdata") #phylo_cor_cort

# Creating a duplicate species variable for the phylogenetic analysis
meta.final_ok_ok_cort$Species_phylo_cort <- meta.final_ok_ok_cort$Species


# VARIANCE-COVARIANCE MATRIX:
VCV_ESVar_cort <- matrix(0, nrow = nrow(meta.final_ok_ok_cort), 
                         ncol = nrow(meta.final_ok_ok_cort))

# Names rows and columns for each obsID
rownames(VCV_ESVar_cort) <- meta.final_ok_ok_cort[, "EffectID"]
colnames(VCV_ESVar_cort) <- meta.final_ok_ok_cort[, "EffectID"]

# Finds effect sizes that come from the same study
shared_coord_cort <- which(meta.final_ok_ok_cort[, "StudyID"] %in% 
                             meta.final_ok_ok_cort[duplicated(meta.final_ok_ok_cort[, "StudyID"]), 
                                                   "StudyID"] == TRUE)

combinations_cort <- do.call("rbind", tapply(shared_coord_cort, 
                                             meta.final_ok_ok_cort[shared_coord_cort, "StudyID"], 
                                             function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations_cort)[1]) {
  p1_cort <- combinations_cort[i, 1]
  p2_cort <- combinations_cort[i, 2]
  p1_p2_cov_cort <- 0.5 * sqrt(meta.final_ok_ok_cort[p1_cort, "cor_var"]) * 
    sqrt(meta.final_ok_ok_cort[p2_cort, "cor_var"])
  VCV_ESVar_cort[p1_cort, p2_cort] <- p1_p2_cov_cort
  VCV_ESVar_cort[p2_cort, p1_cort] <- p1_p2_cov_cort
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(VCV_ESVar_cort) <- meta.final_ok_ok_cort[, "cor_var"]

# # In case you want to visually double check the matrix outside of R
# write.csv(VCV_ESVar_cort, 'data/outputs/variance-covariance_matrices/VCV_ESVar_cort.csv')

# I can now create the new variables for calculating the small-study effect
# size
meta.final_ok_ok_cort$inv_esz <-  1/meta.final_ok_ok_cort$final_n
meta.final_ok_ok_cort$inv_esz.4 <-  1/(meta.final_ok_ok_cort$final_n/4)
meta.final_ok_ok_cort$sqrt_inv_esz  <-  with(meta.final_ok_ok_cort, sqrt(inv_esz))
meta.final_ok_ok_cort$sqrt_inv_esz.4 <-  with(meta.final_ok_ok_cort, sqrt(inv_esz.4))


# STATISTICAL ANALYSIS:
meta.regression.small.study.effects_cort <- rma.mv(cor,
                                                   VCV_ESVar_cort,
                                                   mods = ~ 1 + sqrt_inv_esz,
                                                   random = list(~ 1 | StudyID,
                                                                 ~ 1 | LaboratoryID,
                                                                 ~ 1 | PopulationID,
                                                                 ~ 1 | Species,
                                                                 ~ 1 | Species_phylo_cort,
                                                                 ~ 1 | EffectID),
                                                   method = "REML",
                                                   R = list(Species_phylo_cort = phylo_cor_cort),
                                                   test = "t",
                                                   data = meta.final_ok_ok_cort)

#save(meta.regression.small.study.effects_cort, file = "data/outputs/statistical_models/meta_regression_small_study_effects_cort.RData")
load(file = "data/outputs/statistical_models/meta_regression_small_study_effects_cort.RData") #meta.regression.small.study.effects

# Printing the summary results of the model
print(meta.regression.small.study.effects_cort, digits = 3)
#predict(meta.regression.small.study.effects_cort, digits = 3)

# Calculate marginal R2 with r2_ml
R2.meta.regression.small.study.effects_cort <- r2_ml(meta.regression.small.study.effects_cort) 
round(R2.meta.regression.small.study.effects_cort * 100, 1)


# TABLE WITH RESULTS FOR A CONTINUOUS VARIABLE:
results_meta.regression.small.study.effects_cort <- orchaRd::mod_results(meta.regression.small.study.effects_cort,
                                                                         mod = "sqrt_inv_esz", 
                                                                         group = "StudyID",
                                                                         weights = "prop")

# FIGURE WITH A CONTINOUS VARIABLE:
figure_meta.regression.small.study.effects_cort <- orchaRd::bubble_plot(meta.regression.small.study.effects_cort, 
                                                                        mod = "sqrt_inv_esz", 
                                                                        group = "StudyID",
                                                                        xlab = "sqrt_inv_esz",
                                                                        legend.pos = "bottom.right")


# 
# Thyroid hormones
# I need to create a new phylogenetic tree for this data set.

# # PHYLOGENETIC TREE:
# resolved_names_ths <- tnrs_match_names(as.character(unique(meta.final_ok_ok_ths$Species)))
# 
# # Saving the taxonomic data created on the 5th February 2025 to speed the
# # process in the future and allow full reproducibility
# save(resolved_names_ths, file = "data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_ths.RData")

# Loading the taxonomic data created on the 5th February 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_ths.RData") #resolved_names_ths

# # extracting phylogenetic information
# my_tree_ths <- tol_induced_subtree(ott_ids =
#                                          resolved_names_ths[,"ott_id"],
#                                          label_format = "name")
# # # Quick tree plotting
# # plot(my_tree_ths, no.margin = TRUE)
# 
# # We need to check for the existence of polytomies
# is.binary(my_tree_ths)
# # Yes, meaning there are no polytomies. Let's go on.
# 
# # To confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# my_tree_ths$tip.label <- gsub("_", " ", my_tree_ths$tip.label)
# 
# intersect(as.character(my_tree_ths$tip.label),
#           as.character(meta.final_ok_ok_ths$Species))
# 
# # Listed in our database but not in the tree
# setdiff(as.character(meta.final_ok_ok_ths$Species),
#         as.character(my_tree_ths$tip.label))
# 
# # Listed in the tree but not in our database
# setdiff(as.character(my_tree_ths$tip.label),
#         as.character(meta.final_ok_ok_ths$Species))
# # No error or inconsistencies found.
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# # 5th February 2025
# save(my_tree_ths, file = "data/outputs/phylogenetic_files/tree_ths.Rdata")

# We can now load the saved tree
load("data/outputs/phylogenetic_files/tree_ths.Rdata") #my_tree_ths

# # Compute branch lengths of tree
# phylo_branch_ths <- compute.brlen(my_tree_ths, method = "Grafen", power = 1)
# 
# # Check if tree is ultrametric
# is.ultrametric(phylo_branch_ths)
# # TRUE
# 
# # Matrix to be included in the models
# phylo_cor_ths <- vcv(phylo_branch_ths, cor = T)
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(phylo_cor_ths, file = "data/outputs/phylogenetic_files/phylo_cor_ths.Rdata")

# we can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor_ths.Rdata") #phylo_cor_ths

# Creating a duplicate species variable for the phylogenetic analysis
meta.final_ok_ok_ths$Species_phylo_ths <- meta.final_ok_ok_ths$Species


# VARIANCE-COVARIANCE MATRIX:
VCV_ESVar_ths <- matrix(0, nrow = nrow(meta.final_ok_ok_ths), 
                        ncol = nrow(meta.final_ok_ok_ths))

# Names rows and columns for each obsID
rownames(VCV_ESVar_ths) <- meta.final_ok_ok_ths[, "EffectID"]
colnames(VCV_ESVar_ths) <- meta.final_ok_ok_ths[, "EffectID"]

# Finds effect sizes that come from the same study
shared_coord_ths <- which(meta.final_ok_ok_ths[, "StudyID"] %in% 
                            meta.final_ok_ok_ths[duplicated(meta.final_ok_ok_ths[, "StudyID"]), 
                                                 "StudyID"] == TRUE)

combinations_ths <- do.call("rbind", tapply(shared_coord_ths, 
                                            meta.final_ok_ok_ths[shared_coord_ths, "StudyID"], 
                                            function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations_ths)[1]) {
  p1_ths <- combinations_ths[i, 1]
  p2_ths <- combinations_ths[i, 2]
  p1_p2_cov_ths <- 0.5 * sqrt(meta.final_ok_ok_ths[p1_ths, "cor_var"]) * 
    sqrt(meta.final_ok_ok_ths[p2_ths, "cor_var"])
  VCV_ESVar_ths[p1_ths, p2_ths] <- p1_p2_cov_ths
  VCV_ESVar_ths[p2_ths, p1_ths] <- p1_p2_cov_ths
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(VCV_ESVar_ths) <- meta.final_ok_ok_ths[, "cor_var"]

# # In case you want to visually double check the matrix outside of R
#write.csv(VCV_ESVar_ths, 'data/outputs/variance-covariance_matrices/VCV_ESVar_ths.csv')

# I can now create the new variables for calculating the small-study effect
# size
meta.final_ok_ok_ths$inv_esz <-  1/meta.final_ok_ok_ths$final_n
meta.final_ok_ok_ths$inv_esz.4 <-  1/(meta.final_ok_ok_ths$final_n/4)
meta.final_ok_ok_ths$sqrt_inv_esz  <-  with(meta.final_ok_ok_ths, sqrt(inv_esz))
meta.final_ok_ok_ths$sqrt_inv_esz.4 <-  with(meta.final_ok_ok_ths, sqrt(inv_esz.4))


# STATISTICAL ANALYSIS:
meta.regression.small.study.effects_ths <- rma.mv(cor,
                                                  VCV_ESVar_ths,
                                                  mods = ~ 1 + sqrt_inv_esz,
                                                  random = list(~ 1 | StudyID,
                                                                ~ 1 | LaboratoryID,
                                                                ~ 1 | PopulationID,
                                                                ~ 1 | Species,
                                                                ~ 1 | Species_phylo_ths,
                                                                ~ 1 | EffectID),
                                                  method = "REML",
                                                  R = list(Species_phylo_ths = phylo_cor_ths),
                                                  test = "t",
                                                  data = meta.final_ok_ok_ths)

#save(meta.regression.small.study.effects_ths, file = "data/outputs/statistical_models/meta_regression_small_study_effects_ths.RData")
load(file = "data/outputs/statistical_models/meta_regression_small_study_effects_ths.RData") #meta.regression.small.study.effects

# Printing the summary results of the model
print(meta.regression.small.study.effects_ths, digits = 3)
#predict(meta.regression.small.study.effects_ths, digits = 3)

# Calculate marginal R2 with r2_ml
R2.meta.regression.small.study.effects_ths <- r2_ml(meta.regression.small.study.effects_ths) 
round(R2.meta.regression.small.study.effects_ths * 100, 1)


# TABLE WITH RESULTS FOR A CONTINUOUS VARIABLE:
results_meta.regression.small.study.effects_ths <- orchaRd::mod_results(meta.regression.small.study.effects_ths,
                                                                        mod = "sqrt_inv_esz", 
                                                                        group = "StudyID",
                                                                        weights = "prop")

# FIGURE WITH A CONTINOUS VARIABLE:
figure_meta.regression.small.study.effects_ths <- orchaRd::bubble_plot(meta.regression.small.study.effects_ths, 
                                                                       mod = "sqrt_inv_esz", 
                                                                       group = "StudyID",
                                                                       xlab = "sqrt_inv_esz",
                                                                       legend.pos = "bottom.right")


################################################################################
# Time-lag bias or decline effects
################################################################################

# Mean-centring year before including it in the model
meta.final_ok_ok$Year.c <- as.vector(scale(as.numeric(as.character(meta.final_ok_ok$Year)), 
                                           scale = F))

# STATISTICAL ANALYSIS:
meta.regression.decline.effects <- rma.mv(cor,
                                          VCV_ESVar,
                                          mods = ~ 1 + Year.c,
                                          random = list(~ 1 | StudyID,
                                                        ~ 1 | LaboratoryID,
                                                        ~ 1 | PopulationID,
                                                        ~ 1 | Species,
                                                        ~ 1 | Species_phylo,
                                                        ~ 1 | EffectID),
                                          method = "REML",
                                          R = list(Species_phylo = phylo_cor),
                                          test = "t",
                                          data = meta.final_ok_ok)

#save(meta.regression.decline.effects, file = "data/outputs/statistical_models/meta_regression_decline_effects.RData")
load(file = "data/outputs/statistical_models/meta_regression_decline_effects.RData") #meta.regression.decline.effects

# Printing the summary results of the model
print(meta.regression.decline.effects, digits = 3)

# Calculate marginal R2 with r2_ml
R2.meta.regression.decline.effects <- r2_ml(meta.regression.decline.effects) 
round(R2.meta.regression.decline.effects * 100, 1)


# TABLE WITH RESULTS FOR A CONTINUOUS VARIABLE:
results_meta.regression.decline.effects <- orchaRd::mod_results(meta.regression.decline.effects,
                                                                mod = "Year.c", 
                                                                group = "StudyID",
                                                                weights = "prop")

# FIGURE WITH A CONTINOUS VARIABLE:
figure_meta.regression.decline.effects <- orchaRd::bubble_plot(meta.regression.decline.effects, 
                                                               mod = "Year.c", 
                                                               group = "StudyID",
                                                               xlab = "Year.c",
                                                               legend.pos = "bottom.left")


# STATISTICAL ANALYSIS: testing the decline effect per hormone
meta.regression.decline.effects.by.hormone <- rma.mv(cor,
                                                     VCV_ESVar,
                                                     mods = ~ 1 + Hormone_measured_general * Year.c,
                                                     random = list(~ 1 | StudyID,
                                                                   ~ 1 | LaboratoryID,
                                                                   ~ 1 | PopulationID,
                                                                   ~ 1 | Species,
                                                                   ~ 1 | Species_phylo,
                                                                   ~ 1 | EffectID),
                                                     method = "REML",
                                                     R = list(Species_phylo = phylo_cor),
                                                     test = "t",
                                                     data = meta.final_ok_ok)

#save(meta.regression.decline.effects.by.hormone, file = "data/outputs/statistical_models/meta_regression_decline_effects_by_hormone.RData")
load(file = "data/outputs/statistical_models/meta_regression_decline_effects_by_hormone.RData") #meta.regression.decline.effects.by.hormone

# Printing the summary results of the model
print(meta.regression.decline.effects.by.hormone, digits = 3)

# Calculate marginal R2 with r2_ml
R2.meta.regression.decline.effects.by.hormone <- r2_ml(meta.regression.decline.effects.by.hormone) 
round(R2.meta.regression.decline.effects.by.hormone * 100, 1)


# # TABLE WITH RESULTS FOR A CONTINUOUS VARIABLE:
results_meta.regression.decline.effects.by.hormone <- orchaRd::mod_results(meta.regression.decline.effects.by.hormone,
                                                                           mod = "Hormone_measured_general", 
                                                                           group = "StudyID",
                                                                           weights = "prop",
                                                                           by = "Year.c")


round(as.data.frame(results_meta.regression.decline.effects.by.hormone[[1]])[,c(2:7)], 3)

# FIGURE WITH A CONTINOUS VARIABLE:
figure_results_meta.regression.decline.effects.by.hormone.1 <- orchaRd::bubble_plot(meta.regression.decline.effects.by.hormone,
                                                                                    group = "StudyID",
                                                                                    mod = "Year.c",
                                                                                    xlab = "Year.c",
                                                                                    by = "Hormone_measured_general",
                                                                                    legend.pos = "bottom.right")

# # Showing the effect per hormone too
figure_meta.regression.decline.effects.by.hormone <- orchaRd::orchard_plot(results_meta.regression.decline.effects.by.hormone,
                                                                           mod = "Year.c",
                                                                           group = "StudyID",
                                                                           xlab = "Effect size",
                                                                           trunk.size = 2,
                                                                           branch.size = 2,
                                                                           twig.size = 1)


################################################################################
# Reporting effects
################################################################################

# In our pre-registration we mentioned that we wanted to test for reporting.
# This is, for how 1) data completeness and 2) selective reporting might influence 
# the results obtained. We also pre-registered to test for 3) blinding effects

# 1) DATA COMPLETENESS
table(meta.final_ok_ok$DataReporting)

table(meta.final_ok_ok$DataReporting, meta.final_ok_ok$Hormone_measured_general)
# In total, we have 407 effect sizes for which we have complete information, 
# and 36 for which we don't. It is important to keep in mind that the number of 
# incomplete studies is underestimated because these are the studies for which, 
# even if the data was still incomplete, we were able to obtain the missing data 
# to perform the analysis (36 effect sizes). Unfortunately, for the majority of 
# the incomplete studies (230 effect sizes) were not able to obtain this data 
# and therefore they were not included in this meta-analysis. 

# Besides this, and following our pre-registration, we do not have enough data to 
# test this (> 5 points per category). Hence, I will only test this for androgens. 

table(meta.final_ok_ok_androgens$DataReporting)


# We have 3 variables, but this is because there is a typo in one 'complete'
# I will modify this: 
meta.final_ok_ok_androgens$DataReporting <- str_replace(meta.final_ok_ok_androgens$DataReporting,
                                                        "complete ",
                                                        "complete")
table(meta.final_ok_ok_androgens$DataReporting)
# For androgens, we have 304 complete rows and 34 incomplete rows.


meta.regression.completeness <- rma.mv(cor,
                                       VCV_ESVar_androgens,
                                       mods = ~ 1 + DataReporting,
                                       random = list(~ 1 | StudyID,
                                                     ~ 1 | LaboratoryID,
                                                     ~ 1 | PopulationID,
                                                     ~ 1 | Species,
                                                     ~ 1 | Species_phylo_androgens,
                                                     ~ 1 | EffectID),
                                       method = "REML",
                                       #verbose = TRUE,
                                       control = list(rel.tol = 1e-8),
                                       R = list(Species_phylo_androgens = phylo_cor_androgens),
                                       test = "t",
                                       data = meta.final_ok_ok_androgens)
#
# We had a Convergence Problem with the rma.mv() Function. 
# By using the command "verbose = TRUE" we obtain the information that starting 
# around the iteration 297 there is no further change in the log likelihood when 
# rounded to 4 decimal placed. 
# However, because the default settings of the model had a very small threshold, 
# the model could not for determined when convergence occurred.
# We therefore modified this threshold with the commands "control = list(...)".
# Results remain the same and we no longer have a error message.

#save(meta.regression.completeness, file = "data/outputs/statistical_models/meta.regression.completeness.RData")
load(file = "data/outputs/statistical_models/meta.regression.completeness.RData") #meta.regression.completeness

# Printing the summary results of the model 
print(meta.regression.completeness, digits = 3)

# Calculate marginal R2 with r2_ml
R2.meta.regression.completeness <- r2_ml(meta.regression.completeness) 
round(R2.meta.regression.completeness * 100, 1)


# TABLE:
results_meta.regression.completeness <- orchaRd::mod_results(meta.regression.completeness, 
                                                             mod = "DataReporting", 
                                                             group = "StudyID", 
                                                             subset = TRUE) 


round(as.data.frame(results_meta.regression.completeness[[1]])[,c(2:6)], 3)


# FIGURE: effect of each group of hormones on maternal fitness traits.
fig_meta.regression.completeness <- orchaRd::orchard_plot(meta.regression.completeness, 
                                                          mod = "DataReporting", 
                                                          group = "StudyID", 
                                                          xlab = "Effect size",
                                                          trunk.size = 2,
                                                          branch.size = 2,
                                                          twig.size = 1)

# Complete and incomplete studies provide similar results. Complete studies 
# tend to show a stronger negative effect than incomplete studies.



## 2) SELECTIVE REPORTING
meta.final_ok_ok$Partial_or_selective_data_rep <- as.factor(meta.final_ok_ok$Partial_or_selective_data_rep)
table(meta.final_ok_ok$Partial_or_selective_data_rep)

# I need to rename 2 variables:
meta.final_ok_ok$Partial_or_selective_data_rep <- recode(meta.final_ok_ok$Partial_or_selective_data_rep, 
                                                         "yes (backwards selection)" = "yes")
meta.final_ok_ok$Partial_or_selective_data_rep <- recode(meta.final_ok_ok$Partial_or_selective_data_rep, 
                                                         "yes (some eggs excluded because of bacterial infections)" = "yes")

table(meta.final_ok_ok$Partial_or_selective_data_rep, meta.final_ok_ok$Hormone_measured_general)
# We have enough data for testing this in the three hormone groups.

# It is important to keep in mind that for the majority of the effect sizes
# extracted, we could not confidently decide whether researchers did or did not
# selectively report data (N = 171/443).


meta.regression.selective_data_rep <- rma.mv(cor,
                                             VCV_ESVar,
                                             mods = ~ Hormone_measured_general * Partial_or_selective_data_rep,
                                             random = list(~ 1 | StudyID,
                                                           ~ 1 | LaboratoryID,
                                                           ~ 1 | PopulationID,
                                                           ~ 1 | Species,
                                                           ~ 1 | Species_phylo,
                                                           ~ 1 | EffectID),
                                             method = "REML",
                                             R = list(Species_phylo = phylo_cor),
                                             test = "t",
                                             data = meta.final_ok_ok)

#save(meta.regression.selective_data_rep, file = "data/outputs/statistical_models/meta.regression.selective_data_rep.RData")
load(file = "data/outputs/statistical_models/meta.regression.selective_data_rep.RData") #meta.regression.selective_data_rep

# Printing the summary results of the model
print(meta.regression.selective_data_rep, digits = 3)

# Calculate marginal R2 with r2_ml
R2.method_selective_data_rep <- r2_ml(meta.regression.selective_data_rep)
round(R2.method_selective_data_rep * 100, 1)


# In order to better understand and plot the results, I will create an artificial
# variable that includes this interaction.
# But first, I need to exclude those rows that have 'NAs'. Else, I will have
# three categories per hormone (yes, no, NA)

meta.final_ok_ok_sr <- droplevels(subset(meta.final_ok_ok, 
                                         meta.final_ok_ok$Partial_or_selective_data_rep != "NA"))

nrow(meta.final_ok_ok_sr)

unique(meta.final_ok_ok_sr$StudyID)


# # I need to create a new phylogenetic tree because this subset has only 12 species.
# # PHYLOGENETIC TREE:
# resolved_names_selective_rep <- tnrs_match_names(as.character(unique(meta.final_ok_ok_sr$Species)))
# 
# # Saving the taxonomic data created on the 5th February 2025 to speed the
# # process in the future and allow full reproducibility
# save(resolved_names_selective_rep, file = "data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_selective_rep.RData")

# Loading the taxonomic data created on the 5th February 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_selective_rep.RData") #resolved_names_selective_rep

# # extracting phylogenetic information
# my_tree_selective_rep <- tol_induced_subtree(ott_ids =
#                                          resolved_names_selective_rep[,"ott_id"],
#                                          label_format = "name")
# # # Quick tree plotting
# # plot(my_tree_selective_rep, no.margin = TRUE)
# 
# # We need to check for the existence of polytomies
# is.binary(my_tree_selective_rep)
# # Yes, meaning there are no polytomies. Let's go on.
# 
# # To confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# my_tree_selective_rep$tip.label <- gsub("_", " ", my_tree_selective_rep$tip.label)
# 
# intersect(as.character(my_tree_selective_rep$tip.label),
#           as.character(meta.final_ok_ok_sr$Species))
# 
# # Listed in our database but not in the tree
# setdiff(as.character(meta.final_ok_ok_sr$Species),
#         as.character(my_tree_selective_rep$tip.label))
# 
# # Listed in the tree but not in our database
# setdiff(as.character(my_tree_selective_rep$tip.label),
#         as.character(meta.final_ok_ok_sr$Species))
# # No error or inconsistencies found.
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# # 5th February 2025
# save(my_tree_selective_rep, file = "data/outputs/phylogenetic_files/tree_selective_rep.Rdata")

# We can now load the saved tree
load("data/outputs/phylogenetic_files/tree_selective_rep.Rdata") #my_tree_selective_rep

# # Compute branch lengths of tree
# phylo_branch_selective_rep <- compute.brlen(my_tree_selective_rep, method = "Grafen", power = 1)
# 
# # Check if tree is ultrametric
# is.ultrametric(phylo_branch_selective_rep)
# # TRUE
# 
# # Matrix to be included in the models
# phylo_cor_selective_rep <- vcv(phylo_branch_selective_rep, cor = T)
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(phylo_cor_selective_rep, file = "data/outputs/phylogenetic_files/phylo_cor_selective_rep.Rdata")

# we can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor_selective_rep.Rdata") #phylo_cor_selective_rep

# Creating a duplicate species variable for the phylogenetic analysis
meta.final_ok_ok_sr$Species_phylo_selective_rep <- meta.final_ok_ok_sr$Species

# VARIANCE-COVARIANCE MATRIX:
VCV_ESVar_selective_rep <- matrix(0, nrow = nrow(meta.final_ok_ok_sr), 
                                  ncol = nrow(meta.final_ok_ok_sr))

# Names rows and columns for each obsID
rownames(VCV_ESVar_selective_rep) <- meta.final_ok_ok_sr[, "EffectID"]
colnames(VCV_ESVar_selective_rep) <- meta.final_ok_ok_sr[, "EffectID"]

# Finds effect sizes that come from the same study
shared_coord_selective_rep <- which(meta.final_ok_ok_sr[, "StudyID"] %in% 
                                      meta.final_ok_ok_sr[duplicated(meta.final_ok_ok_sr[, "StudyID"]), 
                                                          "StudyID"] == TRUE)

combinations_selective_rep <- do.call("rbind", tapply(shared_coord_selective_rep, 
                                                      meta.final_ok_ok_sr[shared_coord_selective_rep, "StudyID"], 
                                                      function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations_selective_rep)[1]) {
  p1_selective_rep <- combinations_selective_rep[i, 1]
  p2_selective_rep <- combinations_selective_rep[i, 2]
  p1_p2_cov_selective_rep <- 0.5 * sqrt(meta.final_ok_ok_sr[p1_selective_rep, "cor_var"]) * 
    sqrt(meta.final_ok_ok_sr[p2_selective_rep, "cor_var"])
  VCV_ESVar_selective_rep[p1_selective_rep, p2_selective_rep] <- p1_p2_cov_selective_rep
  VCV_ESVar_selective_rep[p2_selective_rep, p1_selective_rep] <- p1_p2_cov_selective_rep
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(VCV_ESVar_selective_rep) <- meta.final_ok_ok_sr[, "cor_var"]

# # In case you want to visually double check the matrix outside of R
# write.csv(VCV_ESVar_selective_rep, 'data/outputs/variance-covariance_matrices/VCV_ESVar_selective_rep.csv')


# Now I create an artificial variable that combines the hormone measured and information on
# whether the authors selectively reported or not data.
meta.final_ok_ok_sr$hormone_selective_rep <- paste(meta.final_ok_ok_sr$Hormone_measured_general,
                                                   meta.final_ok_ok_sr$Partial_or_selective_data_rep,
                                                   sep = "_")

# we will now run a new model with this artificial variable.
meta.regression.selective_data_rep_artificial <- rma.mv(cor,
                                                        VCV_ESVar_selective_rep,
                                                        mods = ~ hormone_selective_rep - 1,
                                                        random = list(~ 1 | StudyID,
                                                                      ~ 1 | LaboratoryID,
                                                                      ~ 1 | PopulationID,
                                                                      ~ 1 | Species,
                                                                      ~ 1 | Species_phylo_selective_rep,
                                                                      ~ 1 | EffectID),
                                                        method = "REML",
                                                        R = list(Species_phylo_selective_rep = phylo_cor_selective_rep),
                                                        test = "t",
                                                        data = meta.final_ok_ok_sr)


#save(meta.regression.selective_data_rep_artificial, file = "data/outputs/statistical_models/meta.regression.selective_data_rep_artificial.RData")
load(file = "data/outputs/statistical_models/meta.regression.selective_data_rep_artificial.RData") #meta.regression.selective_data_rep_artificial

# Printing the summary results of the model
print(meta.regression.selective_data_rep_artificial, digits = 3)

# Calculate marginal R2 with r2_ml
R2.method_selective_data_rep_artificial <- r2_ml(meta.regression.selective_data_rep_artificial)
round(R2.method_selective_data_rep_artificial * 100, 1)

# TABLE WITH RESULTS:
results_selective_data_rep_artificial <- orchaRd::mod_results(meta.regression.selective_data_rep_artificial, 
                                                              mod = "hormone_selective_rep", 
                                                              group = "StudyID", 
                                                              subset = TRUE)

round(as.data.frame(results_selective_data_rep_artificial[[1]])[,c(2:6)], 3)


# FIGURE: 
fig_selective_data_rep_artificial <- orchaRd::orchard_plot(meta.regression.selective_data_rep_artificial, 
                                                           mod = "hormone_selective_rep", 
                                                           group = "StudyID", 
                                                           xlab = "Effect size",
                                                           trunk.size = 2,
                                                           branch.size = 2,
                                                           twig.size = 1)


# Pair-wise comparisons including post-hoc Wald tests
meta.regression.selective_data_rep_artificial_pc <- rma.mv(cor,
                                                           VCV_ESVar_selective_rep,
                                                           mods = ~ hormone_selective_rep,
                                                           random = list(~ 1 | StudyID,
                                                                         ~ 1 | LaboratoryID,
                                                                         ~ 1 | PopulationID,
                                                                         ~ 1 | Species,
                                                                         ~ 1 | Species_phylo_selective_rep,
                                                                         ~ 1 | EffectID),
                                                           method = "REML",
                                                           R = list(Species_phylo_selective_rep = phylo_cor_selective_rep),
                                                           test = "t",
                                                           data = meta.final_ok_ok_sr)

#save(meta.regression.meh1_pc, file = "data/outputs/statistical_models/meta.regression.selective_data_rep_artificial_pc.RData")
load(file = "data/outputs/statistical_models/meta.regression.selective_data_rep_artificial_pc.RData") #meta.regression.meh1_pc

print(meta.regression.selective_data_rep_artificial_pc, digits = 3)


# Information on hormone group comparison between studies that partially reported data
# and studies that did no:

# Glucocorticoids:
car::linearHypothesis(meta.regression.selective_data_rep_artificial_pc, rbind(c(0,0,1,-1,0,0))) 


# Thyroids:
car::linearHypothesis(meta.regression.selective_data_rep_artificial_pc, rbind(c(0,0,0,0,1,-1))) 

#######################

#### 3) BLINDNESS
table(meta.final_ok_ok$Blind_data_recording, meta.final_ok_ok$Hormone_measured_general)
# Unfortunately, we do not have enough data point per groups to perform this
# test.


################################################################################
# Pseudo-ALL-IN model
################################################################################

# Since we had to use different subsets for each of our questions mostly due to
# lack of data, this pseudo-ALL-IN model, rather than including all moderators 
# we tested, includes those for which we had data for the entire dataset

# mean-centring Study_type to help with interpretations
meta.final_ok_ok$Study_type.num.c <- as.vector(scale(ifelse(meta.final_ok_ok$Study_type == "experimental",
                                                            -1,
                                                            1), 
                                                     scale = F))


meta.regression.all.in <- rma.mv(cor,
                                 VCV_ESVar,
                                 mods = ~ 1 +
                                   #Study_type.num.c +
                                   Hormone_measured_general * sqrt_inv_esz +
                                   Hormone_measured_general * Year.c,
                                 random = list(~ 1 | StudyID,
                                               ~ 1 | LaboratoryID,
                                               ~ 1 | PopulationID,
                                               ~ 1 | Species,
                                               ~ 1 | Species_phylo,
                                               ~ 1 | EffectID),
                                 method = "REML",
                                 R = list(Species_phylo = phylo_cor),
                                 test = "t",
                                 data = meta.final_ok_ok)

# save(meta.regression.all.in, file = "data/outputs/statistical_models/meta_regression_all_in.RData")
load(file = "data/outputs/statistical_models/meta_regression_all_in.RData") #meta.regression.all.in


# Printing the summary results of the model
print(meta.regression.all.in, digits = 3)

# Calculate marginal R2 with r2_ml
R2.meta.regression.all.in <- r2_ml(meta.regression.all.in) 
round(R2.meta.regression.all.in * 100, 1)


figure_results_meta.regression.all.in <- orchaRd::bubble_plot(meta.regression.all.in,
                                                              group = "StudyID",
                                                              mod = "sqrt_inv_esz",
                                                              xlab = "sqrt_inv_esz",
                                                              by = "Hormone_measured_general",
                                                              legend.pos = "bottom.right")

# TABLE WITH RESULTS FOR A CONTINUOUS VARIABLE:
results_meta.regression.all.in <- orchaRd::mod_results(meta.regression.all.in,
                                                       mod = "Hormone_measured_general",
                                                       group = "StudyID",
                                                       weights = "prop",
                                                       by = "sqrt_inv_esz")


# Showing the effect per hormone too
figure_results_meta.regression.all.in.by.hormone <- orchaRd::orchard_plot(results_meta.regression.all.in,
                                                                           mod = "sqrt_inv_esz",
                                                                           group = "StudyID",
                                                                           xlab = "Effect size",
                                                                           trunk.size = 2,
                                                                           branch.size = 2,
                                                                           twig.size = 1)




################################################################################
# Exploratory analysis - non pre-registered - TO BETTER UNDERSTAND THE 
# HETEROGENEINITY FOUND IN THE RANDOM EFFECT 'WITHIN-STUDY'

# We found that 'within-study' explains a high amount of heterogeneity. We think
# this is because researchers generally measure several fitness traits and discrepancies
# between egg hormones and fitness proxies might arise (e.g., some traits are not
# influenced by egg hormones while others are)

# With the idea of testing if this is the case and if there are some traits that
# are particularly influenced by egg hormones, we ran this new analysis. For this
# we use the entire data set and those fitness proxies that have at least 
# 5 data poits (as stated in our pre-registration).

################################################################################

# I first create a new column that combines maternal and offspring fitness
# traits into one.
meta.final_ok_ok$Fitness_trait_all <- meta.final_ok_ok$Fitness_offs

meta.final_ok_ok <- meta.final_ok_ok %>% 
  mutate(Fitness_trait_all = coalesce(Fitness_trait_all, 
                                      Fitness_mother))

meta.final_ok_ok$Fitness_trait_all <- as.factor(meta.final_ok_ok$Fitness_trait_all)

# I now check the information per level:
levels(meta.final_ok_ok$Fitness_trait_all)

# I now group fitness proxies in a way in which all fitness proxies related
# to similar traits are evaluated together. For example, all proxies related
# to offspring growth.

# For this, I first create a new column
meta.final_ok_ok$Fitness_trait_all_group <- meta.final_ok_ok$Fitness_trait_all

meta.final_ok_ok$Fitness_trait_all_group <- as.character(meta.final_ok_ok$Fitness_trait_all_group)

# Growth
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "growth rate (mass)"] <- "growth"
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "growth (PC1)"] <- "growth"
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "growth (PC2)"] <- "growth"
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "growth (tarsus)"] <- "growth"
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "growth rate (mass)"] <- "growth"
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "growth rate"] <- "growth"
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "growth rate (ulna)"] <- "growth"
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "growth rate (flipper length)"] <- "growth"
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "mass gain"] <- "growth"
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "growth (mass)"] <- "growth"
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "growth rate (tarsus)"] <- "growth"


# Wing
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "flipper length"] <- "wing"

# Structural body size
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "structural body size (mass and tarsus)"] <- "structural body size"
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "structural body size (mass and bill length)"] <- "structural body size"

# Hatching success
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "hatching failure"] <- "hatching success"
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "hatching probability"] <- "hatching success"


# Hatching number
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "egg mortality"] <- "hatching number"


# Fledging success
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "fledgling probability"] <- "fledgling success"
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "pre-fledgling survival probability"] <- "fledgling success"
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "fledgling success (fledglings/hatchings)"] <- "fledgling success"

# Offspring survival 
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "offspring mortality"] <- "offspring survival"

# Offspring survival years
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "offspring survival years"] <- "offspring recruit"
meta.final_ok_ok$Fitness_trait_all_group[meta.final_ok_ok$Fitness_trait_all_group == 
                                           "recruitment"] <- "offspring recruit"

meta.final_ok_ok$Fitness_trait_all_group <- as.factor(meta.final_ok_ok$Fitness_trait_all_group)
levels(meta.final_ok_ok$Fitness_trait_all_group)
# 19 levels

# I now check if all levels have at least 5 effect sizes (minimum number of 
# effect sizes that we stated in our pre-registration that we need for running
# statistical models)
table(meta.final_ok_ok$Fitness_trait_all_group)
# Levels for which we do not have enough effect sizes: head-bill length, head length,
# maternal life time reproductive success, and maternal longevity.


# I retain those traits for which we have more than 5 effect sizes:

meta.final_ok_ok_withinstudyID <- droplevels(subset(meta.final_ok_ok, 
                                                    meta.final_ok_ok$Fitness_trait_all_group == "beak flank width" |
                                                      meta.final_ok_ok$Fitness_trait_all_group == "clutch size" |  
                                                      meta.final_ok_ok$Fitness_trait_all_group == "culmen" |
                                                      meta.final_ok_ok$Fitness_trait_all_group == "fledgling number" |
                                                      meta.final_ok_ok$Fitness_trait_all_group == "fledgling success" |
                                                      meta.final_ok_ok$Fitness_trait_all_group == "gape width" |
                                                      meta.final_ok_ok$Fitness_trait_all_group == "hatching number" |
                                                      meta.final_ok_ok$Fitness_trait_all_group == "hatching success" |
                                                      meta.final_ok_ok$Fitness_trait_all_group == "mass" |  
                                                      meta.final_ok_ok$Fitness_trait_all_group == "offspring recruit" |
                                                      meta.final_ok_ok$Fitness_trait_all_group == "offspring survival" |
                                                      meta.final_ok_ok$Fitness_trait_all_group == "structural body size"| 
                                                      meta.final_ok_ok$Fitness_trait_all_group == "tarsus" |
                                                      meta.final_ok_ok$Fitness_trait_all_group == "growth" |
                                                      meta.final_ok_ok$Fitness_trait_all_group == "wing"))

# I will reorder variables to better understand results
meta.final_ok_ok_withinstudyID$Fitness_trait_all_group <- factor(meta.final_ok_ok_withinstudyID$Fitness_trait_all_group, 
                                                                 levels = c("clutch size", 
                                                                            "hatching number",
                                                                            "hatching success",
                                                                            "beak flank width",
                                                                            "culmen",
                                                                            "gape width",
                                                                            "mass",
                                                                            "tarsus",
                                                                            "wing",
                                                                            "structural body size",
                                                                            "growth",
                                                                            "fledgling number",
                                                                            "fledgling success",
                                                                            "offspring survival",
                                                                            "offspring recruit"))

# # I now need to create a new phylogenetic tree because this subset has only 12 species.
# # PHYLOGENETIC TREE:
# resolved_names_withinstudyID <- tnrs_match_names(as.character(unique(meta.final_ok_ok_withinstudyID$Species)))
# 
# # Saving the taxonomic data created on the 5th February 2025 to speed the
# # process in the future and allow full reproducibility
# save(resolved_names_withinstudyID, file = "data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_withinstudyID.RData")

# Loading the taxonomic data created on the 5th February 2025
load("data/outputs/phylogenetic_files/taxa_Open_Tree_of_Life_withinstudyID.RData") #resolved_names_withinstudyID

# # extracting phylogenetic information
# my_tree_withinstudyID <- tol_induced_subtree(ott_ids =
#                                                resolved_names_withinstudyID[,"ott_id"],
#                                              label_format = "name")
# # # Quick tree plotting
# plot(my_tree_withinstudyID, no.margin = TRUE)
# 
# # We need to check for the existence of polytomies
# is.binary(my_tree_withinstudyID)
# # Yes, meaning there are no polytomies. Let's go on.
# 
# # To confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# my_tree_withinstudyID$tip.label <- gsub("_", " ", my_tree_withinstudyID$tip.label)
# 
# intersect(as.character(my_tree_withinstudyID$tip.label),
#           as.character(meta.final_ok_ok_withinstudyID$Species))
# 
# # Listed in our database but not in the tree
# setdiff(as.character(meta.final_ok_ok_withinstudyID$Species),
#         as.character(my_tree_withinstudyID$tip.label))
# 
# # Listed in the tree but not in our database
# setdiff(as.character(my_tree_withinstudyID$tip.label),
#         as.character(meta.final_ok_ok_withinstudyID$Species))
# # No error or inconsistencies found.
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# # 5th February 2025
# save(my_tree_withinstudyID, file = "data/outputs/phylogenetic_files/tree_withinstudyID.Rdata")

# We can now load the saved tree
load("data/outputs/phylogenetic_files/tree_withinstudyID.Rdata") #my_tree_withinstudyID

# # Compute branch lengths of tree
# phylo_branch_withinstudyID <- compute.brlen(my_tree_withinstudyID, method = "Grafen", power = 1)
# 
# # Check if tree is ultrametric
# is.ultrametric(phylo_branch_withinstudyID)
# # TRUE
# 
# # Matrix to be included in the models
 # phylo_cor_withinstudyID <- vcv(phylo_branch_withinstudyID, cor = T)
# 
# # Finally, save matrix for future analyses to speed up and allow full reproducibility
# save(phylo_cor_withinstudyID, file = "data/outputs/phylogenetic_files/phylo_cor_withinstudyID.Rdata")

# we can now load the saved matrix
load("data/outputs/phylogenetic_files/phylo_cor_withinstudyID.Rdata") #phylo_cor_withinstudyID

# Creating a duplicate species variable for the phylogenetic analysis
meta.final_ok_ok_withinstudyID$Species_phylo_withinstudyID <- meta.final_ok_ok_withinstudyID$Species

# VARIANCE-COVARIANCE MATRIX:
VCV_ESVar_withinstudyID <- matrix(0, nrow = nrow(meta.final_ok_ok_withinstudyID), 
                                  ncol = nrow(meta.final_ok_ok_withinstudyID))

# Names rows and columns for each obsID
rownames(VCV_ESVar_withinstudyID) <- meta.final_ok_ok_withinstudyID[, "EffectID"]
colnames(VCV_ESVar_withinstudyID) <- meta.final_ok_ok_withinstudyID[, "EffectID"]

# Finds effect sizes that come from the same study
shared_coord_withinstudyID <- which(meta.final_ok_ok_withinstudyID[, "StudyID"] %in% 
                                      meta.final_ok_ok_withinstudyID[duplicated(meta.final_ok_ok_withinstudyID[, "StudyID"]), 
                                                                     "StudyID"] == TRUE)

combinations_withinstudyID <- do.call("rbind", tapply(shared_coord_withinstudyID, 
                                                      meta.final_ok_ok_sr[shared_coord_withinstudyID, "StudyID"], 
                                                      function(x) t(utils::combn(x, 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1 : dim(combinations_withinstudyID)[1]) {
  p1_withinstudyID <- combinations_withinstudyID[i, 1]
  p2_withinstudyID <- combinations_withinstudyID[i, 2]
  p1_p2_cov_withinstudyID <- 0.5 * sqrt(meta.final_ok_ok_withinstudyID[p1_withinstudyID, "cor_var"]) * 
    sqrt(meta.final_ok_ok_withinstudyID[p2_withinstudyID, "cor_var"])
  VCV_ESVar_withinstudyID[p1_withinstudyID, p2_withinstudyID] <- p1_p2_cov_withinstudyID
  VCV_ESVar_withinstudyID[p2_withinstudyID, p1_withinstudyID] <- p1_p2_cov_withinstudyID
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(VCV_ESVar_withinstudyID) <- meta.final_ok_ok_withinstudyID[, "cor_var"]

# # In case you want to visually double check the matrix outside of R
# write.csv(VCV_ESVar_withinstudyID, 'data/outputs/variance-covariance_matrices/VCV_ESVar_withinstudyID.csv')


# STATISTICAL MODEL:
meta.model_fitness_all <- rma.mv(cor,
                                 VCV_ESVar_withinstudyID,
                                 mods = ~ Fitness_trait_all_group - 1,
                                 random = list(~ 1 | StudyID,
                                               ~ 1 | LaboratoryID,
                                               ~ 1 | PopulationID,
                                               ~ 1 | Species,
                                               ~ 1 | Species_phylo_withinstudyID,
                                               ~ 1 | EffectID),
                                 method = "REML",
                                 R = list(Species_phylo_withinstudyID = phylo_cor_withinstudyID),
                                 test = "t",
                                 data = meta.final_ok_ok_withinstudyID)

save(meta.model_fitness_all, file = "data/outputs/statistical_models/meta.model_fitness_all.Rdata")
load("data/outputs/statistical_models/meta.model_fitness_all.Rdata") 

# Printing the summary results of the model
print(meta.model_fitness_all, digits = 3)

# Calculate marginal R2 with r2_ml
R2.meta.model_fitness_all <- r2_ml(meta.model_fitness_all)
round(R2.meta.model_fitness_all * 100, 1)

# TABLE WITH RESULTS:
results_meta.model_fitness_all <- orchaRd::mod_results(meta.model_fitness_all, 
                                                       mod = "Fitness_trait_all_group", 
                                                       group = "StudyID", 
                                                       subset = TRUE)

round(as.data.frame(results_meta.model_fitness_all[[1]])[, c(2:6)], 3)


fig_meta.model_fitness_all <- orchaRd::orchard_plot(meta.model_fitness_all, 
                                                    mod = "Fitness_trait_all_group", 
                                                    group = "StudyID", 
                                                    xlab = "Effect size",
                                                    trunk.size = 2,
                                                    branch.size = 2,
                                                    twig.size = 1, 
                                                    colour = FALSE)



################################################################################
# Exploratory analysis - non pre-registered - TO BETTER UNDERSTAND THE 
# HETEROGENEINITY FOUND IN THE RANDOM EFFECT 'PHYLOGENETIC RELATIONSHIPS"

# We found that 'phylogeny' explains a high amount of heterogeneity. On option
# is that this is explained because species differ in some key aspects that are
# phylogenetically correlated.

################################################################################

# STATISTICAL MODEL:
meta.model_fitness_altricial_precocial <- rma.mv(cor, 
                                                 VCV_ESVar, 
                                                 mods = ~ Developmental_mode - 1,
                                                 random = list(~ 1 | StudyID,
                                                               ~ 1 | LaboratoryID,
                                                               ~ 1 | PopulationID,
                                                               ~ 1 | Species,
                                                               ~ 1 | Species_phylo,
                                                               ~ 1 | EffectID), 
                                                 method = "REML",
                                                 R = list(Species_phylo = phylo_cor),
                                                 test = "t", 
                                                 data = meta.final_ok_ok)

#save(meta.model_fitness_altricial_precocial, file = "data/outputs/statistical_models/meta_model_fitness_altricial_precocial.Rdata")
load("data/outputs/statistical_models/meta_model_fitness_altricial_precocial.Rdata") 

# Printing the summary results of the model
print(meta.model_fitness_altricial_precocial, digits = 3)

# Calculate marginal R2 with r2_ml
R2.meta.model_fitness_altricial_precocial <- r2_ml(meta.model_fitness_altricial_precocial)
round(R2.meta.model_fitness_altricial_precocial * 100, 1)

# TABLE WITH RESULTS:
results_meta.model_fitness_altricial_precocial <- orchaRd::mod_results(meta.model_fitness_altricial_precocial, 
                                                                       mod = "Developmental_mode", 
                                                                       group = "StudyID", 
                                                                       subset = TRUE)

round(as.data.frame(results_meta.model_fitness_altricial_precocial[[1]])[, c(2:6)], 3)


fig_meta.model_fitness_altricial_precocial <- orchaRd::orchard_plot(meta.model_fitness_altricial_precocial, 
                                                                    mod = "Developmental_mode", 
                                                                    group = "StudyID", 
                                                                    xlab = "Effect size",
                                                                    trunk.size = 2,
                                                                    branch.size = 2,
                                                                    twig.size = 1, 
                                                                    colour = FALSE)

