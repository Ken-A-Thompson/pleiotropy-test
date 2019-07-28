# to do: different parent SD.

# Fisher's model test
# author: ken a. thompson
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### install & load packages ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# note for users during review process: I am saving Figures to a local overleaf directory because it's synced to my manuscript preparation system. this won't work on your machine but if you do really want to save rather than view the figures just change to a directory that exists on your machine.

# uncomment to install
# install.packages('devtools')
# install.packages('ggfortify')
# install.packages('tidyverse')
# install.packages('phytools')
# install.packages('pspearman')

# load 'em
library(cowplot)
library(ggfortify)
library(pspearman) # for the spearman test of correlation
# library(hablar) # rationalize function
library(tidyverse)
# for phytools
library(devtools)
# install_github("liamrevell/phytools")
library(phytools)

#%%%%%%%%%%%%%%%%%#
#### load data ####
#%%%%%%%%%%%%%%%%%#

GenDist_All <- read.csv(file = 'data/GenDist_Data_All.csv')
nis_traits_std <- read_csv(file = 'data/nis_traits_std.csv')
nis_traits_SD <- read_csv(file = 'data/nis_traits_SD.csv')
nis_traits_different_or_not <- read_csv('data/nis_traits_different_or_not.csv')
nis_species_intra_inter <- read_csv('data/nis_species_intra_inter.csv')
phyloT_tree <- read.newick('data/2019-03-02_phyloT_generated_tree_1551546803_newick.txt') # tree from PhyloT; more spp no BL
# timetree_tree
timetree_dt_data <- read.csv(file = 'data/FILLED_timetree_data.csv')


# simulation data
no_p_sim_data <- read_csv(file = 'simulations/data/no_p_phenotypes_n2_K1000_pmut0.1_alpha0.1_B2_u0.001.csv', col_names = c('ind_no', 'rep', 'nmuts', 'dist', 't1', 't2'))
p_lg_sim_data <- read_csv(file = 'simulations/data/p_lg_phenotypes_n2_K1000_pmut0.1_alpha0.1_B2_u0.001.csv', col_names = c('ind_no', 'rep', 'nmuts', 'dist', 't1', 't2'))
p_sm_sim_data <- read_csv(file = 'simulations/data/p_sm_phenotypes_n2_K1000_pmut0.1_alpha0.0_B2_u0.100.csv', col_names = c('ind_no', 'rep', 'nmuts', 'dist', 't1', 't2'))

#%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%#
#### functions ####
#%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%#

# creates summary statistics for particular trait types binning by stats AND SD
alternative_filtering_function_stats_SD <- function(trait_type){
  
  # create a dataset that has studies with divergent and non divergent traits 
  df1 <- nis_traits_std %>% 
    dplyr::filter(Species_or_CrossType == "F1") %>% # start by looking only at F1; Parent data are associated
    left_join(., nis_traits_different_or_not) %>% # now bring in the dataset asking if traits are different or not
    filter(TraitType %in% trait_type) %>% # only include morphological & pigment traits
    hablar::rationalize(max_SD_diff) %>%  # convert infinite to NA
    group_by(StudyID, Cross_ID) %>% # grouping variables for pipe
    select(max_SD_diff, parents_different_stats_and_SD) %>% 
    group_by(StudyID, Cross_ID, parents_different_stats_and_SD) %>% 
    summarise(mean_SD_diff = mean(log1p(max_SD_diff), na.rm = T)) %>% 
    spread(key = parents_different_stats_and_SD, value = mean_SD_diff) %>% 
    select(-`<NA>`) %>% 
    rename(mean_sds_diff_non_divergent_traits = `FALSE`,
           mean_sds_diff_divergent_traits = `TRUE`) %>% 
    na.omit()
  
  # we need to reduce this dataset to studies that have F1s and F2s.
  df2 <- nis_traits_SD %>% 
    left_join(., nis_traits_different_or_not) %>%
    filter(StudyID %in% df1$StudyID) %>% 
    filter(TraitType %in% trait_type) %>% # only include morphological & pigment traits
    mutate(Species_or_CrossType = ifelse(!Species_or_CrossType %in% c("F1", "F2", "BC"), NA, Species_or_CrossType)) %>% 
    mutate(Species_or_CrossType = coalesce(Species_or_CrossType, ParentID)) %>% 
    filter(parents_different_stats_and_SD == F) %>% 
    group_by(StudyID, Cross_ID, Cross_Dir, Sex, TraitNo) %>% 
    select(Species_or_CrossType, TraitNo, SD) %>% 
    # filter(Species_or_CrossType %in% c("F1", "F2")) %>%
    group_by(StudyID, Cross_ID, Sex, Species_or_CrossType, TraitNo) %>% 
    summarise(SD = mean(SD, na.rm = T)) %>% 
    group_by(StudyID, Cross_ID, Species_or_CrossType, TraitNo) %>% 
    summarise(SD = mean(SD, na.rm = T)) %>% # average across sexes
    select(Species_or_CrossType, SD, TraitNo) %>%
    spread(key = Species_or_CrossType, value = SD) %>% 
    mutate(diff_segregation_variance_alternative = F2^2 - F1^2) %>%
    mutate(diff_segregation_variance = 4 * F2^2 - ((2 * (F1^2) + A^2 + B^2 ))) %>% 
    mutate(segregation_variance_alternative = F2^2 / F1^2) %>%
    mutate(segregation_variance = 4 * F2^2 / ((2 * (F1^2) + A^2 + B^2 ))) %>% 
    group_by(StudyID, Cross_ID) %>% 
    summarise(mean_segvar_non_diff = mean(log1p(segregation_variance), na.rm = T),
              mean_segvar_non_diff_alt = mean(log1p(segregation_variance_alternative), na.rm = T),
              mean_segvar_eq1_diff = mean(log1p(diff_segregation_variance), na.rm = T),
              mean_segvar_alt_diff = mean(log1p(diff_segregation_variance_alternative), na.rm = T)
    ) %>% 
    left_join(., df1) %>% 
    select(StudyID, Cross_ID, mean_segvar_non_diff, mean_segvar_non_diff_alt, mean_sds_diff_divergent_traits, mean_sds_diff_non_divergent_traits, mean_segvar_alt_diff, mean_segvar_eq1_diff)  %>% 
    na.omit() %>% 
    mutate(studycross = paste(StudyID, Cross_ID, sep = " "))
  return(df2)
}

# creates summary statistics for particular trait types binning by stats ONLY
alternative_filtering_function_stats <- function(trait_type){
  
  # create a dataset that has studies with divergent and non divergent traits 
  df1 <- nis_traits_std %>% 
    dplyr::filter(Species_or_CrossType == "F1") %>% # start by looking only at F1; Parent data are associated
    left_join(., nis_traits_different_or_not) %>% # now bring in the dataset asking if traits are different or not
    filter(TraitType %in% trait_type) %>% # only include morphological & pigment traits
    hablar::rationalize(max_SD_diff) %>%  # convert infinite to NA
    group_by(StudyID, Cross_ID) %>% # grouping variables for pipe
    select(max_SD_diff, parents_different_stats) %>% 
    group_by(StudyID, Cross_ID, parents_different_stats) %>% 
    summarise(mean_SD_diff = mean(log1p(max_SD_diff), na.rm = T)) %>% 
    spread(key = parents_different_stats, value = mean_SD_diff) %>% 
    select(-`<NA>`) %>% 
    rename(mean_sds_diff_non_divergent_traits = `FALSE`,
           mean_sds_diff_divergent_traits = `TRUE`) %>% 
    na.omit()
  
  # we need to reduce this dataset to studies that have F1s and F2s.
  df2 <- nis_traits_SD %>% 
    left_join(., nis_traits_different_or_not) %>%
    filter(StudyID %in% df1$StudyID) %>% 
    filter(TraitType %in% trait_type) %>% # only include morphological & pigment traits
    mutate(Species_or_CrossType = ifelse(!Species_or_CrossType %in% c("F1", "F2", "BC"), NA, Species_or_CrossType)) %>% 
    mutate(Species_or_CrossType = coalesce(Species_or_CrossType, ParentID)) %>% 
    filter(parents_different_stats == F) %>% 
    group_by(StudyID, Cross_ID, Cross_Dir, Sex, TraitNo) %>% 
    select(Species_or_CrossType, TraitNo, SD) %>% 
    # filter(Species_or_CrossType %in% c("F1", "F2")) %>%
    group_by(StudyID, Cross_ID, Sex, Species_or_CrossType, TraitNo) %>% 
    summarise(SD = mean(SD, na.rm = T)) %>% 
    group_by(StudyID, Cross_ID, Species_or_CrossType, TraitNo) %>% 
    summarise(SD = mean(SD, na.rm = T)) %>% # average across sexes
    select(Species_or_CrossType, SD, TraitNo) %>%
    spread(key = Species_or_CrossType, value = SD) %>% 
    mutate(diff_segregation_variance_alternative = F2^2 - F1^2) %>%
    mutate(diff_segregation_variance = 4 * F2^2 - ((2 * (F1^2) + A^2 + B^2 ))) %>% 
    mutate(segregation_variance_alternative = F2^2 / F1^2) %>%
    mutate(segregation_variance = 4 * F2^2 / ((2 * (F1^2) + A^2 + B^2 ))) %>% 
    group_by(StudyID, Cross_ID) %>% 
    summarise(mean_segvar_non_diff = mean(log1p(segregation_variance), na.rm = T),
              mean_segvar_non_diff_alt = mean(log1p(segregation_variance_alternative), na.rm = T),
              mean_segvar_eq1_diff = mean(log1p(diff_segregation_variance), na.rm = T),
              mean_segvar_alt_diff = mean(log1p(diff_segregation_variance_alternative), na.rm = T)
    ) %>% 
    left_join(., df1) %>% 
    select(StudyID, Cross_ID, mean_segvar_non_diff, mean_segvar_non_diff_alt, mean_sds_diff_divergent_traits, mean_sds_diff_non_divergent_traits, mean_segvar_alt_diff, mean_segvar_eq1_diff)  %>% 
    na.omit() %>% 
    mutate(studycross = paste(StudyID, Cross_ID, sep = " "))
  return(df2)
}
#%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%#
#### data processing ####
#%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%#

# create a dataset that has studies with divergent and non divergent traits 
divergence_and_transgression_df <- nis_traits_std %>% 
  dplyr::filter(Species_or_CrossType == "F1") %>% # start by looking only at F1; Parent data are associated
  left_join(., nis_traits_different_or_not) %>% # now bring in the dataset asking if traits are different or not
  # filter(parents_different_stats == T) %>% filter to restrict dataset to traits that are different at P < 0.05
  # filter(TraitType %in% c("Morphology", "Pigment")) %>% # only include morphological & pigment traits
  filter(TraitType %in% c("Morphology")) %>% # only include morphological & pigment traits
  # filter(TraitType != "Physiology") %>%
  hablar::rationalize(max_SD_diff) %>%  # convert infinite to NA
  group_by(StudyID, Cross_ID) %>% # grouping variables for pipe
  select(max_SD_diff, parents_different_stats) %>% 
  group_by(StudyID, Cross_ID, parents_different_stats) %>% 
  summarise(mean_SD_diff = mean(log1p(max_SD_diff), na.rm = T)) %>% 
  spread(key = parents_different_stats, value = mean_SD_diff) %>% 
  select(-`<NA>`) %>% 
  rename(mean_sds_diff_non_divergent_traits = `FALSE`,
         mean_sds_diff_divergent_traits = `TRUE`) %>% 
  na.omit()

# we need to reduce this dataset to studies that have F1s and F2s.
divergence_and_transgression_df_segvar <- nis_traits_SD %>% 
  # filter(Species_or_CrossType %in% c("F1", "F2")) %>%
  left_join(., nis_traits_different_or_not) %>%
  filter(StudyID %in% divergence_and_transgression_df$StudyID) %>% 
  # filter(TraitType %in% c("Morphology", "Pigment")) %>% # only include morphological & pigment traits
  filter(TraitType %in% c("Morphology")) %>% # only include morphological & pigment traits
  # filter(TraitType != "Physiology") %>%
  # filter(TraitType != "Behaviour") %>%
  mutate(Species_or_CrossType = ifelse(!Species_or_CrossType %in% c("F1", "F2", "BC"), NA, Species_or_CrossType)) %>% 
  mutate(Species_or_CrossType = coalesce(Species_or_CrossType, ParentID)) %>% 
  filter(parents_different_stats == F) %>% 
  group_by(StudyID, Cross_ID, Cross_Dir, Sex, TraitNo) %>% 
  select(Species_or_CrossType, TraitNo, SD) %>% 
  # filter(Species_or_CrossType %in% c("F1", "F2")) %>%
  group_by(StudyID, Cross_ID, Sex, Species_or_CrossType, TraitNo) %>% 
  summarise(SD = mean(SD, na.rm = T)) %>% 
  group_by(StudyID, Cross_ID, Species_or_CrossType, TraitNo) %>% 
  summarise(SD = mean(SD, na.rm = T)) %>% # average across sexes
  select(Species_or_CrossType, SD, TraitNo) %>%
  spread(key = Species_or_CrossType, value = SD) %>% 
  mutate(segregation_variance_alternative = F2^2 / F1^2) %>%
  mutate(segregation_variance = 4 * F2^2 / ((2 * (F1^2) + A^2 + B^2 ))) %>% # I want to try this other metric
  group_by(StudyID, Cross_ID) %>% 
  summarise(mean_segvar_non_diff = mean(log1p(segregation_variance), na.rm = T),
            mean_segvar_non_diff_alt = mean(log1p(segregation_variance_alternative), na.rm = T)) %>% 
  left_join(., divergence_and_transgression_df) %>% 
  select(StudyID, Cross_ID, mean_segvar_non_diff, mean_segvar_non_diff_alt, mean_sds_diff_divergent_traits, mean_sds_diff_non_divergent_traits)  %>% 
  na.omit() %>% 
  mutate(crosscat = "F2") %>% # add a label because these are all data from F2s
  mutate(studycross = paste(StudyID, Cross_ID, sep = " "))
  
# vector containing list of studies.
study_list_master <- divergence_and_transgression_df_segvar$StudyID

# add in the genetic distance data 
divergence_and_transgression_df_segvar_gendist <- divergence_and_transgression_df_segvar %>% 
  mutate(studycross = paste(StudyID, Cross_ID, sep = " ")) %>% 
  left_join(nis_species_intra_inter, by = "studycross")
  
# create a dataset with each trait treated differently
data_for_dryad <- nis_traits_std %>% 
  filter(TraitType == "Morphology") %>%
  left_join(., nis_traits_different_or_not) %>% 
  filter(StudyID %in% study_list_master)

# generate a list of traits in the studies
study_traits <- data_for_dryad %>% 
  filter(TraitType %in% c("Morphology")) %>%
  select(StudyID, TraitDesc, TraitType) %>% 
  unique()

# list of crosses in the main study
crosses_list <- nis_species_intra_inter %>% 
  filter(studycross %in% divergence_and_transgression_df_segvar$studycross)

# crosses to upload to timetree to get divergence time estimates
crosses_for_timetree <- nis_species_intra_inter %>% 
  filter(intra_inter == "inter") %>% 
  gather(`sp 1`:`sp 2`, key = "species", value = "name")

write.table(crosses_for_timetree$name, file = 'data/timetree.txt',
            row.names = F,
            quote = F,
            col.names = F)

crosses_for_timetree_entry <- nis_species_intra_inter %>% 
  filter(intra_inter == "inter") %>% 
  select(-intra_inter)

write_csv(crosses_for_timetree_entry, path = 'data/timetree_data_entry_template.csv')

# only retain plant/animal info
PA_data <- GenDist_All %>% 
  select(PA, StudyID)

GenDist_all_forMerge <- GenDist_All %>% 
  select(StudyCross, StudyID, GenDist, first, second) %>% 
  unique() %>% 
  mutate(unique_sp = paste(first, second, sep = "_"))

# Study_Species DF
Study_Species_DF <- nis_traits_std %>% 
  ungroup() %>% 
  filter(Parent_Hybrid == "Parent") %>% # only need to consider parents here
  select(StudyID, Cross_ID, ParentID, Species_or_CrossType) %>% 
  unique() %>% 
  mutate(Species_or_CrossType.Short = sub("_", " ", Species_or_CrossType)) %>% 
  separate(Species_or_CrossType.Short, into = c('species', 'col2', 'col3'), sep = "_") %>% 
  group_by(StudyID, Cross_ID) %>% 
  select(-c(col2, col3)) %>% 
  arrange(species) %>%  # alphabetize
  mutate(id = 1:n()) %>% 
  mutate(New_ParentID = ifelse(id == 1, "A", "B")) %>% # recode the numbers to letters for nicer spread
  select(StudyID, Cross_ID, species, New_ParentID) %>% 
  mutate(species = sub(" ", "_", species)) %>% # add underscore
  dplyr::filter(New_ParentID == "A") %>% 
  ungroup() %>% 
  select(StudyID, species)

# divergence time and phenotypic divergence
divergence_time_trait_DF <- nis_traits_std %>% 
  dplyr::filter(Species_or_CrossType == "F1") %>% # start by looking only at F1; Parent data are associated
  left_join(., nis_traits_different_or_not) %>% # now bring in the dataset asking if traits are different or not
  # filter(parents_different_stats == T) %>% filter to restrict dataset to traits that are different at P < 0.05
  hablar::rationalize(max_SD_diff) %>%  # convert infinite to NA
  group_by(StudyID, Cross_ID) %>% # grouping variables for pipe
  summarise(mean_sds_diff_divergent_traits = 
              mean(max_SD_diff[match('TRUE', parents_different_stats)], na.rm = T), # mean divergence of traits that are / are not different
            mean_sds_diff_non_divergent_traits = 
              mean(max_SD_diff[match('FALSE', parents_different_stats)], na.rm = T)) %>% 
  # na.omit() %>% 
  mutate(studycross = paste(StudyID, Cross_ID, sep = " ")) %>% # gendist df has underscore 
  left_join(timetree_dt_data)
  
# all crosses
All_Crosses_GROUPS <- nis_traits_std %>% 
  ungroup() %>% 
  filter(Parent_Hybrid == "Parent") %>% # only need to consider parents here
  select(StudyID, Cross_ID, ParentID, Species_or_CrossType) %>% 
  unique() %>% 
  mutate(Species_or_CrossType.Short = sub("_", " ", Species_or_CrossType)) %>% 
  separate(Species_or_CrossType.Short, into = c('species', 'col2', 'col3'), sep = "_") %>% 
  group_by(StudyID, Cross_ID) %>% 
  select(-c(col2, col3)) %>% 
  arrange(species) %>%  # alphabetize
  mutate(id = 1:n()) %>% 
  mutate(New_ParentID = ifelse(id == 1, "A", "B")) %>% # recode the numbers to letters for nicer spread
  select(StudyID, Cross_ID, species, New_ParentID) %>% 
  spread(value = species, key = New_ParentID) %>% 
  mutate(StudyCross = paste(StudyID, Cross_ID, sep = "_")) %>% # gendist df has underscore 
  mutate(Both_parents_alpha = paste(A, B, sep = " & ")) %>% # join the parents for grouping factor
  group_by(Both_parents_alpha) %>% 
  mutate(group_index = group_indices()) %>% # group indices
  ungroup() %>% 
  select(StudyCross, group_index)

# create the same version of the data for all studies
nis_traits_diff_int_allstudies <- nis_traits_std %>% 
  mutate(studycross = paste(StudyID, Cross_ID, sep = " ")) %>% 
  filter(Species_or_CrossType == "F1") %>% # start by looking only at F1; Parent data are associated
  left_join(., nis_traits_different_or_not) %>% # now bring in the dataset asking if traits are different or not
  # filter(parents_different_stats == T) %>% filter to restrict dataset to traits that are different at P < 0.05
  filter(TraitType %in% c("Morphology")) %>% # only include morphological data & life history
  hablar::rationalize(max_SD_diff) %>%  # convert infinite to NA
  group_by(studycross, StudyID, Cross_ID) %>% # grouping variables for pipe
  summarise(mean_sds_diff_divergent_traits = 
              mean(max_SD_diff[match('TRUE', parents_different_stats)], na.rm = T)) %>% 
  left_join(nis_species_intra_inter)  %>% 
  group_by(StudyID, Cross_ID, intra_inter) %>% 
  mutate(StudyCross = paste(StudyID, Cross_ID, sep = "_")) %>% # gendist df has underscore 
  left_join(., GenDist_All, by = "StudyCross") %>% 
  group_by(StudyCross, StudyID.x, intra_inter) %>% 
  summarise(mean_sds_diff_divergent_traits = mean(mean_sds_diff_divergent_traits, na.rm=T),
            mean_GenDist = mean(GenDist)) %>%
  rename(StudyID = "StudyID.x") %>% 
  left_join(PA_data) %>% # add in plant animal
  mutate(PA = as.character(PA)) %>% 
  mutate(PA = ifelse(StudyID == "Kohn_2001", "Animal", PA)) %>% 
  left_join(., All_Crosses_GROUPS) %>% 
  group_by(group_index, intra_inter, PA) %>% 
  summarise(mean_sds_diff_divergent_traits = mean(mean_sds_diff_divergent_traits, na.rm=T),
            mean_GenDist = mean(mean_GenDist))

#%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%#
#### main text stats ####
#%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%#

# # Spearman's rank-order correlation
spearman.test(divergence_and_transgression_df_segvar$mean_sds_diff_divergent_traits, divergence_and_transgression_df_segvar$mean_segvar_non_diff)
spearman.test(divergence_and_transgression_df_segvar$mean_sds_diff_non_divergent_traits, divergence_and_transgression_df_segvar$mean_segvar_non_diff)

# here are some alternative models that can be uncommented to run
# # main linear models; univariate
# # running with and without transformation so I can show diagnostics
# hypothesis_lm_log <- lm((mean_segvar_non_diff) ~ (mean_sds_diff_divergent_traits), data = divergence_and_transgression_df_segvar)
# # 
# ## summarize the main model
# summary(hypothesis_lm_log)
# 
# # # run a model with non-divergent traits as predictor
# alternative_lm <- lm(mean_segvar_non_diff ~ (mean_sds_diff_non_divergent_traits), data = divergence_and_transgression_df_segvar)
# summary(alternative_lm) # analysis
# 
# # multiple regression
# multiple_lm <-  lm((mean_segvar_non_diff) ~ (mean_sds_diff_non_divergent_traits) + (mean_sds_diff_divergent_traits), data = divergence_and_transgression_df_segvar)
# summary(multiple_lm)
# # 
# multiple_lm <-  lm(mean_sds_diff_divergent_traits ~  mean_sds_diff_non_divergent_traits, data = divergence_and_transgression_df_segvar)
# 
# summary(multiple_lm)

#%%%%%%%%%%%%%%%%%%%%%%%%#
# genetic distance tests # 
#%%%%%%%%%%%%%%%%%%%%%%%%#

# in the dataset currently being considered, is there a difference between inter and intra-specific crosses in phenotypic divergence?
summary(aov(mean_sds_diff_divergent_traits ~ intra_inter, data = divergence_and_transgression_df_segvar_gendist))

#%%%%%%%%%%%%%%%%%%%%%%%%#
# pylogenetic signal # 
#%%%%%%%%%%%%%%%%%%%%%%%%#

# # phytools uses SPECIES as rownames (ick)
#
Study_Species_DF_FULL <- divergence_and_transgression_df_segvar %>%
  left_join(., Study_Species_DF) %>%
  unique() %>%
  ungroup() %>%
  mutate(rn = row_number()) %>% # this is not great right now but works; deletes duplicates selz
  filter(rn %in% c(1:12, 15:17)) %>%
  ungroup()

# trying treeplyr; will move up if I like it
# install.packages('treeplyr')
# library(treeplyr)

# load TIMETREE
# timetree_tree <- read.tree('!FGM_Test_Offshoot/data/timetree.nwk')
# try make.treedata from treeplyr
# td <- make.treedata(tree = timetree_tree, data = )

mean_segvar_non_diff_DF <- Study_Species_DF_FULL %>%
  select(mean_segvar_non_diff)

mean_p_div_DF <- Study_Species_DF_FULL %>%
  select(mean_sds_diff_divergent_traits)

# stupid but have to change rownames
rownames(mean_segvar_non_diff_DF) <- Study_Species_DF_FULL$species
rownames(mean_p_div_DF) <- Study_Species_DF_FULL$species

species_to_retain <- as.vector(rownames(mean_segvar_non_diff_DF))

# THIS SEEMS TO WORK!
# create data input for phytools::phloysig
trait_segvar <- as.vector(t(log(mean_segvar_non_diff_DF[,1])))
names(trait_segvar) <- rownames(mean_segvar_non_diff_DF)

# divergence
# mean_p_div_DF2 <-F
# clean this up; I get sloppy in base!
trait_div <- as.vector(t(log(Study_Species_DF_FULL[,4])))
names(trait_div) <- rownames(mean_p_div_DF)


pruned.tree <- drop.tip(phyloT_tree, setdiff(phyloT_tree$tip.label, species_to_retain));
# # set branch lengths to random
pruned.tree$edge.length <- rep(0.2, 2 * phyloT_tree$Nnode)

phylosig(pruned.tree, trait_segvar, test = TRUE)
phylosig(pruned.tree, trait_div, test = TRUE)

# also can do lambda
# phylogenetic signal for trait segregation variance (log)
phylosig(pruned.tree, trait_segvar, method="lambda",test=TRUE)

# phylogenetic signal for parental divergence
phylosig(pruned.tree, trait_div, method="lambda",test=TRUE)

#%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%#
#### figures ####
#%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%#

# load the ggplot figure theme
theme_KT_FGM <-
  theme(
    aspect.ratio = 1.0,
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(size = 1),
    axis.line.x = element_line(color = "black", size = 1),
    axis.line.y = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.title.y = element_text(vjust = 0.2, size = 14),
    axis.title.x = element_text(vjust = 0.1, size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none")

#%%%%%%%%%%%%%%%#
# MAIN ANALYSIS #
#%%%%%%%%%%%%%%%#

# here's the plot; appears it doesn't predict F1 mean phenotype... what about segregation variance?
segregation_variance_phenotype <- 
  ggplot(divergence_and_transgression_df_segvar, 
         aes(x = mean_sds_diff_divergent_traits, y = mean_segvar_non_diff)) +
  xlab("mean log(phenotypic distance [SDs]\nof divergent traits) between parents") +
  ylab("mean log(segregation variance) of\nnon-divergent traits in hybrids") +
  geom_smooth(method = "loess", se = F, colour = "black") +
  geom_point() +
  theme_KT_FGM + 
  theme(aspect.ratio = 3/4)
# segregation_variance_phenotype

segregation_variance_phenotype_non_divergent_fig <- ggplot(divergence_and_transgression_df_segvar, 
         aes(x = log(mean_sds_diff_non_divergent_traits), y = log(mean_segvar_non_diff))) +
  xlab("ln(phenotypic distance (SDs)\nof non-divergent traits + 1)") +
  ylab("ln(segregation variance of\nnon-divergent traits)") +
  geom_point() +
  # geom_smooth(method = "lm", se = T, colour = "black") +
  theme_KT_FGM + 
  theme(aspect.ratio = 3/4) 
  # annotate("text", label = expression(paste(italic("P"), " = 0.0265; ",
  #                                           italic("r")^2, " = 0.374")), x = 0.2, y = 1.5, size = 6)
# segregation_variance_phenotype_non_divergent_fig

# summary(lm(log(mean_segvar_non_diff) ~ log1p(mean_sds_diff_non_divergent_traits), divergence_and_transgression_df_segvar))

ggsave(segregation_variance_phenotype, filename = '../../../../Apps/Overleaf/pleiotropy_ms/Figures/pleiotropy_Figure_2.pdf', height = 4, width = 5)
ggsave(segregation_variance_phenotype_non_divergent_fig, filename = '../../../../Apps/Overleaf/pleiotropy_ms/Figures/non_divergent_traits.pdf', height = 4, width = 5)

#%%%%%%%%%%%%%#
# SIMULATIONS #
#%%%%%%%%%%%%%#

no_p_sim_data_summary <- no_p_sim_data %>% 
  group_by(rep, dist) %>% 
  summarise(var_t1 = var(t1),
            var_t2 = var(t2))

p_lg_sim_data_summary <- p_lg_sim_data %>% 
  group_by(rep, dist) %>% 
  summarise(var_t1 = var(t1),
            var_t2 = var(t2))

p_sm_sim_data_summary <- p_sm_sim_data %>% 
  group_by(rep, dist) %>% 
  summarise(var_t1 = var(t1),
            var_t2 = var(t2))

divergent_p_lg <- 
  ggplot(p_lg_sim_data_summary, 
         aes(x = dist, y = var_t1)) +
  xlab("distance to optimum") +
  ylab("segregation variance\n(divergent trait)") +
  geom_point() +
  geom_smooth(method = "loess", se = T, colour = "black") +
  ggtitle('pleiotropy; finite') +
  theme_KT_FGM

NONdivergent_p_lg <- 
  ggplot(p_lg_sim_data_summary, 
         aes(x = dist, y = var_t2)) +
  xlab("distance to optimum") +
  ylab("segregation variance\n(non-divergent trait)") +
  geom_point() +
  geom_smooth(method = "loess", se = T, colour = "black") +
  ggtitle('pleiotropy; finite') +
  theme_KT_FGM


divergent_p_sm <- 
  ggplot(p_sm_sim_data_summary, 
         aes(x = dist, y = var_t1)) +
  xlab("distance to optimum") +
  ylab("segregation variance\n(divergent trait)") +
  geom_point() +
  # ylim(c(-0.001, 0.001)) +
  # geom_smooth(method = "loess", se = T, colour = "black") +
  ggtitle('pleiotropy; "infinitesimal"') +
  theme_KT_FGM


NONdivergent_p_sm <- 
  ggplot(p_sm_sim_data_summary, 
         aes(x = dist, y = var_t2)) +
  xlab("distance to optimum") +
  ylab("segregation variance\n(non-divergent trait)") +
  geom_point() +
  # geom_smooth(method = "loess", se = T, colour = "black") +
  ggtitle('pleiotropy; "infinitesimal"') +
  theme_KT_FGM

divergent_no_p <- 
  ggplot(no_p_sim_data_summary, 
         aes(x = dist, y = var_t1)) +
  xlab("distance to optimum") +
  ylab("segregation variance\n(divergent trait)") +
  geom_point() +
  geom_smooth(method = "loess", se = T, colour = "black") +
  ggtitle('no pleiotropy; finite') +
  theme_KT_FGM

NONdivergent_no_p <- 
  ggplot(no_p_sim_data_summary, 
         aes(x = dist, y = var_t2)) +
  xlab("distance to optimum") +
  ylab("segregation variance\n(non-divergent trait)") +
  geom_point() +
  # geom_smooth(method = "loess", se = T, colour = "black") +
  ggtitle('no pleiotropy; finite') +
  theme_KT_FGM

# plot it
sims_fig <- plot_grid(NONdivergent_p_sm, NONdivergent_no_p, NONdivergent_p_lg, 
                      divergent_p_sm, divergent_no_p, divergent_p_lg, labels = "AUTO")
# save it
ggsave(sims_fig, filename = '../../../../Apps/Overleaf/pleiotropy_ms/Figures/sims_fig.pdf', height = 7, width = 10)

#%%%%%%%%%%%%%#
# DIAGNOSTICS #
#%%%%%%%%%%%%%#

hypothesis_lm <- lm(mean_segvar_non_diff ~ mean_sds_diff_divergent_traits, data = divergence_and_transgression_df_segvar)
summary(hypothesis_lm)
# install.packages('lmtest')
library(lmtest)

bptest(hypothesis_lm)

# diagnostic plots
diagnostics_plot <- autoplot(hypothesis_lm)

# save em
ggsave(diagnostics_plot, filename = '../../../../Apps/Overleaf/pleiotropy_ms/Figures/diagnostics_fig.pdf', height = 5, width = 5)

#%%%%%%%%%%%%%%%%%%%%#
# NON-DIVERGENT DATA #
#%%%%%%%%%%%%%%%%%%%%#

segregation_variance_phenotype_non_divergent_fig <- 
  ggplot(divergence_and_transgression_df_segvar, 
         aes(x = (mean_sds_diff_non_divergent_traits), y = (mean_segvar_non_diff))) +
  xlab("phenotypic distance (SDs)\nof non-divergent traits") +
  ylab("segregation variance of\nnon-divergent traits") +
  geom_point() +
  # geom_smooth(method = "lm", se = F, colour = "black") +
  theme_KT_FGM + 
  theme(aspect.ratio = 3/4)
# segregation_variance_phenotype_non_divergent_fig_RAW

ggsave(segregation_variance_phenotype_non_divergent_fig, filename = '../../../../Apps/Overleaf/pleiotropy_ms/Figures/non_divergent_fig.pdf', height = 4, width = 5)

#%%%%%%%%%%#
# GEN DIST #
#%%%%%%%%%%#
Cross.Names <- read.csv('NIS_Analysis/data/List_of_Species_Crosses.csv')
# violin plot reduced dataset
Intra_Inter_Reduced_Fig <- ggplot(divergence_and_transgression_df_segvar_gendist, aes(x = intra_inter, y = log(mean_sds_diff_divergent_traits), fill = intra_inter)) +
  geom_violin() +
  geom_jitter(width = 0.05, alpha = 0.3) +
  stat_summary(fun.y=mean, geom="point", shape=19, size=3) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", size = 1, aes(width = 0.5)) +
  scale_fill_manual(values=c("tomato", "lightseagreen")) +
  ylab("ln(phenotypic distance [SDs]\nof divergent traits)") + 
  xlab("cross type") +
  theme_KT_FGM

# some studies (e.g., David_2002) use the same species for all crosses and 
# plot the thing for all studies
Intra_Inter_All_Fig <- ggplot(unique(nis_traits_diff_int_allstudies), aes(x = intra_inter, y = log(mean_sds_diff_divergent_traits), fill = intra_inter)) +
  geom_violin() +
  geom_jitter(width = 0.05, alpha = 0.3) +
  stat_summary(fun.y=mean, geom="point", shape=19, size=3) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", size = 1, aes(width = 0.5)) +
  scale_fill_manual(values=c("tomato", "lightseagreen")) +
  ylab("ln(phenotypic distance [SDs]\nof divergent traits)") + 
  xlab("cross type") +
  theme_KT_FGM 

# 
summary(aov(log(mean_sds_diff_divergent_traits) ~ intra_inter, data = unique(nis_traits_diff_int_allstudies)))


Gendist_Cont_Fig <- ggplot(nis_traits_diff_int_allstudies, aes(x = log(mean_GenDist), y = log(mean_sds_diff_divergent_traits), colour = PA)) +
  geom_point() +
  ylab("log(genetic divergence)") + 
  xlab("ln(phenotypic distance [SDs]\nof divergent traits)") +
  scale_colour_manual(values = c("brown", "green")) +
  theme_KT_FGM

summary(lm(log(mean_sds_diff_divergent_traits) ~ log(mean_GenDist), nis_traits_diff_int_allstudies))

Divergence_Time_Fig <- ggplot(divergence_time_trait_DF, aes(x = estimated_time_timetree, y = log(mean_sds_diff_divergent_traits))) +
  geom_point() + 
  geom_smooth() + 
  ylab("ln(phenotypic distance [SDs]\nof divergent traits)") + 
  xlab("divergence time (MYA)") + 
  theme_KT_FGM

summary(lm(log(mean_sds_diff_divergent_traits) ~ estimated_time_timetree, divergence_time_trait_DF))

gendist_fig <- plot_grid(Intra_Inter_Reduced_Fig, Intra_Inter_All_Fig, Gendist_Cont_Fig, Divergence_Time_Fig, labels = "AUTO", ncol = 2)
ggsave(gendist_fig, filename = '../../../../Apps/Overleaf/pleiotropy_ms/Figures/gendist_fig.pdf', height = 6, width = 6)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### ALTERNATIVE CHOICES ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# create a dataset that has studies with divergent and non divergent traits 
stats_or_SD_morph_df <- alternative_filtering_function_stats_SD(trait_type = c("Morphology"))
stats_morph_df <- alternative_filtering_function_stats(trait_type = c("Morphology"))
stats_all_df <- alternative_filtering_function_stats(trait_type = c("Morphology", "Life_history", "Behaviour", "Chemical", "Pigment", "Physiology"))


# results for table S1
## reference
spearman.test(stats_morph_df$mean_sds_diff_divergent_traits, stats_morph_df$mean_segvar_non_diff)

## different binning
spearman.test(stats_or_SD_morph_df$mean_sds_diff_divergent_traits, stats_or_SD_morph_df$mean_segvar_non_diff)

##  eq1 but difference
spearman.test(stats_morph_df$mean_sds_diff_divergent_traits, stats_morph_df$mean_segvar_eq1_diff)

## ratio f2 to f1
spearman.test(stats_morph_df$mean_sds_diff_divergent_traits, stats_morph_df$mean_segvar_non_diff_alt)

## difference f2 - f1
spearman.test(stats_morph_df$mean_sds_diff_divergent_traits, stats_morph_df$mean_segvar_alt_diff)

# all traits
spearman.test(stats_all_df$mean_sds_diff_divergent_traits, stats_all_df$mean_segvar_alt_diff)


