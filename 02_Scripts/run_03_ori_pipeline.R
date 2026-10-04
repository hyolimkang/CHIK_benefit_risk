# Driver: runs the prerequisite chain 03_brazil_all_draws_ori_v3.R expects
# (its own header says "source setup.R and setup_age_props.R before this
# script" -- it has no library()/pacman::p_load() of its own, so a cold
# `Rscript 03_brazil_all_draws_ori_v3.R` fails on missing package functions
# like tibble()). Sourcing 01_setup.R here also means lhs_sample is already
# fresh in memory before 03 checks lhs_has_age_props(), so it never falls
# back to loading 01_Data/lhs_sample.RData at all in this run.
setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
source("02_Scripts/01_setup.R")
source("02_Scripts/02_setup_age_props.R")
source("02_Scripts/03_brazil_all_draws_ori_v3.R")
