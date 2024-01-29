# Script name: 002_select_participants.R
# Project: WCST project
# Script purpose: select the participants incuded in the project.
# @author: Corrado Caudek <corrado.caudek@unifi.it>
# Date Created: Wed Jun  1 17:05:59 2022
# Last Modified Date: Thu Jun  2 09:02:44 2022
#
# 👉 
# This script is used for reading the raw data of both controls and patients.
# Change GROUP as "controls" or "patients" to generate the RDS file for each
# group.


# Prelims
library("here")
library("tidyverse")
library("stringi")
library("readxl")
library("janitor")

# source(here("scripts", "01_prelims.R"))
# source(here("libraries", "func_eating.R"))


patients_df <- readRDS(
  here("data", "interim", "task_switching", "patients_task_switching_data.rds")
)

controls_df <- readRDS(
  here("data", "interim", "task_switching", "controls_task_switching_data.rds")
)

# Import the subject codes for the actual participants to the project.
source(here::here("data", "list_participants.R"))

an_df <- patients_df[patients_df$subj_id %in% participants_list[[1]], ]
length(unique(an_df$subj_id))
an_df$group <- "AN"

hc_df <- controls_df[controls_df$subj_id %in% participants_list[[2]], ]
length(unique(hc_df$subj_id))
hc_df$group <- "HC"

ri_df <- controls_df[controls_df$subj_id %in% participants_list[[3]], ]
length(unique(ri_df$subj_id))
ri_df$group <- "RI"

temp <- rbind(an_df, hc_df)
task_switch_df <- rbind(temp, ri_df)

task_switch_df |> 
  group_by(group) |> 
  summarize(
    n = n_distinct(subj_id)
  )

rio::export(
  task_switch_df,
  here::here("data", "processed", "task_switching", "task_switching_data.csv")
)


#---- eof 


