# Script name: 001_read_data.R
# Project: Eating disorders Montecatini
# Script purpose: read task switching data.
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


# Remove from the folder with the files of each participant these two files 
# because two subjects have repeated the task twice:
# "taskswitching1.2021-03-15-2022.data.5f7ff571-0c46-468b-8a65-73f1a1af3b7c.txt
# "taskswitching1.2020-12-02-1054.data.00517945-bf6d-454c-b1ca-cf9ca94c173a.txt"

GROUP <- "controls" # "patients" "controls"


# Get participant's identifier --------------------------------------------

d <- read_excel(here::here("data", "raw", "task_switching", GROUP, "data.xlsx"))

d$mese_c <- ifelse(
  d$`mese:1` < 10, stri_join("0", as.character(d$`mese:1`), sep=''), as.character(d$`mese:1`)
  )

d$giorno_c <- ifelse(
  d$`giorno:1` < 10, 
  stri_join("0", as.character(d$`giorno:1`), sep=''), 
  as.character(d$`giorno:1`)
)

d$cellulare_c <- ifelse(
  d$`cellulare:1` < 100, 
  stri_join("0", as.character(d$`cellulare:1`), sep=''), 
  as.character(d$`cellulare:1`)
)

d$sex <- ifelse(d$`sesso:1` == 1, "f",
                ifelse(d$`sesso:1` == 2, "m", NA))

d$subj_id <- tolower(
  stri_join(d$`nome:1`, d$`cognome:1`, d$`anno:1`, 
            d$mese_c, d$giorno_c, d$cellulare_c, d$sex, 
            sep='_')
)


# Get list for all participants -------------------------------------------

subj_file_dat <- d$`esperimento:1`

df_fake <- data.frame(matrix(ncol = 10, nrow = 0))
col_names <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10")
colnames(df_fake) <- col_names
df_fake[1, ] <- rep(NA, 10)

n_subj <- length(d$subj_id)

subj_list <- list()

for(i in 1:n_subj) {
  
  if (!is.na(subj_file_dat[i])) {
    one_subj <- read.table(
      here("data", "raw", "task_switching", GROUP, subj_file_dat[i])
    )
    one_subj$subj_id <- d$subj_id[i]
    one_subj$trial <- 1:120
    subj_list[[i]] <- one_subj
  } else {
    one_subj <- df_fake
    one_subj$subj_id <- d$subj_id[i]
    subj_list[[i]] <- one_subj
  }
}


# Create data.frame -------------------------------------------------------

df <- do.call(rbind.data.frame, subj_list)
df$V10 <- NULL

df1 <- df %>% 
  dplyr::rename(
    stim_category = V1, 
    stim_position = V2, # 1 = upper left, 2 = upper right, 3 = low right, 4 low left
    task_type = V3, # 1 = food, 2 = plant
    food_img = V4, # 1:4 iper-caloric food, 4:8 ipo-caloric food
    plant_img = V5, # 1:4 leaves, 4:8 flower
    block_type = V6, # 1, 2 single task, 0 mixed
    is_task_switch = V7, # 1 yes, 0 no
    is_correct = V8, # 1 correct, 2 error, 3 too slow
    rt = V9
  )


# Save data in RDS files --------------------------------------------------


if (GROUP == "patients") {
  saveRDS(df1, here("data", "interim", "task_switching", "patients_task_switching_data.rds"))
} else {
  saveRDS(df1, here("data", "interim", "task_switching", "controls_task_switching_data.rds"))
}


foo <- df1 %>%
  group_by(subj_id) %>%
  summarise(
    n = n_distinct(trial)
  )
foo %>% as.data.frame()

#---- eof 


