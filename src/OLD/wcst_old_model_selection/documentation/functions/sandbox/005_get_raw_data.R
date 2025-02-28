library(here)
library(tidyverse)


d <- rio::import(
  here::here(
    "data", "processed", "wcst", "raw_data_project.csv"
  )
) 

length(unique(d$subj_name))


d$group <- factor(d$group)
summary(d$group)

# Remove At-Risk participants.
# Remove unnecessary variables.
df <- d |> 
  dplyr::select(-starts_with("V")) |> 
  dplyr::filter(group != "ri")

df |> 
  group_by(group) |> 
  summarize(
    n = n_distinct(subj_name)
  )
#   group     n
# 1 an       45
# 2 hc       46

rio::export(
  df,
  here::here(
    "data", "processed", "wcst", "raw_2grps_4classification.csv"
  )
)

# eof ----

