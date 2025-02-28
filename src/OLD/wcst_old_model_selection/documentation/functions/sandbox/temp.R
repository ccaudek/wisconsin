

# Read raw data
raw_data_patients_df <- rio::import(
  here::here("src", "wcst_model_selection", "data", 
             "raw_data_patients.csv")
)

length(unique(raw_data_patients_df$subj_name))
# [1] 60

raw_data_controls_df <- rio::import(
  here::here("src", "wcst_model_selection", "data", 
             "raw_data_controls.csv")
)

# join
raw_data_df <- rbind(raw_data_patients_df, raw_data_controls_df) |> 
  dplyr::rename(
    code_psytoolkit = subj_name
  )
length(unique(raw_data_df$code_psytoolkit))
# [1] 338

# Add user_id
patients_wcst_look_up_tbl <- gen_correspondence_table_codes("patients")
controls_wcst_look_up_tbl <- gen_correspondence_table_codes("controls")
look_up_tbl <- rbind(patients_wcst_look_up_tbl, controls_wcst_look_up_tbl)
dim(look_up_tbl)
# [1] 346   2

raw_data_df <- left_join(raw_data_df, look_up_tbl, by = "code_psytoolkit")
length(unique(raw_data_df$subj_name))
# [1] 335

final_data <- rio::import(
  here::here("data", "processed", "wcst", "raw_data_project.csv")
)

final_data |> 
  group_by(group) |> 
  summarize(
    n = n_distinct(subj_name)
  )
#   group     n
# 1 an       45
# 2 hc       46
# 3 ri       24

code_psytoolkit_selected_an <- 
  final_data[final_data$group == "an", ]$subj_name |> 
  unique()

code_psytoolkit_selected_hc <- 
  final_data[final_data$group == "hc", ]$subj_name |> 
  unique()

code_psytoolkit_selected_ri <- 
  final_data[final_data$group == "ri", ]$subj_name |> 
  unique()

# Vector of code_psytoolkit codes for subjects who have completed 
# the WCST and the PRL task
keep_psytoolkit_codes <- c(
  code_psytoolkit_selected_an, code_psytoolkit_selected_hc, 
  code_psytoolkit_selected_ri
)

# Select only the subjects who have completed both tasks
raw_data_three_groups_df <- 
  raw_data_df[raw_data_df$code_psytoolkit %in% keep_psytoolkit_codes, ]

length(unique(raw_data_three_groups_df$code_psytoolkit))
# [1] 115

# Create the 'group' column
raw_data_three_groups_df <- raw_data_three_groups_df |> 
  mutate(group = case_when(
    code_psytoolkit %in% code_psytoolkit_selected_an ~ "an",
    code_psytoolkit %in% code_psytoolkit_selected_hc ~ "hc",
    code_psytoolkit %in% code_psytoolkit_selected_ri ~ "ri",
    TRUE ~ NA_character_
  ))

raw_data_three_groups_df |> 
  group_by(group) |> 
  summarise(
    n = n_distinct(subj_name)
  )
#  group     n
# 1 an       45
# 2 hc       46
# 3 ri       24

raw_df <- raw_data_three_groups_df |> 
  dplyr::select(
    !starts_with("V")
  )

# Generate the input for the Stan models
stan_data <- compile_data_for_stan(raw_df)
