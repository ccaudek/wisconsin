

get_one_subj_data_for_stan <- function(df, subj_name) {
  # Filter data for the given subject
  df_subj <- dplyr::filter(df, subj_name == !!subj_name)
  
  # Process the data for the filtered subject (df_subj)
  df_subj <- df_subj %>%
    dplyr::mutate(
      rule_choice = dplyr::case_when(
        name_of_task == "color"  ~ 1,
        name_of_task == "shape"  ~ 2,
        name_of_task == "number" ~ 3,
        TRUE  ~ NA_integer_
      ),
      resp_choice = correct_card,
      resp_color = dplyr::case_when(
        card_color == "red"    ~ 1,
        card_color == "green"  ~ 2,
        card_color == "blue"   ~ 3,
        card_color == "yellow" ~ 4,
        TRUE  ~ NA_integer_
      ),
      resp_shape = dplyr::case_when(
        card_shape == "circle"    ~ 1,
        card_shape == "triangle"  ~ 2,
        card_shape == "cross"     ~ 3,
        card_shape == "star"      ~ 4,
        TRUE  ~ NA_integer_
      ),
      resp_number = as.integer(card_number),
      rew = ifelse(is_correct == 1, 1, 0),
      los = ifelse(is_error == 1, -1, 0)
    ) %>%
    dplyr::select(rew, los, rule_choice, resp_choice, resp_color, resp_shape, resp_number)
  
  return(df_subj)
}


compile_data_for_stan <- function(df) {
  df <- df %>%
    dplyr::mutate(group = dplyr::case_when(
      group == "an" ~ 1,
      group == "hc" ~ 2,
      group == "ri" ~ 3,
      TRUE ~ NA_integer_
    ))
  
  subject_names <- unique(df$subj_name)
  N <- length(subject_names)
  T_max <- 60  # As each subject has 60 trials
  
  # Initialize variables
  Tsubj <- rep(T_max, N)  # Each subject has 60 trials
  group <- integer(N)
  holdout <- rep(0, N)  # Assuming no holdout
  
  # Initialize matrices for each response type
  rew_matrix <- matrix(NA_integer_, N, T_max)
  los_matrix <- matrix(NA_integer_, N, T_max)
  rule_choice_matrix <- matrix(NA_integer_, N, T_max)
  resp_choice_matrix <- matrix(NA_integer_, N, T_max)
  resp_color_matrix <- matrix(NA_integer_, N, T_max)
  resp_shape_matrix <- matrix(NA_integer_, N, T_max)
  resp_number_matrix <- matrix(NA_integer_, N, T_max)
  
  for (i in seq_along(subject_names)) {
    subj_data <- get_one_subj_data_for_stan(df, subject_names[i])
    group[i] <- unique(df %>% dplyr::filter(subj_name == subject_names[i]) %>% dplyr::pull(group))[1]
    
    # Fill matrices
    rew_matrix[i, 1:T_max] <- as.integer(subj_data$rew)
    los_matrix[i, 1:T_max] <- as.integer(subj_data$los)
    rule_choice_matrix[i, 1:T_max] <- as.integer(subj_data$rule_choice)
    resp_choice_matrix[i, 1:T_max] <- as.integer(subj_data$resp_choice)
    resp_color_matrix[i, 1:T_max] <- as.integer(subj_data$resp_color)
    resp_shape_matrix[i, 1:T_max] <- as.integer(subj_data$resp_shape)
    resp_number_matrix[i, 1:T_max] <- as.integer(subj_data$resp_number)
  }
  
  # Return the compiled data list
  return(list(
    N = N, 
    T = T_max, 
    Tsubj = Tsubj, 
    holdout = holdout, 
    group = group,
    rew = rew_matrix, 
    los = los_matrix, 
    rule_choice = rule_choice_matrix, 
    resp_choice = resp_choice_matrix,
    resp_color = resp_color_matrix, 
    resp_shape = resp_shape_matrix, 
    resp_number = resp_number_matrix
  ))
}



# Load necessary libraries
library("here")
library("tidyverse")
library("stringi")
library("rio")

# Source functions
source(here::here("src", "functions", "funs_wcst.R"))
source(here::here("src", "functions", "funs_input_for_stan_wcst.R"))

# Function to generate RDS raw data for patients or controls
# Input: group - a string, either "patients" or "controls"
# Output: Generates a CSV file in the specified directory if it does not already exist
# Requirements: The function requires the 'here', 'tidyverse', 'stringi', and 'rio' packages to be installed and loaded.
#               It also requires the source functions to be available in the specified paths.
#
# Example usage
# generate_csv("controls")
# generate_csv("patients")

generate_csv <- function(group) {
  
  # Validate group argument
  if (!(group %in% c("patients", "controls"))) {
    stop("Invalid group. Choose either 'patients' or 'controls'.")
  }
  
  # Define file path for the output CSV
  file_path <- here("src", "wcst_model_selection", "data", paste0("raw_data_", group, ".csv"))
  
  # Check if the file already exists
  if (file.exists(file_path)) {
    message("File already exists: ", file_path)
    return(NULL)
  }
  
  # Define directory and file pattern based on the group argument
  dir <- here("data", "raw", group)
  pattern <- if (group == "patients") "wcst_pazienti" else "wcst_eds1"
  
  # List all files in the directory that match the pattern
  file_names <- as.character(list.files(path = dir, pattern = pattern))
  n_files <- length(file_names)
  
  # Initialize list to store data from each file
  d_list <- list()
  
  # Loop over each file, read the data, and process it
  for (i in 1:n_files) {
    d <- read.table(here("data", "raw", group, file_names[i]), header = FALSE)
    d$subj_name <- file_names[i]
    d$block <- rep(1:6, each = 10)
    
    # Rename columns for clarity
    d <- d %>%
      rename(
        card_shown = V1,
        correct_card = V2,
        card_chosen_if_perseveration = V3,
        trial_in_a_sequence = V4,
        name_of_task = V5,
        card_shape = V6,
        card_number_of_symbols = V7,
        card_color = V8,
        rt = V9,
        is_correct = V10, 
        chosen_card = V11,
        is_error = V12, 
        is_perseverative_error = V13, 
        is_non_perseverative_error = V14
      )
    
    # Append the processed data to the list
    d_list[[i]] <- d
  }
  
  # Combine all data frames in the list into a single data frame
  df0 <- do.call(rbind.data.frame, d_list)
  
  # Calculate the error rate for each subject and identify subjects with error rate > 0.5
  bysubj_acc <- df0 %>%
    group_by(subj_name) %>%
    summarise(error_rate = mean(is_error)) %>%
    as.data.frame()
  
  bad_subj_df <- bysubj_acc %>%
    filter(error_rate > 0.5)
  
  bad_subjects <- bad_subj_df$subj_name
  
  # Remove data from subjects with high error rates
  df1 <- df0[!(df0$subj_name %in% bad_subjects), ]
  
  # Remove rows with NA in chosen_card column and rename columns as needed
  df <- df1[!is.na(df1$chosen_card), ] %>%
    rename(card_number = card_number_of_symbols)
  
  # Export the cleaned data to a CSV file
  rio::export(df, file_path)
  message("File created: ", file_path)
}





process_and_prepare_stan_data <- function(DEBUGGING = FALSE) {
  # This function imports, processes, and prepares data for further analysis using Stan models.
  # 1. Import Libraries: The function uses libraries like `rio` for data import, `here` for constructing file paths, and `dplyr` for data manipulation.
  # 2. Read Data: It reads raw data from CSV files for patients and controls.
  # 3. Combine Data: The data from patients and controls are combined into a single dataframe.
  # 4. Rename Columns: It renames a specific column for clarity.
  # 5. Add User IDs: Correspondence tables for patients and controls are generated and merged into the main data.
  # 6. Final Data Processing: The final processed data is imported and subjects are filtered based on specific groups.
  # 7. Filter and Select Data: It filters subjects who completed both tasks and creates a 'group' column.
  # 8. Prepare Data for Stan: Finally, it prepares the data for use in Stan models.
  
  # Requirements:
  # CSV files: `raw_data_patients.csv`, `raw_data_controls.csv`, and `raw_data_project.csv` located in specific directories.
  
  # Load necessary libraries
  library(rio)
  library(here)
  library(dplyr)
  
  # Read raw data
  raw_data_patients_df <- rio::import(here::here("src", "wcst_model_selection", "data", "raw_data_patients.csv"))
  raw_data_controls_df <- rio::import(here::here("src", "wcst_model_selection", "data", "raw_data_controls.csv"))
  
  # Combine data and rename columns
  raw_data_df <- rbind(raw_data_patients_df, raw_data_controls_df) %>%
    dplyr::rename(code_psytoolkit = subj_name)
  
  # Generate correspondence tables
  patients_wcst_look_up_tbl <- gen_correspondence_table_codes("patients")
  controls_wcst_look_up_tbl <- gen_correspondence_table_codes("controls")
  
  # Combine correspondence tables
  look_up_tbl <- rbind(patients_wcst_look_up_tbl, controls_wcst_look_up_tbl)
  
  # Merge correspondence tables with raw data
  raw_data_df <- left_join(raw_data_df, look_up_tbl, by = "code_psytoolkit")
  
  # Import final data
  final_data <- rio::import(here::here("data", "processed", "wcst", "raw_data_project.csv"))
  
  # Group subjects by their group and select unique subjects for each group
  code_psytoolkit_selected_an <- unique(final_data[final_data$group == "an", ]$subj_name)
  code_psytoolkit_selected_hc <- unique(final_data[final_data$group == "hc", ]$subj_name)
  code_psytoolkit_selected_ri <- unique(final_data[final_data$group == "ri", ]$subj_name)
  
  # If DEBUGGING is TRUE, take only two subjects from each group
  if (DEBUGGING) {
    code_psytoolkit_selected_an <- head(code_psytoolkit_selected_an, 2)
    code_psytoolkit_selected_hc <- head(code_psytoolkit_selected_hc, 2)
    code_psytoolkit_selected_ri <- head(code_psytoolkit_selected_ri, 2)
  }
  
  # Vector of subjects who have completed both tasks.
  # Select only participants of the AN and HC groups.
  keep_psytoolkit_codes <- c(
    code_psytoolkit_selected_an, code_psytoolkit_selected_hc
    # code_psytoolkit_selected_ri
  )
  
  # Filter data to keep only the subjects who completed both tasks
  raw_data_two_groups_df <- 
    raw_data_df[raw_data_df$code_psytoolkit %in% keep_psytoolkit_codes, ]
  
  # Create the 'group' column
  raw_data_two_groups_df <- raw_data_two_groups_df %>%
    mutate(group = case_when(
      code_psytoolkit %in% code_psytoolkit_selected_an ~ "an",
      code_psytoolkit %in% code_psytoolkit_selected_hc ~ "hc",
      # code_psytoolkit %in% code_psytoolkit_selected_ri ~ "ri",
      TRUE ~ NA_character_
    ))
  
  # Remove unwanted columns
  raw_df <- raw_data_two_groups_df %>%
    dplyr::select(!starts_with("V"))
  
  # Generate the input for the Stan models
  stan_data <- compile_data_for_stan(raw_df)
  
  return(stan_data)
}


