


# This function will be used by 101_get_mle_mod7.R.
# get_one_subj_data_for_stan <- function(df) {
#   
#   ## ------------------------------------------------------------------
#   ## rule_choice
#   ## ------------------------------------------------------------------
#   
#   # rule_choice indicates the category (color = 1, shape = 2, 
#   # number = 3) which is rewarded in each trial.
#   
#   df <- df %>%
#     mutate(
#       rule_choice = case_when(
#         name_of_task == "color"  ~ 1,
#         name_of_task == "shape"  ~ 2,
#         name_of_task == "number" ~ 3,
#         TRUE  ~ 999
#       )
#     )
#   
#   ## ------------------------------------------------------------------
#   ## resp_choice
#   ## ------------------------------------------------------------------
#   
#   # resp_choice indicates the key which the subject must press in order 
#   # to provide a correct response in each trial. 
#   # The four possible response keys are:
#   # 1 -> one red circle
#   # 2 -> two green triangles
#   # 3 -> three blue crosses
#   # 4 -> four yellow stars
#   
#   # in the psytoolbox output, this seems to correspond to correct_card.
#   df$resp_choice <- df$correct_card
#   
#   ## ------------------------------------------------------------------
#   ## resp_color
#   ## ------------------------------------------------------------------
#   
#   # resp_color indicates the key that the subject would press if she
#   # would chose the color dimension.
#   
#   df <- df %>%
#     mutate(
#       resp_color = case_when(
#         card_color == "red"    ~ 1,
#         card_color == "green"  ~ 2,
#         card_color == "blue"   ~ 3,
#         card_color == "yellow" ~ 4,
#         TRUE  ~ 999
#       )
#     )
#   
#   ## ------------------------------------------------------------------
#   ## resp_shape
#   ## ------------------------------------------------------------------
#   
#   # resp_shape indicates the key that the subject would press if she
#   # would chose the shape dimension.
#   
#   df <- df %>%
#     mutate(
#       resp_shape = case_when(
#         card_shape == "circle"    ~ 1,
#         card_shape == "triangle"  ~ 2,
#         card_shape == "cross"     ~ 3,
#         card_shape == "star"      ~ 4,
#         TRUE  ~ 999
#       )
#     )
#   
#   ## ------------------------------------------------------------------
#   ## resp_number
#   ## ------------------------------------------------------------------
#   
#   # resp_number indicates the key that the subject would press if she
#   # would chose the number dimension.
#   
#   df <- df %>%
#     mutate(
#       resp_number = case_when(
#         card_number == 1  ~ 1,
#         card_number == 2  ~ 2,
#         card_number == 3  ~ 3,
#         card_number == 4  ~ 4,
#         TRUE  ~ 999
#       )
#     )
#   
#   ## ------------------------------------------------------------------
#   ## rew
#   ## ------------------------------------------------------------------
#   
#   # rew indicate whether a reward has been provided (1) or not (0)
#   
#   df$rew <- ifelse(df$is_correct == 1, 1, 0)
#   
#   
#   ## ------------------------------------------------------------------
#   ## los
#   ## ------------------------------------------------------------------
#   
#   # low indicate whether a punishment has been provided (-1) or not (0)
#   
#   df$los <- ifelse(df$is_error == 1, -1, 0)
#   
#   
#   ## ------------------------------------------------------------------
#   ## holdout
#   ## ------------------------------------------------------------------
#   
#   # holdout == 1 means that the subject is excluded from the analysis.
#   # Here, holdout is always equal to 0.
#   
#   df$holdout <- rep(0, nrow(df))
#   
#   
#   ## ------------------------------------------------------------------
#   ## N
#   ## ------------------------------------------------------------------
#   
#   # N is the number of subjects
#   
#   N <- length(unique(df$subj_name))
#   
#   
#   ## ------------------------------------------------------------------
#   ## T
#   ## ------------------------------------------------------------------
#   
#   # T is the maximum number of trials for each subject. In our case
#   # is always equal to 60.
#   
#   T <- 60
#   
#   
#   ## ------------------------------------------------------------------
#   ## Tsubj
#   ## ------------------------------------------------------------------
#   
#   # Tsubj is the effective number of trials for each subject. In our
#   # case is always equal to 60.
#   
#   Tsubj <- 60
#   
#   
#   ## ------------------------------------------------------------------
#   ## odd_response
#   ## ------------------------------------------------------------------
#   
#   # odd_response seems to be the keypress of the subject (the chosen
#   # category), with -1 indicating no response. Here, there are no -1.
#   
#   # In the psytoolbox output, this seems to correspond to chosen_card.
#   df$odd_response <- df$chosen_card
#   
#   
#   ## ------------------------------------------------------------------
#   # Rearange data for modeling.
#   ## ------------------------------------------------------------------
#   
#   Tsubj <- 60
#   
#   # each subject is in a single row.
#   rew <- matrix(
#     df$rew,
#     nrow = N,
#     byrow = TRUE
#   )
#   
#   los <- matrix(
#     df$los,
#     nrow = N,
#     byrow = TRUE
#   )
#   
#   odd_response <- matrix(
#     df$odd_response,
#     nrow = N,
#     byrow = TRUE
#   )
#   
#   rule_choice <- matrix(
#     df$rule_choice,
#     nrow = N,
#     byrow = TRUE
#   )
#   
#   resp_choice <- matrix(
#     df$resp_choice,
#     nrow = N,
#     byrow = TRUE
#   )
#   
#   resp_color <- matrix(
#     df$resp_color,
#     nrow = N,
#     byrow = TRUE
#   )
#   
#   resp_shape <- matrix(
#     df$resp_shape,
#     nrow = N,
#     byrow = TRUE
#   )
#   
#   resp_number <- matrix(
#     df$resp_number,
#     nrow = N,
#     byrow = TRUE
#   )
#   
#   holdout = rep(0, N) # No data is hold out
#   
#   # Wrap-up data
#   data <- list(
#     N      = N,
#     T      = T,
#     Tsubj  = as.vector(rep(60, N)),
#     rew    = rew,     
#     los    = los,
#     odd_response = odd_response,
#     rule_choice = rule_choice,
#     resp_choice = resp_choice,
#     resp_color  = resp_color,
#     resp_shape  = resp_shape,
#     resp_number = resp_number,
#     holdout     = rep(0, N) # No data is hold out
#   )
#   
#   data
# }

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







