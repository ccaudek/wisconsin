
# gen_correspondence_table_codes() ----------------------------------------


#' @description 
#' by using the information in the Excel file, generates a correspondence table
#' which associates the codes used by psytoolkit with the codes used in the
#' questionnaires
#' @return data.frame.
gen_correspondence_table_codes <- function(GROUP) {
  
  d <- read_excel_code(GROUP)
  d <- gen_subj_name(d)
  
  d_clean <- d %>% 
    dplyr::rename("subj_name" = "subj_id") %>% 
    dplyr::select(subj_name, code_psytoolkit)
  
  d_clean2 <- d_clean[!is.na(d_clean$code_psytoolkit), ] 
  # d_food_clean2$code_psytoolkit
  d_clean2
}


# read_excel_code() -------------------------------------------------------


#' @description 
#' read Excel file.
#' @return data.frame.
#' @param 
#' GROUP = "patients", EXCEL_FILE = "misc_food".
read_excel_code <- function(GROUP) {
  d <- readxl::read_excel(
    here("data", "raw", GROUP, "data.xlsx")
  )
  d
}


# gen_subj_name() ---------------------------------------------------------


#' @description 
#' generate subject code
#' @return data.frame.
#' @param data.frame.
gen_subj_name <- function(d) {
  
  library("stringi")
  
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
  
  d$code_psytoolkit <- d$`esperimento:1`
  
  d
}



