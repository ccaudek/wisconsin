library(dplyr)
library(tidyr)


file <- file.path("src", "wcst", "models", "stan", "steinke_7rev.stan") 

params_mod7 <-  c("mu_MB_Arew", "mu_MB_Apun", "mu_MB_gamma", "mu_MF_Arew",
                  "mu_MF_Apun", "mu_MF_gamma","mu_temp","mu_w", "MB_Arew",
                  "MB_Apun", "MB_gamma","MF_Arew", "MF_Apun", "MF_gamma","temp",
                  "w","log_lik","y_pred")

# compile model
mod <- cmdstan_model(file)

# Variational inference ---
fit_vi <- mod$variational(
  data = stan_data,
  seed = 1234
)

draws <- fit_vi$draws(variables = params_mod7, format = "data.frame")

# Trasforma i dati in un formato lungo per facilitare il calcolo delle medie per soggetto
draws_long <- draws %>%
  pivot_longer(cols = !starts_with("mu_"), names_to = "parameter", values_to = "value")


# Estrai l'indice del soggetto da ogni parametro e calcola la media
draws_long <- draws_long %>%
  mutate(subject = str_extract(parameter, "\\[\\d+\\]") %>% str_remove_all("\\[|\\]"),
         parameter = str_remove(parameter, "\\[\\d+\\]")) %>%
  group_by(subject, parameter) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd_value = sd(value, na.rm = TRUE),
            lower_95 = quantile(value, 0.025, na.rm = TRUE),
            upper_95 = quantile(value, 0.975, na.rm = TRUE),
            .groups = 'drop')

# View the result
draws_long %>% tbl_df %>% print(n=200)


# Supponiamo di aver già effettuato la stima con VI:
# fit_vi <- mod$variational(data = stan_data, seed = 1234)

summ <- fit_vi$summary()

# Supponiamo di conoscere il numero di soggetti N usato nei dati
N <- length(stan_data$Tsubj) # O un altro modo per sapere N

# Filtriamo solo i parametri indicizzati
subj_params <- summ %>%
  dplyr::filter(str_detect(variable, "\\[")) %>%
  # Escludiamo parametri che sappiamo non essere a livello soggetto:
  dplyr::filter(!str_detect(variable, "^mu_p"),
         !str_detect(variable, "^sigma"),
         !str_detect(variable, "^log_lik"),
         !str_detect(variable, "^y_pred"))

# Estrazione del nome parametro e indice soggetto
# Questo presuppone che i parametri soggetto abbiano una sola dimensione di indicizzazione,
# del tipo param[i]
subj_params <- subj_params %>%
  mutate(
    param = sub("\\[.*", "", variable),
    subj_id = as.integer(str_extract(variable, "(?<=\\[)\\d+(?=\\])"))
  )

# Ora verifichiamo che ogni (subj_id, param) appaia una sola volta.
# Se per qualche motivo c'è più di una riga per ciascun soggetto e parametro, 
# potremmo dover usare summarise per ottenere un'unica media.
# Idealmente, ogni (subj_id, param) dovrebbe essere unico.
# Controlliamo:
dups <- subj_params %>%
  group_by(subj_id, param) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)

if(nrow(dups) > 0){
  # Se qui entri significa che ci sono parametri con più righe per lo stesso soggetto.
  # In tal caso, potresti dover aggregare (ad es. prendere la media) prima di fare pivot.
  subj_params <- subj_params %>%
    group_by(subj_id, param) %>%
    summarise(mean = mean(mean), .groups = "drop")
} else {
  # Altrimenti abbiamo già una riga per soggetto e parametro
  subj_params <- subj_params %>%
    dplyr::select(subj_id, param, mean)
}

# Creiamo la tabella larga
wide_df <- subj_params %>%
  pivot_wider(names_from = param, values_from = mean) %>%
  arrange(subj_id)

wide_df$is_patient <- c(rep(1, 45), rep(0, 46))

# wide_df ora dovrebbe avere una colonna 'subj_id' e le altre colonne come i parametri,
# tutte numeric non-list

# Infine salviamo il risultato
rio::export(
  wide_df, 
  here::here("src", "stan_auc", "models_params", "wcst_steinke_rev_params.csv")
)


t.test(MF_gamma_pr ~ is_patient, wide_df)

