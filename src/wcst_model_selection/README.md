The scripts in this folder are used for model selection with the data of the WCST. 
The models are the 5 proposed by Seinke et al. (2020) -- the AUC (attention based) models are  not implemented.

In the script `01_models_comparisons.R`, the `generate_csv()` function is used for data wrangling and to generate the CVS files with the raw data of the patients (46) and controls (45).

By using as input the CSV files created by `generate_csv()`, the function `process_and_prepare_stan_data()` generates the list used as input for the Stan models created by Seinke et al. (2020). The same data list is used as input for all Steinke's models.

The `stan_data` list contains the data of 91 participants (46 patients and 45 controls) who completed the WCST task.

In the `Chose model` section of the script, by setting the `model_flag` it is possible to choose one of the five models that are considered. This operation must be repeated five times, for all models.

Once chosen the model and the list of its parameters, the model is compiled. 

Sampling is done with the `pathfinder` variational inference method. However, `sample()` should be used for the model `07_PRL_weighting.stan` in order to obtain a better estimate of the model's parameters, for each subject.

The purpose of the script is to estimate WAIC for each of the five models. This is achieved by running the script five times. The results are saved at the end of the script.

The estimates of the model's parameters, obtained with `sample()`, for the model `07_PRL_weighting.stan` have been saved in `src/stan_auc/models_parameters/wcst_steinke_params.csv`

