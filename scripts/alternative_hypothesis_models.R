# Load necessary libraries
library(mgcv)        # For GAM modeling
library(glmnet)      # For Lasso feature selection
library(ggplot2)     # For plotting
library(tidyverse)   # For data wrangling
library(here)
library(lubridate)   # For handling date and 
library(MASS)
library(gratia)
library(patchwork)

  # Step 1: Data Preparation
  physical_forcings_biology_table <- read.csv(here("data/MB_Wildfire_Obs/processed_data/joined_physical_drivers/physical_forcings_biology_table_alltime.csv"))
  data <- physical_forcings_biology_table
  data$datetime <- dmy(data$datetime)  # If date is day-month-year
  data$yearmonth <- year(data$datetime) * 100 + month(data$datetime)
  
  # Define predictor and response variables
  predictor_vars <- c('sst', 'BEUTI', 'wspd_along','wspd_across', 'surface_downwelling_photosynthetic_photon_flux_in_air', 
                      'pm2_5_sc_pm2_5_daily_slv_tt', 'Nitrate', 'Phosphate', 'Silicate', 'Ammonium', 'outflow')
  response_vars <- c('Asterionellopsis', 'Centric', 'Hemiaulus', 'Leptocylindrus', 'Skeletonema', 'Thalassionema', 'Thalassiosira','total_cells')
  
  # Clean data: remove rows with missing values and filter values as needed
  data_clean <- na.omit(data[, c(predictor_vars, response_vars)])
  data_clean <- filter(data_clean, surface_downwelling_photosynthetic_photon_flux_in_air >= 0)
  data_clean$pm2_5_sc_pm2_5_daily_slv_tt <- ifelse(data_clean$pm2_5_sc_pm2_5_daily_slv_tt < 0, 0, data_clean$pm2_5_sc_pm2_5_daily_slv_tt)
  
  # Step 2: Transformation and Scaling
  skewed_predictors <- c('Ammonium', 'Nitrate', 'outflow', 'Phosphate', 'pm2_5_sc_pm2_5_daily_slv_tt', 'Silicate')
  for (var in skewed_predictors) data_clean[[var]] <- log(data_clean[[var]] + 1)
  
  # Log-transform responses
  for (var in response_vars) data_clean[[var]] <- log(data_clean[[var]] + 1)
  
  # Center and scale predictors
  for (var in predictor_vars) data_clean[[var]] <- scale(data_clean[[var]], center = TRUE, scale = TRUE)
  
  # Step 3: Plot distributions after log transformations
  data_long <- data_clean %>% 
    pivot_longer(cols = c(predictor_vars, response_vars), names_to = "variable", values_to = "value")
  
  ggplot(data_long, aes(x = value)) +
    geom_histogram(bins = 30, fill = "blue", color = "white") +
    facet_wrap(~ variable, scales = "free") +
    labs(title = "Distributions of Predictors and Responses After Log Transformation", 
         x = "Value", y = "Count") +
    theme_minimal()
  
  # Step 4: Feature Selection and GAM Fitting
  # Initialize storage for models, feature importance, and AIC tracking
  models <- list()
  best_model_info <- data.frame(ResponseTaxa = character(), AIC = numeric(), ModelFormula = character(), stringsAsFactors = FALSE)
  feature_importance_table <- data.frame()
  
  # Loop through each response variable to build models
  for (response_var in response_vars) {
    # Define initial formula with all predictors and select = TRUE
    initial_formula <- as.formula(paste(response_var, "~", paste(paste0("s(", predictor_vars, ", k = 5)"), collapse = " + ")))
    initial_model <- gam(initial_formula, data = data_clean, select = TRUE, method = "REML")
    
    # Extract summary and identify significant predictors
    initial_summary <- summary(initial_model)
    predictor_summ <- initial_summary$s.table
    significant_predictors <- rownames(predictor_summ)[predictor_summ[, "p-value"] < 0.05]
    
    # Remove 's()' for final refined formula by extracting only variable names from significant predictors
    refined_predictors <- gsub("s\\((.*)\\)", "\\1", significant_predictors)  # Extract variable name from s() terms
    
    # If no predictors are significant, fall back to the initial model with all predictors
    if (length(refined_predictors) == 0) {
      refined_predictors <- predictor_vars  # Fallback to use all predictors if no significant predictors found
    }
    
    # Create refined formula without s() around all predictors
    refined_formula <- as.formula(paste(response_var, "~", paste(paste0("s(", refined_predictors, ", k = 5)"), collapse = " + ")))
    refined_model <- gam(refined_formula, data = data_clean, method = "REML")
    
    # Calculate AIC for refined model
    refined_model_aic <- AIC(refined_model)
    
    # Store the model and AIC if it is the best model for this response variable
    models[[response_var]] <- refined_model
    best_model_info <- rbind(best_model_info, data.frame(ResponseTaxa = response_var, AIC = refined_model_aic, ModelFormula = as.character(refined_formula)))
    
    # Extract feature importance for significant predictors in the final model
    refined_summary <- summary(refined_model)
    predictor_summ <- refined_summary$s.table
    stat_col <- if ("F" %in% colnames(predictor_summ)) "F"
    
    for (j in 1:nrow(predictor_summ)) {
      predictor_name <- rownames(predictor_summ)[j]
      predictor_stat_value <- predictor_summ[j, stat_col]
      predictor_p_value <- predictor_summ[j, "p-value"]
      
      # Determine direction based on the predictor effect in the model
      predictor_effect <- predict(refined_model, type = "terms")[, predictor_name]
      direction <- ifelse(mean(predictor_effect) > 0, "Positive", "Negative")
      
      feature_importance_table <- rbind(feature_importance_table,
                                        data.frame(ResponseTaxa = response_var,
                                                   Predictor = predictor_name,
                                                   Fvalue = predictor_stat_value,
                                                   p_value = predictor_p_value))
    }
  }
  
  # Save the best model information to CSV
  write.csv(best_model_info, here("C:/Users/Dante Capone/OneDrive/Desktop/Scripps_PhD/Wildfire_Obs/MB_Wildfire_Obs/tables/GAM_best_model_info.csv"), row.names = FALSE)
  
  # Save the feature importance table to CSV
  write.csv(feature_importance_table, here("C:/Users/Dante Capone/OneDrive/Desktop/Scripps_PhD/Wildfire_Obs/MB_Wildfire_Obs/tables/GAM_feature_importance_statistics.csv"), row.names = FALSE)
  
  
  # Optional: Print summaries of each refined model
  for (response_var in response_vars) {
    cat("\nRefined Model for", response_var, ":\n")
    print(summary(models[[response_var]]))
  }
  
  
  # Step 5: Partial Dependence Plots for GAMs
  output_dir <- "C:/Users/Dante Capone/OneDrive/Desktop/Scripps_PhD/Wildfire_Obs/MB_Wildfire_Obs/figures/physical_drivers/GAM_partial_dependence/R_GAMs/"
  font_size <- 10
  line_width <- 2  # Thicker line width for visibility
  
  # Define color scheme
  taxa_colors <- c(
    rgb(0.8, 0.8, 0.8), #Grey
    rgb(0.3, 0.7, 0.6),   # Soft blue-green
    rgb(0.2, 0.6, 0.8),  # Light blue
    rgb(0.9, 0.5, 0.1),  # Orange
    rgb(0.4, 0.7, 0.2),  # Light green
    rgb(0.8, 0.2, 0.4),  # Pink
    rgb(0.2, 0.4, 0.8),  # Dark blue
    rgb(0.7, 0.2, 0.9)  # Purple
    
  )
  
  predictor_var_titles <- c('Temperature', 'BEUTI', 'Alongshore Windspeed','Across-shore Windspeed', 'PAR', 
                      'PM 2.5', 'Nitrate', 'Phosphate', 'Silicate', 'Ammonium', 'Outflow')
  
  
  
  predictor_var_labels <- c(
    'Temperature [°C]', 
    'BEUTI', 
    'Along-shore windspeed [m s^-1]',
    'Across-shore windspeed [m s^-1]',
    'PAR [μmol photons m⁻² s⁻¹]', 
    'PM 2.5 [μg/m³]', 
    'Nitrate [µmol L⁻¹]', 
    'Phosphate [µmol L⁻¹]', 
    'Silicate [µmol L⁻¹]', 
    'Ammonium [µmol L⁻¹]', 
    'Outflow [cfs]'
  )
  
  
  output_dir <- "C:/Users/Dante Capone/OneDrive/Desktop/Scripps_PhD/Wildfire_Obs/MB_Wildfire_Obs/figures/physical_drivers/GAM_partial_dependence/R_GAMs/"
  
  for (response_var in response_vars) {
    gam_model <- models[[response_var]]
    
    retained_predictors <- purrr::map_chr(gam_model$smooth, ~ .x$term)
    plot_list <- list()
    
    for (j in seq_along(retained_predictors)) {
      
      pred <- retained_predictors[j]
      select_term <- j
      
      plt <- draw(gam_model, select = select_term, residuals = TRUE) +
        geom_point(data = data_clean, aes_string(x = pred, y = response_var),
                   alpha = 0.3, size = 1.2, color = taxa_colors[j]) +
        labs(
          x = predictor_var_labels[match(pred, predictor_vars)],
          y = paste("Partial Effect on log(", response_var, sep = ""),
          title = predictor_var_labels[match(pred, predictor_vars)]
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 12)
        )
      
      plot_list[[j]] <- plt
    }
    
    combined_plot <- wrap_plots(plot_list, ncol = 2) + plot_annotation(title = paste("GAM Effects for", response_var))
    
    ggsave(filename = paste0(output_dir, "GAM_combined_", response_var, ".png"),
           plot = combined_plot, width = 12, height = 24, dpi = 400)
  }
  
  
  


# 
# 
# # Random Forest -----------------------------------------------------------
# 
# library(randomForest) # For Random Forest regression
# 
# data_clean <- na.omit(data[, c(predictor_vars, response_vars)])
# 
# # Alternatively, Impute missing values with mean (example for predictors)
# for (var in c(predictor_vars, response_vars)) {
#   if (any(is.na(data_clean[[var]]))) {
#     data_clean[[var]][is.na(data_clean[[var]])] <- mean(data_clean[[var]], na.rm = TRUE)
#   }
# }
# 
# # Check for remaining missing values
# if (any(is.na(data_clean))) {
#   stop("Dataset still contains missing values after preprocessing.")
# }
# 
# # Step 2: Transformation and Scaling
# skewed_predictors <- c('Ammonium', 'Nitrate', 'outflow', 'Phosphate', 'pm2_5_sc_pm2_5_daily_slv_tt', 'Silicate')
# for (var in skewed_predictors) data_clean[[var]] <- log(data_clean[[var]] + 1)
# 
# # Log-transform responses
# for (var in response_vars) data_clean[[var]] <- log(data_clean[[var]] + 1)
# 
# # Train and scale predictors
# for (var in predictor_vars) data_clean[[var]] <- scale(data_clean[[var]], center = TRUE, scale = TRUE)
# 
# # Step 3: Random Forest Regression
# models <- list()
# feature_importance_table <- data.frame()
# 
# # Train a Random Forest model for each response variable
# for (response_var in response_vars) {
#   formula <- as.formula(paste(response_var, "~", paste(predictor_vars, collapse = " + ")))
#   rf_model <- randomForest(formula, data = data_clean, importance = TRUE, ntree = 500)
#   models[[response_var]] <- rf_model
#   
#   # Extract feature importance
#   importance_metrics <- importance(rf_model)
#   importance_df <- data.frame(
#     ResponseTaxa = response_var,
#     Predictor = rownames(importance_metrics),
#     Importance = importance_metrics[, "IncNodePurity"]
#   )
#   feature_importance_table <- rbind(feature_importance_table, importance_df)
# }
# 
# # Save the feature importance table to CSV
# # write.csv(feature_importance_table, here("C:/Users/Dante Capone/OneDrive/Desktop/Scripps_PhD/Wildfire_Obs/MB_Wildfire_Obs/tables/RF_feature_importance_statistics.csv"), row.names = FALSE)
# 
# # Step 4: Partial Dependence Plots for Random Forest
# output_dir <- "C:/Users/Dante Capone/OneDrive/Desktop/Scripps_PhD/Wildfire_Obs/MB_Wildfire_Obs/figures/physical_drivers/RFR/"
# 
# for (response_var in response_vars) {
#   if (!is.null(models[[response_var]])) {
#     rf_model <- models[[response_var]]
#     
#     for (predictor in predictor_vars) {
#       pdp_plot <- partial(rf_model, pred.var = predictor, train = data_clean)
#       
#       # Save the plot
#       g <- ggplot(pdp_plot, aes_string(x = predictor, y = "yhat")) +
#         geom_line(color = "blue", size = 1.2) +
#         labs(
#           title = paste("Partial Dependence of", response_var, "on", predictor),
#           x = predictor,
#           y = "Partial Dependence"
#         ) +
#         theme_minimal()
#       
#       ggsave(filename = paste0(output_dir, response_var, "_", predictor, "_pdp.png"), plot = g, width = 7, height = 5)
#     }
#   }
# }
