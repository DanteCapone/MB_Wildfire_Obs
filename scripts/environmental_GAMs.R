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
library(ggtext)

  # Data Preparation
  physical_forcings_biology_table <- read.csv(here("data/joined_physical_drivers/physical_forcings_biology_table_alltime.csv"))
  data <- physical_forcings_biology_table
  data$datetime <- dmy(data$datetime)  # If date is day-month-year
  data$yearmonth <- year(data$datetime) * 100 + month(data$datetime)
  
  # Define predictor and response variables
  predictor_vars <- c('sst', 'BEUTI', 'wspd_along','wspd_across', 'surface_downwelling_photosynthetic_photon_flux_in_air', 
                      'pm2_5_sc_pm2_5_daily_slv_tt', 'Nitrate', 'Phosphate', 'Silicate', 'Ammonium', 'outflow')
  response_vars <- c('Asterionellopsis', 'Centric', 'Hemiaulus', 'Leptocylindrus', 'Skeletonema', 'Thalassionema', 'Thalassiosira','total_cells')
  
  # Clean data: remove rows with missing values and filter values as needed
  data_clean <- na.omit(data[, c(predictor_vars, response_vars)])
  data_clean <- filter(data_clean, surface_downwelling_photosynthetic_photon_flux_in_air >= 0,
                       Phosphate>= 0.0, Nitrate >= 0.0, Silicate >= 0, Ammonium >= 0,
                       Asterionellopsis >=0, Centric>=0, Hemiaulus>=0, Leptocylindrus>=0, Skeletonema>=0, Thalassionema>=0, Thalassiosira >=0)
  data_clean$pm2_5_sc_pm2_5_daily_slv_tt <- ifelse(data_clean$pm2_5_sc_pm2_5_daily_slv_tt < 0, 0, data_clean$pm2_5_sc_pm2_5_daily_slv_tt)
  
  #Transformation and Scaling
  skewed_predictors <- c('Ammonium', 'Nitrate', 'outflow', 'Phosphate', 'pm2_5_sc_pm2_5_daily_slv_tt', 'Silicate')
  for (var in skewed_predictors) data_clean[[var]] <- log(data_clean[[var]] + 1)
  # Log-transform responses
  for (var in response_vars) data_clean[[var]] <- log(data_clean[[var]] + 1)
  
  # Center and scale predictors
  for (var in predictor_vars) data_clean[[var]] <- scale(data_clean[[var]], center = TRUE, scale = TRUE)
  
  #Plot distributions after log transformations
  data_long <- data_clean %>% 
    pivot_longer(cols = c(predictor_vars, response_vars), names_to = "variable", values_to = "value")
  
  ggplot(data_long, aes(x = value)) +
    geom_histogram(bins = 30, fill = "blue", color = "white") +
    facet_wrap(~ variable, scales = "free") +
    labs(title = "Distributions of Predictors and Responses After Log Transformation", 
         x = "Value", y = "Count") +
    theme_minimal()
  
  # Feature Selection and GAM Fitting
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
  
  
  #Partial Dependence Plots for GAMs
  output_dir <- "C:/Users/Dante Capone/OneDrive/Desktop/Scripps_PhD/Wildfire_Obs/MB_Wildfire_Obs/figures/physical_drivers/GAM_partial_dependence/R_GAMs/"
  font_size <- 10
  line_width <- 2  # Thicker line width for visibility
  
  # Define color scheme
  taxa_colors <- c(
    rgb(0.2, 0.6, 0.8),  # Light blue
    rgb(0.9, 0.5, 0.1),  # Orange
    rgb(0.4, 0.7, 0.2),  # Light green
    rgb(0.8, 0.2, 0.4),  # Pink
    rgb(0.2, 0.4, 0.8),  # Dark blue
    rgb(0.7, 0.2, 0.9),  # Purple   
    rgb(0.3, 0.7, 0.6),   # Soft blue-green
    rgb(0.8, 0.8, 0.8) #Grey
    
    
  )
  
  predictor_var_titles <- c('Temperature', 'BEUTI', 'Along-shore\nWindspeed','Across-shore\nWindspeed', 'PAR', 
                      'PM 2.5', 'Nitrate', 'Phosphate', 'Silicate', 'Ammonium', 'Outflow')
  
  
  
  predictor_var_labels <- c(
    'Temperature [°C]', 
    'BEUTI', 
    'Along-shore\nwindspeed [m s⁻¹]',
    'Across-shore\nwindspeed [m s⁻¹]',
    'PAR [μmol photons m⁻² s⁻¹]', 
    'PM 2.5 [μg m⁻³]', 
    'Nitrate [µmol L⁻¹]', 
    'Phosphate [µmol L⁻¹]', 
    'Silicate [µmol L⁻¹]', 
    'Ammonium [µmol L⁻¹]', 
    'Outflow [cfs]'
  )
  
response_vars <- c('Asterionellopsis', 'Centric', 'Hemiaulus', 'Leptocylindrus', 'Skeletonema', 'Thalassionema', 'Thalassiosira',"Total cells")
  
response_var_labels <- ifelse(response_vars %in% c("Centric", "Total cells"),
                              response_vars,
                              paste0(response_vars))

response_var_titles <- ifelse(response_vars %in% c("Centric", "Total cells"),
                              response_vars,
                              paste0("*", response_vars, "*"))
  
  output_dir <- "C:/Users/Dante Capone/OneDrive/Desktop/Scripps_PhD/Wildfire_Obs/MB_Wildfire_Obs/figures/physical_drivers/GAM_partial_dependence/R_GAMs/"
  
  for (i in seq_along(response_vars)) {
    response_var <- response_vars[i]
    gam_model <- models[[response_var]]
    retained_predictors <- purrr::map_chr(gam_model$smooth, ~ .x$term)
    
    # Retrieve and order by F-stat
    preds_ordered <- feature_importance_table %>%
      filter(ResponseTaxa == response_var, Predictor %in% paste0("s(", retained_predictors, ")")) %>%
      arrange(desc(Fvalue)) %>%
      pull(Predictor) %>%
      str_remove_all("s\\(|\\)")
    
    plot_list <- list()
    n_preds <- length(preds_ordered)
    ncol_plots <- 3
    nrow_plots <- ceiling(n_preds / ncol_plots)
    
    for (j in seq_along(preds_ordered)) {
      pred <- preds_ordered[j]
      select_term <- which(retained_predictors == pred)
      
      f_stat <- feature_importance_table %>%
        filter(ResponseTaxa == response_var, Predictor == paste0("s(", pred, ")")) %>%
        pull(Fvalue) %>%
        round(2)
      
      if(length(f_stat) == 0) f_stat <- NA
      
      # Create a complete prediction grid
      pred_grid <- data_clean %>% 
        summarise(across(all_of(predictor_vars), mean, na.rm = TRUE)) %>% 
        slice(rep(1, 200))  # repeat row 200 times
      
      # Now vary only the predictor of interest across its full observed range
      pred_grid[[pred]] <- seq(
        min(data_clean[[pred]], na.rm = TRUE),
        max(data_clean[[pred]], na.rm = TRUE),
        length.out = 200
      )
      
      # Compute smooth estimates explicitly on this full range
      sm <- smooth_estimates(gam_model, data = pred_grid, select = pred, partial_match = TRUE) %>%
        mutate(
          .lower_ci = .estimate - 1.96 * .se,
          .upper_ci = .estimate + 1.96 * .se
        )
      # Calculate partial residuals explicitly
      partial_resid <- data_clean %>%
        mutate(partial_residual = residuals(gam_model, type = "working") + predict(gam_model, type = "terms")[, paste0("s(", pred, ")")])
      
      if ((j - 1) %% ncol_plots == 0) {
        if (response_var_labels[i] %in% c("Centric", "Total cells")) {
          y_label <- bquote("Partial Effect on" ~ ln(.(response_var_labels[i])))  # plain
        } else {
          y_label <- bquote("Partial Effect on" ~ ln(italic(.(response_var_labels[i]))))  # italic
        }
      } else {
        y_label <- NULL
      }
      
      plt <- sm %>%
        ggplot(aes(x = .data[[pred]])) +
        geom_point(data = partial_resid, aes_string(x = pred, y = "partial_residual"),
                   alpha = 0.5, size = 1.5, colour = taxa_colors[i]) +
        geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2) +
        geom_line(aes(y = .estimate), linewidth = 1.2, alpha = 0.8) +
        labs(
          x = predictor_var_labels[match(pred, predictor_vars)],
          y = y_label,
          title = predictor_var_titles[match(pred, predictor_vars)]
        )+
        theme_minimal(base_size = 8) +
        theme(
          plot.title = element_text(size = 8, face = "bold"),
          axis.title = element_text(size = 11),
          plot.margin = margin(5.5, 5.5, 5.5, 20)
        )
      
      plot_list[[j]] <- plt
    }
    
    combined_plot <- wrap_plots(plot_list, ncol = ncol_plots, guides = "collect") &
      theme(axis.title.y = element_text(size = 8))
    
    combined_plot <- combined_plot + plot_annotation(
      title = response_var_titles[i],
      tag_levels = 'a',
      theme = theme(
        plot.title = element_markdown(size = 10, face = "bold"),
        plot.tag = element_text(size = 10, face = "bold")
      )
    )
    
    # Dynamically scale height based on number of rows
    plot_height <- 2.5 * nrow_plots
    
    ggsave(filename = paste0(output_dir, "GAM_combined_", response_var, "_sized.png"),
           plot = combined_plot, width = 5, height = plot_height, dpi = 1200)
    
    ggsave(filename = paste0(output_dir, "GAM_combined_", response_var, "_sized.pdf"),
           plot = combined_plot, width = 5, height = plot_height, dpi = 1200)
  }

