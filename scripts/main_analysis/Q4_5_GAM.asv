%% GAM v5
%% GAM Analysis with Lasso Feature Selection, Partial Dependence, and Residual Diagnostics

% Step 1: Data Preparation
data = physical_forcings_biology_table;  % Load dataset
data.yearmonth = datenum(datetime(data.datetime.Year, data.datetime.Month, 1, 'TimeZone', 'UTC'));  % Add year-month variable

% Define predictor and response variables
predictor_vars = {'temperature', 'BEUTI', 'wspd_along', ...
                  'surface_downwelling_photosynthetic_photon_flux_in_air', 'pm2_5_sc_pm2_5_daily_slv_tt', ...
                  'Nitrate', 'Phosphate', 'Silicate', 'Ammonium', 'outflow', 'yearmonth'};
response_vars = {'Asterionellopsis', 'Centric', 'Skeletonema', 'Hemiaulus', 'Leptocylindrus', 'Thalassionema', 'Thalassiosira'};

% Clean data (remove rows with missing values)
data_clean = timetable2table(rmmissing(data(:, [predictor_vars, response_vars])));  
data_clean.datetime = [];  % Remove datetime column

% Log-transform the response variables (avoid log(0) by adding small constant)
for i = 1:length(response_vars)
    data_clean.(response_vars{i}) = log(data_clean.(response_vars{i}) + 1);
end

%% Step 2: Lasso Feature Selection and GAM Model Fitting
selected_predictors = struct();  % Store selected predictors
models = struct();  % Store fitted models
feature_importance_table = [];  % Initialize table for storing feature importance

for i = 1:length(response_vars)
    response_var = response_vars{i};
    
    % Convert predictor and response data to numeric arrays
    predictor_data = table2array(data_clean(:, predictor_vars));
    response_data = data_clean.(response_var);
    
    % Lasso feature selection (with 10-fold cross-validation)
    [B, FitInfo] = lasso(predictor_data, response_data, 'CV', 10);
    selected_idx = find(B(:, FitInfo.IndexMinMSE) ~= 0);  % Non-zero coefficients
    selected_predictors.(response_var) = predictor_vars(selected_idx);
    
    % Display selected predictors
    fprintf('Selected Predictors for %s:\n', response_var);
    disp(selected_predictors.(response_var));
    
    % Fit the GAM model with only the selected predictors
    if ~isempty(selected_idx)
        reduced_data = predictor_data(:, selected_idx);  % Use only selected predictors
        gam_model = fitrgam(reduced_data, response_data, 'CrossVal', 'on', 'KFold', 10, 'Verbose', 1, 'PredictorNames', predictor_vars(selected_idx));
        models.(response_var) = gam_model.Trained{1};  % Store the trained model

        % Step 5: Extract Lasso Coefficients, Signs, and Partial Dependence Importance
        lasso_coeffs = B(:, FitInfo.IndexMinMSE);  % Lasso coefficients
        selected_idx = find(lasso_coeffs ~= 0);  % Indices of selected features

        % Create a temporary table to hold the feature importance for this taxon
        temp_table = table();
        temp_table.Taxon = repmat({response_var}, length(selected_idx), 1);  % Create rows for each selected predictor
        temp_table.Predictor = predictor_vars(selected_idx)';  % Selected predictors

        % Check the sign of each coefficient and assign "Positive" or "Negative"
        temp_table.Effect = cell(length(selected_idx), 1);
        for j = 1:length(selected_idx)
            if lasso_coeffs(selected_idx(j)) > 0
                temp_table.Effect{j} = 'Positive';
            else
                temp_table.Effect{j} = 'Negative';
            end
        end

        % Compute partial dependence standard deviation as importance
        temp_table.Importance = zeros(length(selected_idx), 1);  % Initialize importance column
        for j = 1:length(selected_idx)
            % Compute the partial dependence for the j-th predictor
            [pdX, pdY] = partialDependence(models.(response_var), j, table2array(data_clean(:, selected_predictors.(response_var))));
            
            % Use standard deviation of the partial dependence output as an importance measure
            temp_table.Importance(j) = std(pdY);
        end

        % Append to the main feature importance table
        feature_importance_table = [feature_importance_table; temp_table];
    end
end

%% Save feature importance as CSV
output_csv_filename = 'GAM_feature_importance_table.csv';
writetable(feature_importance_table, output_csv_filename);
fprintf('Feature importance table saved to %s\n', output_csv_filename);


%% Step 3: Partial Dependence Plots for GAMs
% Define the publication settings
predictor_titles = {'Temperature', 'BEUTI', 'Alongshore Wind Speed', ...
                    'Surface PAR', 'PM2.5 SLV', 'Nitrate', 'Phosphate', 'Silicate', ...
                    'Ammonium', 'Outflow', 'Year-Month'};

% Set color scheme for each taxa
colorz = [
    0.2, 0.6, 0.8;  % Light blue
    0.9, 0.5, 0.1;  % Orange
    0.4, 0.7, 0.2;  % Light green
    0.8, 0.2, 0.4;  % Pink
    0.2, 0.4, 0.8;  % Dark blue
    0.7, 0.2, 0.9;  % Purple
    0.3, 0.7, 0.6   % Soft blue-green 
];

% Font and line width settings for publication
ftsz = 10;  % Font size
ftname = 'Helvetica';  % Font name

% Loop through response variables to create models and plot partial dependence plots
for i = 1:length(response_vars)
    response_var = response_vars{i};

    % Check if the model for this taxa exists
    if isfield(models, response_var)
        % Get the list of retained predictors for this taxa (from Lasso or best model)
        retained_predictors = selected_predictors.(response_var);

        % Subset the data to include only the retained predictors for this model
        retained_predictor_data = data_clean(:, retained_predictors);

        % Create figure for each taxa with publication-quality formatting
        figure;
        plot_idx = 1;  % Initialize subplot counter

        % Loop through each retained predictor to plot partial dependence
        for j = 1:length(retained_predictors)
            predictor_name = retained_predictors{j};

            subplot(ceil(length(retained_predictors)/2), 2, plot_idx);  % Create subplot

            % Select color for the current taxa
            color_idx = mod(i - 1, size(colorz, 1)) + 1;
            line_color = colorz(color_idx, :);

            % Find the index of the current predictor
            predictor_idx = find(strcmp(predictor_name, retained_predictors));

            % Use plotPartialDependence to generate the plot
            plotPartialDependence(models.(response_var), predictor_idx, retained_predictor_data);

            % Manually apply the color by extracting the line from the plot
            line_handle = findobj(gca, 'Type', 'Line');
            set(line_handle, 'Color', line_color, 'LineWidth', 1.5);  % Manually set the color and line width

            % Add titles and labels for the subplot
            title(predictor_titles{strcmp(predictor_name, predictor_vars)}, 'FontSize', ftsz, 'FontName', ftname);
            ylabel('Partial Dependence', 'FontSize', ftsz, 'FontName', ftname);

            % Increment subplot counter
            plot_idx = plot_idx + 1;
        end

        % Set the super title (sgtitle) to the response variable (taxa)
        sgtitle(response_var, 'FontSize', ftsz + 2, 'FontName', ftname, 'FontWeight', 'bold');

        % Set overall figure formatting for L&O
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 7.16, 10]);  % Double-column width
        set(gcf, 'PaperPositionMode', 'auto');  % Ensure correct figure size for paper

        % Save figure if required (uncomment if saving is needed)
        print(gcf, '-dpng', '-r300', ['GAM_partial_dependence_', response_var, '.png']);
        print(gcf, '-dpdf', '-r300', ['GAM_partial_dependence_', response_var, '.pdf']);
    end
end

%% Step 4: Quantitative Feature Importance
feature_importance_table = [];

for i = 1:length(response_vars)
    response_var = response_vars{i};
    
    if isfield(models, response_var)
        % Get the list of selected predictors for this response variable
        selected_predictors_for_var = selected_predictors.(response_var);
        
        % Subset the data to include only the selected predictors
        final_gam = models.(response_var);
        retained_predictor_data = data_clean(:, selected_predictors_for_var);
        
        % Initialize a temporary table to store importance values
        temp_table = table();
        temp_table.Taxon = repmat({response_var}, length(selected_predictors_for_var), 1);  % Rows for each predictor
        temp_table.Predictor = selected_predictors_for_var';  % Predictor names

        % Compute partial dependence standard deviation as importance
        temp_table.Importance = zeros(length(selected_predictors_for_var), 1);  % Initialize importance column
        for j = 1:length(selected_predictors_for_var)
            % Compute the partial dependence for the j-th predictor
            [pdX, pdY] = partialDependence(final_gam, j, retained_predictor_data);
            
            % Use standard deviation of the partial dependence output as an importance measure
            temp_table.Importance(j) = std(pdY);
        end

        % Extract the coefficients from the GAM model to determine the effect sign
        gam_coeffs = final_gam.Beta;  % Coefficients of the final GAM
        
        % Add the sign of the coefficient to indicate Positive or Negative effect
        temp_table.Effect = cell(length(selected_predictors_for_var), 1);
        for j = 1:length(selected_predictors_for_var)
            if gam_coeffs(j) > 0
                temp_table.Effect{j} = 'Positive';
            else
                temp_table.Effect{j} = 'Negative';
            end
        end

        % Append to the main feature importance table
        feature_importance_table = [feature_importance_table; temp_table];
    end
end

% Save feature importance as CSV
output_csv_filename = 'GAM_feature_importance_table.csv';
writetable(feature_importance_table, output_csv_filename);
fprintf('Feature importance table saved to %s\n', output_csv_filename);


%% Step 5: Residual Diagnostics for Each GAM Model
for i = 1:length(response_vars)
    response_var = response_vars{i};
    
    if isfield(models, response_var)
        selected_idx = find(ismember(predictor_vars, selected_predictors.(response_var)));
        reduced_data = table2array(data_clean(:, selected_predictors.(response_var)));
        predicted_values = predict(models.(response_var), reduced_data);
        residuals = data_clean.(response_var) - predicted_values;
        
        % Residual Diagnostics Plots
        figure;
        subplot(2, 2, 1);
        scatter(predicted_values, residuals);
        xlabel('Predicted Values');
        ylabel('Residuals');
        title(['Residual Plot for ', response_var]);
        grid on;

        subplot(2, 2, 2);
        histogram(residuals, 20);
        xlabel('Residuals');
        ylabel('Frequency');
        title('Histogram of Residuals');
        grid on;

        subplot(2, 2, 3);
        qqplot(residuals);
        title('Q-Q Plot of Residuals');

        subplot(2, 2, 4);
        plot(residuals);
        xlabel('Observation Index');
        ylabel('Residuals');
        title('Residuals over Time');
        grid on;
    end
end























%% GAM v4
%% GAM Analysis with Lasso Feature Selection, Partial Dependence, and Residual Diagnostics

% Step 1: Data Preparation
data = physical_forcings_biology_table;  % Load dataset
data.yearmonth = datenum(datetime(data.datetime.Year, data.datetime.Month, 1, 'TimeZone', 'UTC'));  % Add year-month variable

% Define predictor and response variables
predictor_vars = {'temperature', 'BEUTI', 'wspd_along', ...
                  'surface_downwelling_photosynthetic_photon_flux_in_air',  'pm2_5_sc_pm2_5_daily_slv_tt', ...
                  'Nitrate', 'Phosphate', 'Silicate', 'Ammonium', 'outflow', 'yearmonth'};
response_vars = {'Asterionellopsis', 'Centric', 'Skeletonema', 'Hemiaulus', 'Leptocylindrus', 'Thalassionema','Thalassiosira'};


% Clean data (remove rows with missing values)
data_clean = timetable2table(rmmissing(data(:, [predictor_vars, response_vars])));  
data_clean.datetime = [];  % Remove datetime column

% Log-transform the response variables (avoid log(0) by adding small constant)
for i = 1:length(response_vars)
    data_clean.(response_vars{i}) = log(data_clean.(response_vars{i}) + 1);
end

%% Step 2: Lasso Feature Selection and GAM Model Fitting
selected_predictors = struct();  % Store selected predictors
models = struct();  % Store fitted models

for i = 1:length(response_vars)
    response_var = response_vars{i};
    
    % Convert predictor and response data to numeric arrays
    predictor_data = table2array(data_clean(:, predictor_vars));
    response_data = data_clean.(response_var);
    
    % Lasso feature selection (with 10-fold cross-validation)
    [B, FitInfo] = lasso(predictor_data, response_data, 'CV', 10);
    selected_idx = find(B(:, FitInfo.IndexMinMSE) ~= 0);  % Non-zero coefficients
    selected_predictors.(response_var) = predictor_vars(selected_idx);
    
    % Display selected predictors
    fprintf('Selected Predictors for %s:\n', response_var);
    disp(selected_predictors.(response_var));
    
    % Fit the GAM model with only the selected predictors
    if ~isempty(selected_idx)
        reduced_data = predictor_data(:, selected_idx);  % Use only selected predictors
        gam_model = fitrgam(reduced_data, response_data, 'CrossVal', 'on', 'KFold', 10, 'Verbose', 1, 'PredictorNames', predictor_vars(selected_idx));
        models.(response_var) = gam_model.Trained{1};  % Store the trained model
    end
end

%% Step 3: Partial Dependence Plots for GAMs
% Define the publication settings
predictor_titles = {'Temperature', 'BEUTI', 'Alongshore Wind Speed', ...
                    'Surface PAR', 'PM2.5 SLV', 'Nitrate', 'Phosphate', 'Silicate', ...
                    'Ammonium', 'Outflow', 'Year-Month'};

% X-axis labels (with units)
x_labels = {'Temperature [°C]', 'BEUTI', 'Alongshore Wind Speed [m/s]', ...
            'Surface PAR [μmol photons m⁻² s⁻¹]', 'PM2.5 SLV [μg/m³]', 'Nitrate [μmol/L]', ...
            'Phosphate [μmol/L]', 'Silicate [μmol/L]', 'Ammonium [μmol/L]', ...
            'Outflow [m³/s]', 'Year-Month'};

% Set color scheme for each taxa
colorz = [
    0.2, 0.6, 0.8;  % Light blue
    0.9, 0.5, 0.1;  % Orange
    0.4, 0.7, 0.2;  % Light green
    0.8, 0.2, 0.4;  % Pink
    0.2, 0.4, 0.8;  % Dark blue
    0.7, 0.2, 0.9;  % Purple
    0.3, 0.7, 0.6   % Soft blue-green 
];

% Font and line width settings for publication
ftsz = 10;  % Font size
ftname = 'Helvetica';  % Font name

% Response variables (taxa)
response_vars = {'Asterionellopsis', 'Centric', 'Skeletonema', 'Hemiaulus', ...
                 'Leptocylindrus', 'Thalassionema', 'Thalassiosira'};

% Loop through response variables to create models and plot partial dependence plots
for i = 1:length(response_vars)
    response_var = response_vars{i};

    % Check if the model for this taxa exists (in case some taxa are not modeled)
    if isfield(models, response_var)
        % Get the list of retained predictors for this taxa (from Lasso or best model)
        retained_predictors = selected_predictors.(response_var);

        % Subset the data to include only the retained predictors for this model
        retained_predictor_data = data_clean(:, retained_predictors);

        % Create figure for each taxa with publication-quality formatting
        figure;
        plot_idx = 1;  % Initialize subplot counter

        % Loop through each retained predictor to plot partial dependence
        for j = 1:length(retained_predictors)
            predictor_name = retained_predictors{j};

            subplot(ceil(length(retained_predictors)/2), 2, plot_idx);  % Create subplot

            % Select color for the current taxa
            color_idx = mod(i - 1, size(colorz, 1)) + 1;
            line_color = colorz(color_idx, :);

            % Find the index of the current predictor
            predictor_idx = find(strcmp(predictor_name, retained_predictors));

            % Use plotPartialDependence to generate the plot (with the feature and data)
            % Now we pass the predictor data alongside the predictor index
            plotPartialDependence(models.(response_var), predictor_idx, retained_predictor_data);

            % Manually apply the color by extracting the line from the plot
            line_handle = findobj(gca, 'Type', 'Line');
            set(line_handle, 'Color', line_color, 'LineWidth', 1.5);  % Manually set the color and line width

            % Add titles and labels for the subplot
            title(predictor_titles{strcmp(predictor_name, predictor_vars)}, 'FontSize', ftsz, 'FontName', ftname);
            % xlabel(x_labels{strcmp(predictor_name, predictor_vars)}, 'FontSize', ftsz, 'FontName', ftname);

            % Set y-label as Partial Dependence (no specific unit)
            ylabel(["Partial"+newline+"Dependence"], 'FontSize', ftsz, 'FontName', ftname);


            % Format x-axis for 'yearmonth' if applicable
            if strcmp(predictor_name, 'yearmonth')
                datetick('x', 'yyyy');
                xlabel('Year-Month', 'FontSize', ftsz, 'FontName', ftname);
            end

            % Set axis properties for publication
            set(gca, 'FontSize', ftsz, 'FontName', ftname);

            plot_idx = plot_idx + 1;  % Increment subplot counter
        end

        
         % Set the super title (sgtitle) to the response variable (taxa)
        sgtitle(response_var, 'FontSize', ftsz + 2, 'FontName', ftname, 'FontWeight', 'bold');


        % Set overall figure formatting for L&O
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 7.16, 10]);  % Double-column width
        set(gcf, 'PaperPositionMode', 'auto');  % Ensure correct figure size for paper

        % Save figure if required (uncomment if saving is needed)
        print(gcf, '-dpng', '-r300', ['C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\GAM_partial_dependence\GAM_partial_dependence_', response_var, '.png']);
        print(gcf, '-dpdf', '-r300', ['C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\GAM_partial_dependence\GAM_partial_dependence_', response_var, '.pdf']);
    end
end
%% Step 4: Quantitative Feature Importance
feature_importance_table = [];

for i = 1:length(response_vars)
    response_var = response_vars{i};
    
    if isfield(models, response_var)
        % Get the list of selected predictors for this response variable
        selected_predictors_for_var = selected_predictors.(response_var);
        
        % Find the indices of the selected predictors in the full predictor list
        selected_idx = find(ismember(predictor_vars, selected_predictors_for_var));
        
        % Subset the GAM model to use only the selected predictors
        final_gam = models.(response_var);
        
        % Subset the data to include only the selected predictors (ensure it's a table)
        retained_predictor_data = data_clean(:, selected_predictors_for_var);
        
        % Ensure 'Data' is a table and matches the predictor names
        retained_predictor_data.Properties.VariableNames = final_gam.PredictorNames;
        
        % Initialize a temporary table to store importance values
        temp_table = table();
        temp_table.Taxon = repmat({response_var}, length(selected_idx), 1);  % Rows for each predictor
        temp_table.Predictor = selected_predictors_for_var';  % Predictor names

        % Compute partial dependence standard deviation as importance
        temp_table.Importance = zeros(length(selected_idx), 1);  % Initialize importance column
        for j = 1:length(selected_idx)
            % Compute the partial dependence for the j-th predictor
            [pdX, pdY] = partialDependence(final_gam, j, retained_predictor_data);
            
            % Use standard deviation of the partial dependence output as an importance measure
            temp_table.Importance(j) = std(pdY);
        end

        % Extract the coefficients from the GAM model to determine the effect sign
        gam_coeffs = final_gam.Beta;  % Coefficients of the final GAM
        
        % Add the sign of the coefficient to indicate Positive or Negative effect
        temp_table.Effect = cell(length(selected_idx), 1);
        for j = 1:length(selected_idx)
            if gam_coeffs(selected_idx(j)) > 0
                temp_table.Effect{j} = 'Positive';
            else
                temp_table.Effect{j} = 'Negative';
            end
        end

        % Append to the main feature importance table
        feature_importance_table = [feature_importance_table; temp_table];
    end
end

% Save feature importance as CSV
output_csv_filename = 'GAM_feature_importance_table.csv';
writetable(feature_importance_table, output_csv_filename);
fprintf('Feature importance table saved to %s\n', output_csv_filename);

%% %% Initialize tables to store results
parametric_coeffs_table = table();  % Table for parametric coefficients
model_fitness_table = table();      % Table for model fitness
predictor_significance_table = table();  % Table for predictor significance

% Loop through each response variable
for i = 1:length(response_vars)
    response_var = response_vars{i};

    % Check if the model for this response variable exists
    if isfield(models, response_var)
        final_gam = models.(response_var);

        % Table 1: Parametric Coefficients for Intercept
        % Extract the constant term (intercept)
        intercept = final_gam.Intercept;  % Use 'Intercept' instead of 'Bias'
        
        % Append to parametric coefficients table
        parametric_coeffs_table = [parametric_coeffs_table; ...
            table({response_var}, intercept, ...
            'VariableNames', {'ResponseTaxa', 'Intercept'})];

        % Table 2: Model Fitness
        adj_rsq = final_gam.Rsquared.Adjusted;  % Adjusted R-squared
        deviance_explained = 100 * final_gam.DevianceExplained;  % Percent deviance explained
        gcv = final_gam.ModelCriterion.GCV;  % Generalized Cross-Validation score

        % Append to model fitness table
        model_fitness_table = [model_fitness_table; ...
            table({response_var}, adj_rsq, deviance_explained, gcv, ...
            'VariableNames', {'ResponseTaxa', 'AdjustedRsq', 'DevianceExplained', 'GCV'})];

        % Table 3: Statistical Significance for Predictors
        % Get effective degrees of freedom (EDF), F-values, and p-values for each predictor
        edf = final_gam.EffectiveDF;  % Effective Degrees of Freedom
        f_values = final_gam.FTest.FStatistic;  % F-values for each predictor
        p_values = final_gam.FTest.pValue;  % p-values for each predictor
        
        % Loop over each predictor
        for j = 1:length(final_gam.PredictorNames)
            predictor_name = final_gam.PredictorNames{j};
            predictor_edf = edf(j);
            predictor_f_value = f_values(j);
            predictor_p_value = p_values(j);
            
            % Append to predictor significance table
            predictor_significance_table = [predictor_significance_table; ...
                table({response_var}, {predictor_name}, predictor_edf, predictor_f_value, predictor_p_value, ...
                'VariableNames', {'ResponseTaxa', 'Predictor', 'EDF', 'FValue', 'PValue'})];
        end
    end
end

%% Save the tables as CSV files (optional)
writetable(parametric_coeffs_table, 'GAM_parametric_coeffs_table.csv');
writetable(model_fitness_table, 'GAM_model_fitness_table.csv');
writetable(predictor_significance_table, 'GAM_predictor_significance_table.csv');

fprintf('Tables saved successfully.\n');


%% Step 5: Residual Diagnostics for Each GAM Model
for i = 1:length(response_vars)
    response_var = response_vars{i};
    
    if isfield(models, response_var)
        selected_idx = find(ismember(predictor_vars, selected_predictors.(response_var)));
        reduced_data = table2array(data_clean(:, selected_predictors.(response_var)));
        predicted_values = predict(models.(response_var), reduced_data);
        residuals = data_clean.(response_var) - predicted_values;
        
        % Residual Diagnostics Plots
        figure;
        subplot(2, 2, 1);
        scatter(predicted_values, residuals);
        xlabel('Predicted Values');
        ylabel('Residuals');
        title(['Residual Plot for ', response_var]);
        grid on;

        subplot(2, 2, 2);
        histogram(residuals, 20);
        xlabel('Residuals');
        ylabel('Frequency');
        title('Histogram of Residuals');
        grid on;

        subplot(2, 2, 3);
        qqplot(residuals);
        title('Q-Q Plot of Residuals');

        subplot(2, 2, 4);
        plot(residuals);
        xlabel('Observation Index');
        ylabel('Residuals');
        title('Residuals over Time');
        grid on;
    end
end





















%% GAM v3
data = physical_forcings_biology_table;
data.yearmonth = datenum(datetime(data.datetime.Year, data.datetime.Month, 1, 'TimeZone', 'UTC'));

%% Step 1: Define predictor and response variables
predictor_vars = {'sea_water_temperature', 'BEUTI', 'wspd_along', 'AOD_500nm', ...
                  'surface_downwelling_photosynthetic_photon_flux_in_air',  'pm2_5_sc_pm2_5_daily_slv_tt', ...
                  'Nitrate', 'Phosphate', 'Silicate', 'Ammonium', 'outflow', 'yearmonth'};
response_vars = {'Asterionellopsis', 'Centric', 'Skeletonema', 'Hemiaulus', 'Leptocylindrus', 'Thalassionema'};

%% Step 2: Clean data
data_clean = timetable2table(rmmissing(data(:, [predictor_vars, response_vars])));  % Remove rows with missing data
data_clean.datetime = [];  % Remove datetime column if not needed

%% Step 3: Log-transform the response variables
for i = 1:length(response_vars)
    % Add a small constant to avoid log(0)
    data_clean.(response_vars{i}) = log(data_clean.(response_vars{i}) + 1);
end

%% Step 4: Use Lasso to Perform Feature Selection for Each Response Variable
selected_predictors = struct();  % Store selected predictors for each response
models = struct();               % Store models for each response
for i = 1:length(response_vars)
    response_var = response_vars{i};

    % Convert data to numeric format
    predictor_data = table2array(data_clean(:, predictor_vars));
    response_data = data_clean.(response_var);
    
    % Lasso feature selection (with 10-fold cross-validation to find best lambda)
    [B, FitInfo] = lasso(predictor_data, response_data, 'CV', 10);
    
    % Extract non-zero coefficients (selected predictors)
    selected_idx = find(B(:, FitInfo.IndexMinMSE) ~= 0);
    selected_predictors.(response_var) = predictor_vars(selected_idx);
    
    % Display selected predictors
    fprintf('Selected Predictors for %s:\n', response_var);
    disp(selected_predictors.(response_var));
    
    % Refit the GAM with only the selected predictors
    if ~isempty(selected_idx)
        reduced_data = predictor_data(:, selected_idx);  % Use only selected predictors
        gam_model = fitrgam(reduced_data, response_data, 'CrossVal', 'on', 'KFold', 10, 'Verbose', 1, 'PredictorNames', predictor_vars(selected_idx));
        models.(response_var) = gam_model.Trained{1};  % Store the trained model
    end
end

%%
% Define the publication settings
predictor_titles = {'Temperature', 'BEUTI', 'Alongshore Wind Speed', ...
                    'Surface PAR', 'PM2.5 SLV', 'Nitrate', 'Phosphate', 'Silicate', ...
                    'Ammonium', 'Outflow', 'Year-Month'};

% X-axis labels (with units)
x_labels = {'Temperature [°C]', 'BEUTI', 'Alongshore Wind Speed [m/s]', ...
            'Surface PAR [μmol photons m⁻² s⁻¹]', 'PM2.5 SLV [μg/m³]', 'Nitrate [μmol/L]', ...
            'Phosphate [μmol/L]', 'Silicate [μmol/L]', 'Ammonium [μmol/L]', ...
            'Outflow [m³/s]', 'Year-Month'};

% Set color scheme for each taxa
colorz = [
    0.2, 0.6, 0.8;  % Light blue
    0.9, 0.5, 0.1;  % Orange
    0.4, 0.7, 0.2;  % Light green
    0.8, 0.2, 0.4;  % Pink
    0.2, 0.4, 0.8;  % Dark blue
    0.7, 0.2, 0.9;  % Purple
    0.3, 0.7, 0.6   % Soft blue-green 
];

% Font and line width settings for publication
ftsz = 10;  % Font size
ftname = 'Helvetica';  % Font name

% Response variables (taxa)
response_vars = {'Asterionellopsis', 'Centric', 'Skeletonema', 'Hemiaulus', ...
                 'Leptocylindrus', 'Thalassionema', 'Thalassiosira'};

% Loop through response variables to create models and plot partial dependence plots
for i = 1:length(response_vars)
    response_var = response_vars{i};

    % Check if the model for this taxa exists (in case some taxa are not modeled)
    if isfield(models, response_var)
        % Get the list of retained predictors for this taxa (from Lasso or best model)
        retained_predictors = selected_predictors.(response_var);

        % Subset the data to include only the retained predictors for this model
        retained_predictor_data = data_clean(:, retained_predictors);

        % Create figure for each taxa with publication-quality formatting
        figure;
        plot_idx = 1;  % Initialize subplot counter

        % Loop through each retained predictor to plot partial dependence
        for j = 1:length(retained_predictors)
            predictor_name = retained_predictors{j};

            subplot(ceil(length(retained_predictors)/2), 2, plot_idx);  % Create subplot

            % Select color for the current taxa
            color_idx = mod(i - 1, size(colorz, 1)) + 1;
            line_color = colorz(color_idx, :);

            % Find the index of the current predictor
            predictor_idx = find(strcmp(predictor_name, retained_predictors));

            % Use plotPartialDependence to generate the plot (with the feature and data)
            % Now we pass the predictor data alongside the predictor index
            plotPartialDependence(models.(response_var), predictor_idx, retained_predictor_data);

            % Manually apply the color by extracting the line from the plot
            line_handle = findobj(gca, 'Type', 'Line');
            set(line_handle, 'Color', line_color, 'LineWidth', 1.5);  % Manually set the color and line width

            % Add titles and labels for the subplot
            title(predictor_titles{strcmp(predictor_name, predictor_vars)}, 'FontSize', ftsz, 'FontName', ftname);
            xlabel(x_labels{strcmp(predictor_name, predictor_vars)}, 'FontSize', ftsz, 'FontName', ftname);

            % Set y-label as Partial Dependence (no specific unit)
            ylabel('Partial Dependence', 'FontSize', ftsz, 'FontName', ftname);


            % Format x-axis for 'yearmonth' if applicable
            if strcmp(predictor_name, 'yearmonth')
                datetick('x', 'yyyy');
                xlabel('Year-Month', 'FontSize', ftsz, 'FontName', ftname);
            end

            % Set axis properties for publication
            set(gca, 'FontSize', ftsz, 'FontName', ftname);

            plot_idx = plot_idx + 1;  % Increment subplot counter
        end

        
         % Set the super title (sgtitle) to the response variable (taxa)
        sgtitle(response_var, 'FontSize', ftsz + 2, 'FontName', ftname, 'FontWeight', 'bold');


        % Set overall figure formatting for L&O
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 7.16, 10]);  % Double-column width
        set(gcf, 'PaperPositionMode', 'auto');  % Ensure correct figure size for paper

        % Save figure if required (uncomment if saving is needed)
        print(gcf, '-dpng', '-r300', ['C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\GAM_partial_dependence\GAM_partial_dependence_', response_var, '.png']);
        print(gcf, '-dpdf', '-r300', ['C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\GAM_partial_dependence\GAM_partial_dependence_', response_var, '.pdf']);
    end
end


%%
% Step 1: Initialize the table for storing feature importance results
feature_importance_table = [];

% Step 2: Define the phytoplankton taxa (response variables)
response_vars = {'Asterionellopsis', 'Centric', 'Skeletonema', 'Hemiaulus', ...
                 'Leptocylindrus', 'Thalassionema', 'Thalassiosira'};

% Step 3: Define the predictor variables
predictor_vars = {'temperature', 'BEUTI', 'wspd_along', ...
                  'surface_downwelling_photosynthetic_photon_flux_in_air', 'pm2_5_sc_pm2_5_daily_slv_tt', ...
                  'Nitrate', 'Phosphate', 'Silicate', 'Ammonium', 'outflow', 'yearmonth'};

% Loop over each taxon
for i = 1:length(response_vars)
    response_var = response_vars{i};

    % Step 4: Get the selected predictors for this taxon from the Lasso model
    % Rename to avoid overwriting the main 'selected_predictors' structure
    selected_predictors_for_taxon = selected_predictors.(response_var);  % Get selected predictors for this taxon
    
    % Step 5: Extract the Lasso coefficients and their signs
    lasso_coeffs = B(:, FitInfo.IndexMinMSE);  % Lasso coefficients
    selected_idx = find(lasso_coeffs ~= 0);  % Indices of selected features
    
    % Create a temporary table to hold the feature importance for this taxon
    temp_table = table();
    temp_table.Taxon = repmat({response_var}, length(selected_idx), 1);  % Create rows for each selected predictor
    temp_table.Predictor = predictor_vars(selected_idx)';  % Selected predictors
    
    % Check the sign of each coefficient and assign "Positive" or "Negative"
    temp_table.Effect = cell(length(selected_idx), 1);
    for j = 1:length(selected_idx)
        if lasso_coeffs(selected_idx(j)) > 0
            temp_table.Effect{j} = 'Positive';
        else
            temp_table.Effect{j} = 'Negative';
        end
    end
    
    % Step 6: Estimate Predictor Importance using partial dependence plots
    final_gam = fitrgam(data_clean(:, selected_idx), data_clean.(response_var), 'Verbose', 1);
    % Store the GAM model for this taxon
    models.(response_var) = final_gam;

    % Assign true to all rows of Importance column
    temp_table.Importance = repmat(true, height(temp_table), 1);  % All selected variables are important

    % Step 7: Append the results to the overall feature importance table
    feature_importance_table = [feature_importance_table; temp_table];
end

% Step 8: Save the Feature Importance Table as a CSV file
output_csv_filename = 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\tables\feature_importance_table.csv';
writetable(feature_importance_table, output_csv_filename);

% Display confirmation message
fprintf('Feature importance table saved to %s\n', output_csv_filename);



%% Step 6: Residual Diagnostics for Reduced Models
for i = 1:length(response_vars)
    response_var = response_vars{i};
    
    if isfield(models, response_var)
        % Predict and calculate residuals using the final trained model
        selected_idx = find(ismember(predictor_vars, selected_predictors.(response_var)));
        reduced_data = table2array(data_clean(:, selected_predictors.(response_var)));  % Use selected predictors data
        predicted_values = predict(models.(response_var), reduced_data);
        residuals = data_clean.(response_var) - predicted_values;
        
        % Residual Diagnostics Plots
        figure;
        subplot(2, 2, 1);
        scatter(predicted_values, residuals);
        xlabel('Predicted Values');
        ylabel('Residuals');
        title(['Residual Plot for ', response_var]);
        grid on;

        subplot(2, 2, 2);
        histogram(residuals, 20);
        xlabel('Residuals');
        ylabel('Frequency');
        title('Histogram of Residuals');
        grid on;

        subplot(2, 2, 3);
        qqplot(residuals);
        title('Q-Q Plot of Residuals');

        subplot(2, 2, 4);
        plot(residuals);
        xlabel('Observation Index');
        ylabel('Residuals');
        title('Residuals over Time');
        grid on;
    end
end



%%

% Define the publication settings
predictor_titles = {'Temperature', 'BEUTI', 'Alongshore Wind Speed', ...
                    'Surface PAR', 'PM2.5 SLV', 'Nitrate', 'Phosphate', 'Silicate', ...
                    'Ammonium', 'Outflow', 'Year-Month'};

% X-axis labels (with units)
x_labels = {'Temperature [°C]', 'BEUTI', 'Alongshore Wind Speed [m/s]', ...
            'Surface PAR [μmol photons m⁻² s⁻¹]', 'PM2.5 SLV [μg/m³]', 'Nitrate [μmol/L]', ...
            'Phosphate [μmol/L]', 'Silicate [μmol/L]', 'Ammonium [μmol/L]', ...
            'Outflow [m³/s]', 'Year-Month'};

% Set color scheme for each taxa
colorz = [
    0.2, 0.6, 0.8;  % Light blue
    0.9, 0.5, 0.1;  % Orange
    0.4, 0.7, 0.2;  % Light green
    0.8, 0.2, 0.4;  % Pink
    0.2, 0.4, 0.8;  % Dark blue
    0.7, 0.2, 0.9;  % Purple
    0.3, 0.7, 0.6   % Soft blue-green 
];

% Font and line width settings
ftsz = 10;  % Font size
ftname = 'Helvetica';  % Font name
linewdt = 1.5;  % Line width for plot elements

% Response variables (taxa)
response_vars = {'Asterionellopsis', 'Centric', 'Skeletonema', 'Hemiaulus', ...
                 'Leptocylindrus', 'Thalassionema', 'Thalassiosira'};

% Loop through response variables to create models and plot PDPs
for i = 1:length(response_vars)
    response_var = response_vars{i};

    % Check if the model for this taxa exists (in case some taxa are not modeled)
    if isfield(models, response_var)
        % Get the list of retained predictors for this taxa (from Lasso or best model)
        retained_predictors = selected_predictors.(response_var);

        % Subset the data to include only the retained predictors for this model
        retained_predictor_data = data_clean(:, retained_predictors);

        % Convert retained_predictor_data to a numeric array
        predictor_data_array = table2array(retained_predictor_data);

        % Create figure for each taxa with publication-quality formatting
        figure;
        plot_idx = 1;  % Initialize subplot counter

        % Loop through each retained predictor to plot partial dependence
        for j = 1:length(retained_predictors)
            predictor_name = retained_predictors{j};

            subplot(ceil(length(retained_predictors)/2), 2, plot_idx);  % Create subplot

            % Select color for the current taxa
            color_idx = mod(i - 1, size(colorz, 1)) + 1;
            line_color = colorz(color_idx, :);

            % Find the index of the current predictor
            predictor_idx = find(strcmp(predictor_name, retained_predictors));

            % Compute the partial dependence for the retained predictor
            [pdX, pdY] = partialDependence(models.(response_var), predictor_idx, 'Data', predictor_data_array);

            % Plot partial dependence with specified color and line width
            plot(pdX, pdY, 'Color', line_color, 'LineWidth', linewdt);
            hold on;

            % Add titles and labels for the subplot
            title(predictor_titles{strcmp(predictor_name, predictor_vars)}, 'FontSize', ftsz, 'FontName', ftname);
            xlabel(x_labels{strcmp(predictor_name, predictor_vars)}, 'FontSize', ftsz, 'FontName', ftname);

            % Format x-axis for 'yearmonth' if applicable
            if strcmp(predictor_name, 'yearmonth')
                datetick('x', 'yyyy');
                xlabel('Year-Month', 'FontSize', ftsz, 'FontName', ftname);
            end

            % Set axis properties for publication
            set(gca, 'FontSize', ftsz, 'FontName', ftname);

            plot_idx = plot_idx + 1;  % Increment subplot counter
        end

        % Set overall figure formatting for L&O
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 7.16, 10]);  % Double-column width
        set(gcf, 'PaperPositionMode', 'auto');  % Ensure correct figure size for paper

        % Save figure if required (uncomment if saving is needed)
        % print(gcf, '-dpng', '-r300', ['output_path_for_', response_var, '.png']);
        % print(gcf, '-dpdf', '-r300', ['output_path_for_', response_var, '.pdf']);
    end
end




%% Step 6: Residual Diagnostics for Reduced Models
for i = 1:length(response_vars)
    response_var = response_vars{i};
    
    if isfield(models, response_var)
        % Predict and calculate residuals using the final trained model
        selected_idx = find(ismember(predictor_vars, selected_predictors.(response_var)));
        reduced_data = table2array(data_clean(:, selected_predictors.(response_var)));  % Use selected predictors data
        predicted_values = predict(models.(response_var), reduced_data);
        residuals = data_clean.(response_var) - predicted_values;
        
        % Residual Diagnostics Plots
        figure;
        subplot(2, 2, 1);
        scatter(predicted_values, residuals);
        xlabel('Predicted Values');
        ylabel('Residuals');
        title(['Residual Plot for ', response_var]);
        grid on;

        subplot(2, 2, 2);
        histogram(residuals, 20);
        xlabel('Residuals');
        ylabel('Frequency');
        title('Histogram of Residuals');
        grid on;

        subplot(2, 2, 3);
        qqplot(residuals);
        title('Q-Q Plot of Residuals');

        subplot(2, 2, 4);
        plot(residuals);
        xlabel('Observation Index');
        ylabel('Residuals');
        title('Residuals over Time');
        grid on;
    end
end

%% Step 5: Partial Dependence Plots for the Reduced Models
for i = 1:length(response_vars)
    response_var = response_vars{i};
    
    if isfield(models, response_var)
        figure;
        selected_vars = selected_predictors.(response_var);
        
        % Plot Partial Dependence for the selected predictors
        for j = 1:length(selected_vars)
            predictor_idx = find(strcmp(predictor_vars, selected_vars{j}));  % Get the index of the selected predictor
            
            subplot(ceil(length(selected_vars)/2), 2, j);  % Create subplots for partial dependence
            plotPartialDependence(models.(response_var), predictor_idx, predictor_data);
            title(['Partial Dependence: ', selected_vars{j}, ' on ', response_var]);
        end
    end
end

%% Step 6: Residual Diagnostics for Reduced Models
for i = 1:length(response_vars)
    response_var = response_vars{i};
    
    if isfield(models, response_var)
        % Predict and calculate residuals using the final trained model
        selected_idx = find(ismember(predictor_vars, selected_predictors.(response_var)));
        predicted_values = predict(models.(response_var), predictor_data(:, selected_idx));
        residuals = data_clean.(response_var) - predicted_values;
        
        % Residual Diagnostics Plots
        figure;
        subplot(2, 2, 1);
        scatter(predicted_values, residuals);
        xlabel('Predicted Values');
        ylabel('Residuals');
        title(['Residual Plot for ', response_var]);
        grid on;

        subplot(2, 2, 2);
        histogram(residuals, 20);
        xlabel('Residuals');
        ylabel('Frequency');
        title('Histogram of Residuals');
        grid on;

        subplot(2, 2, 3);
        qqplot(residuals);
        title('Q-Q Plot of Residuals');

        subplot(2, 2, 4);
        plot(residuals);
        xlabel('Observation Index');
        ylabel('Residuals');
        title('Residuals over Time');
        grid on;
    end
end






















%% Step 4: Fit GAM Models for Each Response Variable
% Store the models and perform cross-validation
models = struct();
cv_results = struct();  % Store cross-validation results
for i = 1:length(response_vars)
    response_var = response_vars{i};
    
    % Convert the predictor data to numeric arrays
    predictor_data = table2array(data_clean(:, predictor_vars));
    response_data = data_clean.(response_var);
    
    % Fit the GAM model with cross-validation (10-fold)
    gam_model = fitrgam(predictor_data, response_data, 'CrossVal', 'on', 'KFold', 10, 'Verbose', 1, 'PredictorNames', predictor_vars);
    
    % Store the model and cross-validation results
    models.(response_var) = gam_model.Trained{1};  % Store the trained model
    cv_results.(response_var) = kfoldLoss(gam_model);  % Store the cross-validated loss (MSE)
    
    % Display model summary
    fprintf('GAM model for %s (Cross-validated MSE: %f):\n', response_var, cv_results.(response_var));
end

%% Step 5: Partial Dependence Plots Using GAMs
for i = 1:length(response_vars)
    response_var = response_vars{i};
    
    % Plot Partial Dependence for the predictors
    figure;
    for j = 1:length(predictor_vars)
        subplot(ceil(length(predictor_vars)/2), 2, j);  % Create subplots for partial dependence
        
        % Partial dependence plot requires the model, the feature index, and the input data
        plotPartialDependence(models.(response_var), j, predictor_data);
        title(['Partial Dependence: ', predictor_vars{j}, ' on ', response_var]);
    end
end


%% Step 6: Residual Diagnostics
for i = 1:length(response_vars)
    response_var = response_vars{i};
    
    % Predict and calculate residuals using the final trained model
    predicted_values = predict(models.(response_var), predictor_data);
    residuals = data_clean.(response_var) - predicted_values;
    
    % Residual Diagnostics Plots
    figure;
    subplot(2, 2, 1);
    scatter(predicted_values, residuals);
    xlabel('Predicted Values');
    ylabel('Residuals');
    title(['Residual Plot for ', response_var]);
    grid on;

    subplot(2, 2, 2);
    histogram(residuals, 20);
    xlabel('Residuals');
    ylabel('Frequency');
    title('Histogram of Residuals');
    grid on;

    subplot(2, 2, 3);
    qqplot(residuals);
    title('Q-Q Plot of Residuals');

    subplot(2, 2, 4);
    plot(residuals);
    xlabel('Observation Index');
    ylabel('Residuals');
    title('Residuals over Time');
    grid on;
end







