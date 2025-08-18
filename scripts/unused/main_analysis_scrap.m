%% SCRAP 


%% Q4: GAM SCRAP


%% GAM v2
data = physical_forcings_biology_table;
data.yearmonth = datenum(datetime(data.datetime.Year, data.datetime.Month, 1, 'TimeZone', 'UTC'));

%% Step 1: Define predictor and response variables
predictor_vars = {'sea_water_temperature', 'BEUTI', 'wspd_along', 'AOD_500nm', ...
                  'surface_downwelling_photosynthetic_photon_flux_in_air',  'pm2_5_sc_pm2_5_daily_slv_tt','Nitrate', ...
                  'Phosphate', 'Silicate', 'Ammonium','outflow','yearmonth'};
response_vars = {'Asterionellopsis', 'Centric', 'Skeletonema', 'Hemiaulus', 'Leptocylindrus', 'Thalassionema'};

%% Step 2: Clean data
data_clean = timetable2table(rmmissing(data(:, [predictor_vars, response_vars])));  % Remove rows with missing data
data_clean.datetime = [];  % Remove datetime column if not needed

%% Step 3: Log-transform the response variables
% Apply log transformation to avoid issues with negative or zero values
for i = 1:length(response_vars)
    % Add a small constant to avoid log(0)
    data_clean.(response_vars{i}) = log(data_clean.(response_vars{i}) + 1);
end

%% Step 4: Stepwise Feature Selection using AIC
% Fit an initial model with all predictors for each response variable
models = struct();
best_models = struct();  % Store the best models
for i = 1:length(response_vars)
    response_var = response_vars{i};

    % Stepwise regression to select the best predictors based on AIC
    best_models.(response_var) = stepwiselm(data_clean(:, predictor_vars), data_clean.(response_var), ...
        'Criterion', 'aic', 'Verbose', 1, 'PEnter', 0.05, 'PRemove', 0.10);

    % Store the best model
    models.(response_var) = best_models.(response_var);

    % Display model summary
    fprintf('Best stepwise model for %s:\n', response_var);
    disp(best_models.(response_var));
end

%% Manual K-Fold Cross-Validation for Stepwise Models

K = 10;  % Number of folds for cross-validation
n = height(data_clean);  % Number of observations
indices = randperm(n);  % Randomly shuffle the data indices

% Divide the shuffled indices into K roughly equal parts
fold_size = floor(n / K);
fold_indices = cell(K, 1);  % Preallocate cell array to hold indices for each fold

for k = 1:K
    if k == K
        % Last fold takes the remaining observations
        fold_indices{k} = indices((k-1)*fold_size + 1:end);
    else
        fold_indices{k} = indices((k-1)*fold_size + 1:k*fold_size);
    end
end

mse_values = zeros(K, length(response_vars));  % Preallocate array to store MSE for each fold
predictor_importance = zeros(length(predictor_vars), length(response_vars));  % To store predictor importance

for i = 1:length(response_vars)
    response_var = response_vars{i};
    fprintf('Processing response variable: %s\n', response_var);

    for k = 1:K
        % Split data into training and test sets
        test_idx = fold_indices{k};  % Logical index for test data
        train_idx = setdiff(indices, test_idx);  % Indices for training data

        % Training data
        train_data = data_clean(train_idx, :);

        % Test data
        test_data = data_clean(test_idx, :);

        % Fit the stepwise linear model on training data
        cv_model = stepwiselm(train_data(:, predictor_vars), train_data.(response_var), ...
            'Criterion', 'aic', 'Verbose', 0, 'PEnter', 0.05, 'PRemove', 0.10);

        % Predict on test data
        y_pred = predict(cv_model, test_data(:, predictor_vars));

        % Calculate mean squared error (MSE) on test data
        mse_values(k, i) = mse(test_data.(response_var), y_pred);

        % Estimate predictor importance by calculating permutation importance
        for j = 1:length(predictor_vars)
            permuted_data = test_data;
            permuted_data.(predictor_vars{j}) = permuted_data.(predictor_vars{j})(randperm(height(permuted_data)));  % Permute the predictor
            y_perm_pred = predict(cv_model, permuted_data(:, predictor_vars));
            permuted_mse = mse(test_data.(response_var), y_perm_pred);
            predictor_importance(j, i) = predictor_importance(j, i) + (permuted_mse - mse_values(k, i));  % Difference in MSE indicates importance
        end
    end
end

% Average MSE across all folds for each response variable
mean_mse = mean(mse_values, 1);

% Display the cross-validated mean squared error for each response
for i = 1:length(response_vars)
    fprintf('Cross-validated MSE for %s: %f\n', response_vars{i}, mean_mse(i));
end

% Calculate average predictor importance across all folds
avg_predictor_importance = mean(predictor_importance, 2);

% Display predictor importance
importance_table = table(predictor_vars', avg_predictor_importance, 'VariableNames', {'Predictor', 'Importance'});
disp('Predictor Importance Table:');
disp(importance_table);

%% Model Diagnostics for the Best Model
for i = 1:length(response_vars)
    response_var = response_vars{i};

    % Fit the best model on the entire dataset
    final_model = stepwiselm(data_clean(:, predictor_vars), data_clean.(response_var), ...
        'Criterion', 'aic', 'Verbose', 0, 'PEnter', 0.05, 'PRemove', 0.10);

    % Predict and calculate residuals
    predicted_values = predict(final_model, data_clean(:, predictor_vars));
    residuals = data_clean.(response_var) - predicted_values;

    % Residual Diagnostics
    figure;
    subplot(2, 2, 1);
    scatter(predicted_values, residuals);
    xlabel('Predicted Values');
    ylabel('Residuals');
    title(['Residual Plot for ', response_var]);
    grid on;

    % Histogram of residuals
    subplot(2, 2, 2);
    histogram(residuals, 20);
    xlabel('Residuals');
    ylabel('Frequency');
    title('Histogram of Residuals');
    grid on;

    % Q-Q plot for normality
    subplot(2, 2, 3);
    qqplot(residuals);
    title('Q-Q Plot of Residuals');

    % Residuals over time to check for independence
    subplot(2, 2, 4);
    plot(residuals);
    xlabel('Observation Index');
    ylabel('Residuals');
    title('Residuals over Time');
    grid on;

    % Check assumptions
    % 1. Homoscedasticity: Residual plot should show no patterns (random scatter).
    % 2. Normality: Residuals should be roughly normally distributed (histogram and Q-Q plot).
    % 3. Independence: Residuals over time should not show patterns (random scatter).
end



%% GAM
data = physical_forcings_biology_table;
data.yearmonth=datenum(datetime(data.datetime.Year,data.datetime.Month,1,'TimeZone','UTC'));

%% Define predictor and response variables
% Predictor columns
predictor_vars = {'sea_water_temperature', 'BEUTI', 'wspd_along', 'AOD_500nm', ...
                  'surface_downwelling_photosynthetic_photon_flux_in_air',  'pm2_5_sc_pm2_5_daily_slv_tt','Nitrate', ...
             'Phosphate', 'Silicate', 'Ammonium','outflow','yearmonth'};

% Response columns (biological variables)
response_vars = {'Asterionellopsis', 'Centric', 'Skeletonema', 'Hemiaulus', 'Leptocylindrus', 'Thalassionema'};

predictor_vars_labels={'Temperature [deg C]','BEUTI','Wind Speed [m/s]',"Aerosol Optical Dept"+newline+"at 500nm",'PM2.5 [μg/m³]',"PAR"+newline+"[μmol photons m⁻² s⁻¹]",...
    'PM2.5 [μg/m³]','Nitrate [μmol/L]','Phosphate [μmol/L]',...
    'Silicate [μmol/L]','Ammonium [μmol/L]', 'Outflow [m³/s]','Year-month'};


%% Keep only the relevant predictor and response columns
data_clean = rmmissing(data(:, response_vars));  % Remove rows with missing responses
data_clean = data_clean(~any(ismissing(data_clean), 2), :);  % Ensure no rows with missing predictor variables

% Keep only the relevant predictor and response columns
data_clean = data(:, [predictor_vars, response_vars]);





%% Step 1: Exploratory Data Analysis (EDA)
% Plot correlation matrix to check for multicollinearity
data_matrix = table2array(data_clean);
corr_matrix = corr(data_matrix, 'Rows', 'pairwise');
clf; figure(1)
imagesc(corr_matrix);
colorbar;
title('Correlation Matrix');
xticks(1:length([predictor_vars, response_vars]));
xticklabels([predictor_vars, response_vars]);
yticks(1:length([predictor_vars, response_vars]));
yticklabels([predictor_vars, response_vars]);
xtickangle(45);
set(gca, 'FontSize', 8);
colormap('default');
set(gcf,'Position',[0 0 1000 800])

%% Step 1: Remove rows with missing data in response variables
% Remove rows where any of the response variables (Asterionellopsis) are NaN
data_clean = data(~any(ismissing(data(:, response_vars)), 2), :);

% Ensure no rows have missing predictor variables as well
data_clean = data_clean(~any(ismissing(data_clean(:, predictor_vars)), 2), :);

%% Step 2: Keep only the relevant predictor and response columns
data_clean = data_clean(:, [predictor_vars, response_vars]);


%% Convert the table to numeric arrays
% Convert predictor variables to numeric matrix
predictor_data = table2array(data_clean(:, predictor_vars));

% Convert response variables to numeric vector
response_data = table2array(data_clean(:, response_vars{1}));

%% Step 3: Fit the GAM Model for each response variable
% Store models and diagnostics
models = struct();
for i = 1:length(response_vars)
    response_var = response_vars{i};
    
    % Fit the model using numeric data (X as predictors, Y as response)
    models.(response_var) = fitrgam(predictor_data, response_data);
    
    % Display model summary
    fprintf('GAM model for %s:\n', response_var);
    disp(models.(response_var));
end


%% Step 4: Model Diagnostics

% Residual plot for one of the response variables (Asterionellopsis)
response_var = 'Asterionellopsis';
model = models.(response_var);

% Predict values and calculate residuals
predicted_values = predict(model, data_clean{:, predictor_vars});
residuals = data_clean{:, response_var} - predicted_values;
% Plot residuals to check for patterns (should be randomly scattered)
figure(2);
subplot(2,2,1);
scatter(predicted_values, residuals);
xlabel('Predicted Values');
ylabel('Residuals');
title(['Residual Plot for ', response_var]);
grid on;

% Plot histogram of residuals to check for normality
subplot(2,2,2);
histogram(residuals, 20);
xlabel('Residuals');
ylabel('Frequency');
title('Histogram of Residuals');
grid on;

% Plot a Q-Q plot to check for normality of residuals
subplot(2,2,3);
qqplot(residuals);
title('Q-Q Plot of Residuals');

% Plot residuals over time (if applicable) to check for independence
subplot(2,2,4);
plot(residuals);
xlabel('Observation Index');
ylabel('Residuals');
title('Residuals over Time');
grid on;

%% Step 5: Check for multicollinearity (Variance Inflation Factor - VIF)
% Multicollinearity can cause instability in model coefficients.
% Calculate VIF for predictor variables.

vif_vals = zeros(1, length(predictor_vars));  % Preallocate VIF values

% Convert the relevant columns to numeric matrix (predictor variables)
predictor_data = table2array(data_clean(:, predictor_vars));

for i = 1:length(predictor_vars)
    % Create a matrix of all predictor variables except the current one
    predictors_except_i = predictor_data(:, setdiff(1:length(predictor_vars), i));
    
    % Current predictor variable
    current_predictor = predictor_data(:, i);
    
    % Fit a linear model using the current predictor variable as the response
    vif_model = fitlm(predictors_except_i, current_predictor);
    
    % Calculate VIF as 1 / (1 - R^2)
    vif_vals(i) = 1 / (1 - vif_model.Rsquared.Ordinary);
end

% Display VIF values in a table format
vif_table = table(predictor_vars', vif_vals', 'VariableNames', {'Predictor', 'VIF'});
disp(vif_table);

% General rule: VIF > 5 indicates high multicollinearity



%% Step 6: Cross-Validation to avoid overfitting
% Perform 10-fold cross-validation for model validation
cv_model = crossval(models.(response_var), 'KFold', 10);

% Calculate mean squared error for the cross-validated model
mse_cv = kfoldLoss(cv_model);

% Display cross-validated mean squared error
fprintf('Cross-validated MSE for %s: %f\n', response_var, mse_cv);

%% Updated Step 7: Model Interpretation with Proper Titles
% Extract feature names from the model
model = models.(response_var);  % Get the trained model for the response variable
model_feature_names = model.PredictorNames;  % Get the actual feature names

% Plot the effect of each predictor variable using partial dependence plots
figure;
for i = 1:length(model_feature_names)
    subplot(4,4,i);  % Arrange subplots in a grid
    plotPartialDependence(model, model_feature_names{i});
    
    if i==length(model_feature_names)
        datetick('x')
    end
    % Update the title with the correct variable name
    title(predictor_vars{i});  % Use the actual predictor name instead of x1, x2, etc.
    xlabel(predictor_vars{i});
    ylabel('Partial Dependence');
    grid on;
end


%% Calculate permutation importance for predictors
num_permutations = 5;  % Number of permutations for robustness
predictor_importance = zeros(1, length(predictor_vars));  % Preallocate importance array

% Original cross-validated error
original_cv_model = crossval(models.(response_var), 'KFold', 10);
original_mse = kfoldLoss(original_cv_model);

% Loop through each predictor
for i = 1:length(predictor_vars)
    temp_mse = zeros(1, num_permutations);
    
    % Perform multiple permutations to estimate the impact of each predictor
    for p = 1:num_permutations
        permuted_data = predictor_data;  % Copy the predictor data
        permuted_data(:, i) = permuted_data(randperm(size(predictor_data, 1)), i);  % Permute one predictor
        
        % Fit the model with the permuted predictor data
        temp_model = fitrgam(permuted_data, response_data);
        permuted_cv_model = crossval(temp_model, 'KFold', 10);
        temp_mse(p) = kfoldLoss(permuted_cv_model);
    end
    
    % Calculate importance as the difference in MSE after permutation
    predictor_importance(i) = mean(temp_mse) - original_mse;
end

% Display the permutation importance for each predictor
importance_table = table(predictor_vars', predictor_importance', 'VariableNames', {'Predictor', 'Importance'});
disp(importance_table);



% %% Loop through all response taxa
% %% Load Data
% data = physical_forcings_biology_table;
% data.yearmonth = datenum(datetime(data.datetime.Year, data.datetime.Month, 1, 'TimeZone', 'UTC'));
% 
% %% Define predictor and response variables
% predictor_vars = {'sea_water_temperature', 'BEUTI', 'wspd_along', 'AOD_500nm', ...
%                   'surface_downwelling_photosynthetic_photon_flux_in_air', 'pm2_5_sc_pm2_5_daily_slv_tt','Nitrate', ...
%                   'Phosphate', 'Silicate', 'Ammonium','outflow','yearmonth'};
% 
% response_vars = {'Asterionellopsis', 'Centric', 'Skeletonema', 'Hemiaulus', 'Leptocylindrus', 'Thalassionema'};
% 
% predictor_vars_labels = {'Temperature [deg C]', 'BEUTI', 'Wind Speed [m/s]', "Aerosol Optical Depth at 500nm", ...
%                          'PM2.5 [μg/m³]', "PAR [μmol photons m⁻² s⁻¹]", 'Nitrate [μmol/L]', 'Phosphate [μmol/L]', ...
%                          'Silicate [μmol/L]', 'Ammonium [μmol/L]', 'Outflow [m³/s]', 'Year-month'};
% 
% colorz = [0.2, 0.6, 0.8;  % Light blue
%           0.9, 0.5, 0.1;  % Orange
%           0.4, 0.7, 0.2;  % Light green
%           0.8, 0.2, 0.4;  % Pink
%           0.2, 0.4, 0.8;  % Dark blue
%           0.7, 0.2, 0.9;  % Purple
%           0.3, 0.7, 0.6]; % Soft blue-green
% 
% %% Remove rows with missing data in the response variables
% data_clean = data(~any(ismissing(data(:, response_vars)), 2), :);
% data_clean = data_clean(~any(ismissing(data_clean(:, predictor_vars)), 2), :);
% 
% %% Convert the table to numeric arrays
% predictor_data = table2array(data_clean(:, predictor_vars));
% 
% %% Step 1: Fit the GAM model and plot Partial Dependence for each predictor
% models = struct();
% 
% % Set publication settings
% ftsz = 10;  % Font size for publication
% ftname = 'Helvetica';  % Font name
% linewdt = 2;  % Line width for plotting
% figure_size = [0, 0, 7.16, 10];  % Set figure size for publication
% 
% for j = 1:length(predictor_vars)
%     predictor_var = predictor_vars{j};
%     predictor_label = predictor_vars_labels{j};
% 
%     figure('Name', ['Partial Dependence for ', predictor_var], 'NumberTitle', 'off');
%     hold on;
% 
%     % Loop through each response variable to fit the GAM model and plot
%     for i = 1:length(response_vars)
%         response_var = response_vars{i};
% 
%         % Convert current response variable to numeric array
%         response_data = table2array(data_clean(:, response_var));
% 
%         % Fit the model
%         models.(response_var) = fitrgam(predictor_data, response_data);
% 
%         % Plot partial dependence for the current predictor
%         plotPartialDependence(models.(response_var), predictor_var, 'Color', colorz(j, :), 'LineWidth', linewdt);
%     end
% 
%     % Set the title and axis labels
%     title(['Partial Dependence for ', predictor_label]);
%     xlabel(predictor_label);
%     ylabel('Partial Dependence');
%     grid on;
% 
%     % Apply formatting for publication
%     set(gca, 'FontName', ftname, 'FontSize', ftsz);
% 
%     % Set figure size and aspect ratio
%     set(gcf, 'Units', 'inches', 'Position', figure_size);  
%     set(gcf, 'PaperPositionMode', 'auto');  % Ensure the figure fits the paper size
% 
%     % Save the figure
%     saving = 1;
%     if saving == 1
%         png_filename = sprintf('C:\\path_to_directory\\partial_dependence_%s.png', predictor_var);
%         pdf_filename = sprintf('C:\\path_to_directory\\partial_dependence_%s.pdf', predictor_var);
%         print(gcf, '-dpng', '-r300', png_filename);
%         print(gcf, '-dpdf', '-r300', pdf_filename);
%     end
% 
%     hold off;
%     close;  % Close the figure after saving
% end



% %% Load Data
% data = physical_forcings_biology_table;
% data.yearmonth = datenum(datetime(data.datetime.Year, data.datetime.Month, 1, 'TimeZone', 'UTC'));
% 
% %% Define predictor and response variables
% predictor_vars = {'sea_water_temperature', 'BEUTI', 'wspd_along', 'AOD_500nm', ...
%                   'surface_downwelling_photosynthetic_photon_flux_in_air', 'pm2_5_sc_pm2_5_daily_slv_tt','Nitrate', ...
%                   'Phosphate', 'Silicate', 'Ammonium','outflow','yearmonth'};
% 
% response_vars = {'Asterionellopsis', 'Centric', 'Skeletonema', 'Hemiaulus', 'Leptocylindrus', 'Thalassionema'};
% 
% predictor_vars_labels = {'Temperature [deg C]', 'BEUTI', 'Wind Speed [m/s]', "Aerosol Optical Depth at 500nm", ...
%                          'PM2.5 [μg/m³]', "PAR [μmol photons m⁻² s⁻¹]", 'Nitrate [μmol/L]', 'Phosphate [μmol/L]', ...
%                          'Silicate [μmol/L]', 'Ammonium [μmol/L]', 'Outflow [m³/s]', 'Year-month'};
% 
% colorz = [0.2, 0.6, 0.8;  % Light blue
%           0.9, 0.5, 0.1;  % Orange
%           0.4, 0.7, 0.2;  % Light green
%           0.8, 0.2, 0.4;  % Pink
%           0.2, 0.4, 0.8;  % Dark blue
%           0.7, 0.2, 0.9;  % Purple
%           0.3, 0.7, 0.6]; % Soft blue-green
% 
% %% Remove rows with missing data in the response variables
% data_clean = data(~any(ismissing(data(:, response_vars)), 2), :);
% data_clean = data_clean(~any(ismissing(data_clean(:, predictor_vars)), 2), :);
% 
% %% Convert the table to numeric arrays
% predictor_data = table2array(data_clean(:, predictor_vars));
% 
% %% Step 1: Fit the GAM model and plot Partial Dependence for each predictor
% models = struct();
% 
% % Set publication settings
% ftsz = 10;  % Font size for publication
% ftname = 'Helvetica';  % Font name
% linewdt = 2;  % Line width for plotting
% figure_size = [0, 0, 7.16, 10];  % Set figure size for publication
% 
% for j = 1:length(predictor_vars)
%     predictor_var = predictor_vars{j};
%     predictor_label = predictor_vars_labels{j};
% 
%     figure('Name', ['Partial Dependence for ', predictor_var], 'NumberTitle', 'off');
%     hold on;
% 
%     % Loop through each response variable to fit the GAM model and plot
%     for i = 1:length(response_vars)
%         response_var = response_vars{i};
% 
%         % Convert current response variable to numeric array
%         response_data = table2array(data_clean(:, response_var));
% 
%         % Fit the model
%         models.(response_var) = fitrgam(predictor_data, response_data);
% 
%         % Plot partial dependence for the current predictor
%         h = plotPartialDependence(models.(response_var), predictor_var);
% 
%         % Get the current line handle and set the color and line width
%         set(h, 'Color', colorz(j, :), 'LineWidth', linewdt);
%     end
% 
%     % Set the title and axis labels
%     title(['Partial Dependence for ', predictor_label]);
%     xlabel(predictor_label);
%     ylabel('Partial Dependence');
%     grid on;
% 
%     % Apply formatting for publication
%     set(gca, 'FontName', ftname, 'FontSize', ftsz);
% 
%     % Set figure size and aspect ratio
%     set(gcf, 'Units', 'inches', 'Position', figure_size);  
%     set(gcf, 'PaperPositionMode', 'auto');  % Ensure the figure fits the paper size
% 
%     % Save the figure
%     saving = 0;
%     if saving == 1
%         png_filename = sprintf('C:\\path_to_directory\\partial_dependence_%s.png', predictor_var);
%         pdf_filename = sprintf('C:\\path_to_directory\\partial_dependence_%s.pdf', predictor_var);
%         print(gcf, '-dpng', '-r300', png_filename);
%         print(gcf, '-dpdf', '-r300', pdf_filename);
%     end
% 
%     hold off;
%     close;  % Close the figure after saving
% end
% 
% %%
% 
% %% Load Data
% data = physical_forcings_biology_table;
% data.yearmonth = datenum(datetime(data.datetime.Year, data.datetime.Month, 1, 'TimeZone', 'UTC'));
% 
% %% Define predictor and response variables
% predictor_vars = {'sea_water_temperature', 'BEUTI', 'wspd_along', 'AOD_500nm', 'pm2_5_sc_pm2_5_daily_slv_tt', ...
%                   'surface_downwelling_photosynthetic_photon_flux_in_air','Nitrate', ...
%                   'Phosphate', 'Silicate', 'Ammonium','outflow','yearmonth'};
% 
% response_vars = {'Asterionellopsis', 'Centric', 'Skeletonema', 'Hemiaulus', 'Leptocylindrus', 'Thalassionema'};
% 
% predictor_vars_labels = {'Temperature [deg C]', 'BEUTI', 'Wind Speed [m/s]', "Aerosol Optical Depth at 500nm", ...
%                          'PM2.5 [μg/m³]', "PAR [μmol photons m⁻² s⁻¹]", 'Nitrate [μmol/L]', 'Phosphate [μmol/L]', ...
%                          'Silicate [μmol/L]', 'Ammonium [μmol/L]', 'Outflow [m³/s]', 'Year-month'};
% 
% colorz = [0.2, 0.6, 0.8;  % Light blue
%           0.9, 0.5, 0.1;  % Orange
%           0.4, 0.7, 0.2;  % Light green
%           0.8, 0.2, 0.4;  % Pink
%           0.2, 0.4, 0.8;  % Dark blue
%           0.7, 0.2, 0.9;  % Purple
%           0.3, 0.7, 0.6]; % Soft blue-green
% 
% %% Remove rows with missing data in the response variables
% data_clean = data(~any(ismissing(data(:, response_vars)), 2), :);
% data_clean = data_clean(~any(ismissing(data_clean(:, predictor_vars)), 2), :);
% 
% %% Convert the table to numeric arrays (defining predictor_data)
% predictor_data = table2array(data_clean(:, predictor_vars));
% 
% %% Step 1: Fit the GAM model and plot Partial Dependence using Subplots for Each Predictor
% models = struct();
% 
% % Set publication settings
% ftsz = 10;  % Font size for publication
% ftname = 'Helvetica';  % Font name
% linewdt = 2;  % Line width for plotting
% figure_size = [7.16, 10];  % Set figure size for publication (width x height)
% nrows = 3;  % Number of rows in the subplot grid
% ncols = 4;  % Number of columns in the subplot grid (you have 12 predictors)
% 
% % Loop through each response variable to fit the GAM model
% for i = 1:length(response_vars)
%     response_var = response_vars{i};
% 
%     % Convert current response variable to numeric array
%     response_data = table2array(data_clean(:, response_var));
% 
%     % Fit the model using predictor_data and response_data
%     models.(response_var) = fitrgam(predictor_data, response_data);
% end
% 
% %% Plot partial dependence with subplots for each predictor
% for i = 1:length(response_vars)
%     response_var = response_vars{i};
% 
%     % Create a new figure for each response variable
%     figure('Name', ["Partial Dependence for "+ response_var], ...
%            'NumberTitle', 'off', 'Units', 'inches', 'Position', [0, 0, figure_size(1), figure_size(2)]);
% 
%     % Loop through each predictor
%     for j = 1:length(predictor_vars)
%         predictor_var_label = predictor_vars_labels{j};
% 
%         % Create subplot for the current predictor
%         subplot(nrows, ncols, j);
%         hold on;
% 
%         % Use the index directly (since the model names them as x1, x2, etc.)
%         feature_idx = j;  % The index of the predictor in predictor_data
% 
%         % Plot partial dependence for the current predictor by index
%         plotPartialDependence(models.(response_var), feature_idx);
% 
%         % Retrieve the lines in the current plot and set their color and line width
%         lines = findobj(gca, 'Type', 'Line');
%         set(lines, 'Color', colorz(i, :), 'LineWidth', linewdt);  % Use the same color for the response var
% 
%         % Set the title and axis labels
%         title(predictor_var_label);
%         xlabel(predictor_var_label);
%         ylabel('Partial Dependence');
%         grid on;
% 
%         hold off;
%     end
% 
%     % Apply consistent formatting across all subplots
%     set(findall(gcf, '-property', 'FontName'), 'FontName', ftname, 'FontSize', ftsz);
% 
%     % Set figure size and save the figure
%     saving = 0;
%     if saving == 1
%         png_filename = sprintf('C:\\path_to_directory\\partial_dependence_%s.png', response_var);
%         pdf_filename = sprintf('C:\\path_to_directory\\partial_dependence_%s.pdf', response_var);
%         print(gcf, '-dpng', '-r300', png_filename);
%         print(gcf, '-dpdf', '-r300', pdf_filename);
%     end
% 
%     close;  % Close the figure after saving
% end
% 


%%
% %% Try new model
% 
% %% Load the CSV file into MATLAB
% data = readtable('physical_forcings_biology_table.csv');
% 
% %% Define predictor and response variables
% 
% %% Step 1: Define Extreme Event Thresholds
% % Choose thresholds for extreme events based on domain knowledge or percentiles
% threshold_AOD500 = prctile(data.AOD_500nm, 95);  % 95th percentile of AOD500
% threshold_PM25 = prctile(data.pm2_5_sc_pm2_5_daily_slv_tt, 95);  % 95th percentile of PM2.5
% 
% %% Step 2: Create Indicator Variables for Extreme Events
% % Create binary variables indicating whether an extreme event is happening
% data.Extreme_AOD500 = data.AOD_500nm > threshold_AOD500;
% data.Extreme_PM25 = data.pm2_5_sc_pm2_5_daily_slv_tt > threshold_PM25;
% 
% % Combine the two to capture periods when both AOD500 and PM2.5 are extreme
% data.Extreme_Event = data.Extreme_AOD500 | data.Extreme_PM25;
% 
% %% Step 2: Data Normalization (optional based on the scale of data)
% 
% data_clean = rmmissing(data(:, response_vars));  % Remove rows with missing responses
% data_clean = data_clean(~any(ismissing(data_clean), 2), :);  % Ensure no rows with missing predictor variables
% 
% % Keep only the relevant predictor and response columns
% data_clean = data(:, [predictor_vars, response_vars, 'Extreme_Event']);
% 
% % If predictors are on vastly different scales, consider normalizing
% data_clean{:, predictor_vars} = normalize(data_clean{:, predictor_vars});
% 
% %% Step 3: Fit GAM Model with Interaction Terms
% % Include interaction terms between the extreme event indicator and the predictors
% models = struct();
% for i = 1:length(response_vars)
%     response_var = response_vars{i};
% 
%     % Define the formula for GAM, including interaction terms for the extreme event
%     formula = sprintf('%s ~ %s + Extreme_Event', response_var, strjoin(predictor_vars, ' + '));
% 
%     % Fit the model using GAM
%     models.(response_var) = fitrgam(data, formula);
% 
%     % Display model summary
%     fprintf('GAM model for %s:\n', response_var);
%     disp(models.(response_var));
% end
% 
% %% Step 4: Plot the Partial Dependence of Extreme Events
% figure;
% for i = 1:length(response_vars)
%     response_var = response_vars{i};
% 
%     % Plot the partial dependence of the extreme event on the response variable
%     subplot(3,2,i);
%     plotPartialDependence(models.(response_var), 'Extreme_Event');
%     title(['Effect of Extreme Event on ', response_var]);
%     xlabel('Extreme Event (0 = No, 1 = Yes)');
%     ylabel(response_var);
%     grid on;
% end
% 
% %% Step 5: Time Window Focus
% % Create a time window around extreme events (e.g., 1 week before and 1 week after)
% time_window = 60;  % days before and after the extreme event
% data_clean.event_time_window = movmax(data.Extreme_Event, [time_window, time_window]);
% 
% % Fit a GAM model with this focused time window and examine the response
% for i = 1:length(response_vars)
%     response_var = response_vars{i};
% 
%     % Modify formula to focus on time window
%     formula_time_window = sprintf('%s ~ %s + event_time_window', response_var, strjoin(predictor_vars, ' + '));
% 
%     % Fit the model
%     models_time_window.(response_var) = fitrgam(data, formula_time_window);
% 
%     % Display the results
%     fprintf('GAM model (time window) for %s:\n', response_var);
%     disp(models_time_window.(response_var));
% end
% 
% %% Step 6: Visualize Effect of Time Window on Response
% figure;
% for i = 1:length(response_vars)
%     response_var = response_vars{i};
% 
%     % Plot the effect of the event time window on the response variable
%     subplot(3,2,i);
%     plotPartialDependence(models_time_window.(response_var), 'event_time_window');
%     title(['Effect of Event Time Window on ', response_var]);
%     xlabel('Event Time Window');
%     ylabel(response_var);
%     grid on;
% end
% 
% 
% %% Using Lagged PM2.5
% %% Step 1: Create Lagged PM2.5 Variables for Each Response Variable
% % Use the specific lags from the cross-correlation analysis for each taxon
% lags = [7, -1, 6, 6, 5, 9];  % Lags for each taxon
% 
% for i = 1:length(response_vars)
%     lag = lags(i);
%     response_var = response_vars{i};
% 
%     if lag >= 0
%         % Positive lag means PM2.5 leads the response, shift forward
%         data.(['PM25_lag_', response_var]) = [nan(lag, 1); data.pm2_5_sc_pm2_5_daily_slv_tt(1:end-lag)];
%     else
%         % Negative lag means response leads PM2.5, shift backward
%         data.(['PM25_lag_', response_var]) = [data.pm2_5_sc_pm2_5_daily_slv_tt(-lag+1:end); nan(-lag, 1)];
%     end
% end
% 
% %% Step 2: Log Transform Response Variables (optional)
% % Apply log transformation if needed for each response variable
% for i = 1:length(response_vars)
%     response_var = response_vars{i};
%     data.(['log_', response_var]) = log1p(data.(response_var));  % log1p(x) = log(1 + x) to handle zero values
% end
% 
% %% Step 3: Fit GAM Models for Each Taxon with Lagged PM2.5
% % Fit a GAM for each response variable using its corresponding lagged PM2.5 variable
% models = struct();
% 
% for i = 1:length(response_vars)
%     response_var = response_vars{i};
% 
%     % Use the log-transformed response variable
%     log_response_var = ['log_', response_var];
% 
%     % Include lagged PM2.5 for this specific taxon in the model
%     lagged_pm_var = ['PM25_lag_', response_var];
%     extended_predictor_vars = [predictor_vars, lagged_pm_var];
% 
%     % Define the formula for GAM, including interaction terms and lags
%     formula = sprintf('%s ~ %s', log_response_var, strjoin(extended_predictor_vars, ' + '));
% 
%     % Fit the model using GAM
%     models.(response_var) = fitrgam(data, formula);
% 
%     % Display model summary
%     fprintf('GAM model for %s (Log-Transformed):\n', response_var);
%     disp(models.(response_var));
% end
% 
% %% Step 4: Assess Predictor Importance (Partial Dependence Plots)
% figure;
% for i = 1:length(response_vars)
%     response_var = response_vars{i};
%     log_response_var = ['log_', response_var];
% 
%     % Iterate over each predictor and plot its partial dependence
%     for j = 1:length(predictor_vars)
%         predictor = predictor_vars{j};
% 
%         subplot(length(predictor_vars), 1, j);
%         plotPartialDependence(models.(response_var), predictor);
% 
%         % Title the plots with response variable and predictor names
%         title(sprintf('%s vs %s', response_var, predictor));
%         xlabel(predictor);
%         ylabel(log_response_var);
%         grid on;
%     end
% end
% 
% 
% %% Step 5: Assess Importance by Dropping Predictors (One-by-One)
% % Fit reduced models by dropping each predictor, then compare performance
% drop_model_comparisons = struct();
% 
% for i = 1:length(response_vars)
%     response_var = response_vars{i};
%     log_response_var = ['log_', response_var];
% 
%     % Check if log-transformed response variable exists
%     if ~ismember(log_response_var, data.Properties.VariableNames)
%         error('Response variable %s not found in the dataset', log_response_var);
%     end
% 
%     for j = 1:length(predictor_vars)
%         % Drop one predictor at a time
%         reduced_predictors = predictor_vars([1:j-1, j+1:end]);
%         lagged_pm_var = ['PM25_lag_', response_var];
% 
%         % Make sure the lagged PM2.5 variable exists
%         if ~ismember(lagged_pm_var, data.Properties.VariableNames)
%             error('Lagged variable %s not found in the dataset', lagged_pm_var);
%         end
% 
%         % Append the lagged variable to the reduced predictors
%         reduced_predictors = [reduced_predictors, lagged_pm_var];
% 
%         % Check that all reduced predictors exist in the dataset
%         for k = 1:length(reduced_predictors)
%             if ~ismember(reduced_predictors{k}, data.Properties.VariableNames)
%                 error('Variable %s not found in the dataset', reduced_predictors{k});
%             end
%         end
% 
%         % Define the formula for the reduced model
%         reduced_formula = sprintf('%s ~ %s', log_response_var, strjoin(reduced_predictors, ' + '));
% 
%         % Fit the reduced model
%         try
%             reduced_model = fitrgam(data, reduced_formula);
%         catch err
%             fprintf('Error fitting reduced model for %s without %s: %s\n', response_var, predictor_vars{j}, err.message);
%             continue;
%         end
% 
%         % Store reduced model and its performance
%         drop_model_comparisons.(response_var).(['drop_', predictor_vars{j}]) = reduced_model;
% 
%         % Display results for the reduced model
%         fprintf('Reduced GAM model for %s without %s:\n', response_var, predictor_vars{j});
%         disp(reduced_model);
%     end
% end
% 
% %% Step 6: Compare Models Using AIC and Cross-Validation
% 
% % Structure to store AIC and cross-validation performance for each model
% model_comparisons = struct();
% 
% for i = 1:length(response_vars)
%     response_var = response_vars{i};
%     log_response_var = ['log_', response_var];
% 
%     % Full model performance (AIC and cross-validation)
%     fprintf('Evaluating Full Model for %s:\n', response_var);
%     full_model = models.(response_var);
% 
%     % Get AIC for full model
%     try
%         aic_full = aic(full_model);
%     catch
%         aic_full = NaN; % Handle case where AIC cannot be calculated
%     end
% 
%     % Perform 10-fold cross-validation for the full model
%     cv_full = crossval(full_model, 'KFold', 10);
%     mse_full = kfoldLoss(cv_full);
% 
%     % Store full model performance
%     model_comparisons.(response_var).full_model.aic = aic_full;
%     model_comparisons.(response_var).full_model.mse = mse_full;
% 
%     % Now evaluate reduced models (by dropping one predictor at a time)
%     for j = 1:length(predictor_vars)
%         dropped_predictor = predictor_vars{j};
%         reduced_model_name = ['drop_', dropped_predictor];
% 
%         if isfield(drop_model_comparisons.(response_var), reduced_model_name)
%             reduced_model = drop_model_comparisons.(response_var).(reduced_model_name);
% 
%             % Get AIC for reduced model
%             try
%                 aic_reduced = aic(reduced_model);
%             catch
%                 aic_reduced = NaN;
%             end
% 
%             % Perform cross-validation on the reduced model
%             cv_reduced = crossval(reduced_model, 'KFold', 10);
%             mse_reduced = kfoldLoss(cv_reduced);
% 
%             % Store reduced model performance
%             model_comparisons.(response_var).(reduced_model_name).aic = aic_reduced;
%             model_comparisons.(response_var).(reduced_model_name).mse = mse_reduced;
% 
%             % Display results
%             fprintf('Reduced Model for %s without %s:\n', response_var, dropped_predictor);
%             fprintf('AIC: %f, Cross-Validation MSE: %f\n', aic_reduced, mse_reduced);
%         end
%     end
% end
% 
% %% Step 7: Select the Best Model Based on AIC or Cross-Validation MSE
% for i = 1:length(response_vars)
%     response_var = response_vars{i};
% 
%     fprintf('Best Model Selection for %s:\n', response_var);
% 
%     % Get the AIC and MSE for the full model
%     aic_full = model_comparisons.(response_var).full_model.aic;
%     mse_full = model_comparisons.(response_var).full_model.mse;
% 
%     fprintf('Full Model - AIC: %f, MSE: %f\n', aic_full, mse_full);
% 
%     best_aic = aic_full;
%     best_mse = mse_full;
%     best_model = 'full_model';
% 
%     % Compare reduced models
%     for j = 1:length(predictor_vars)
%         dropped_predictor = predictor_vars{j};
%         reduced_model_name = ['drop_', dropped_predictor];
% 
%         if isfield(model_comparisons.(response_var), reduced_model_name)
%             aic_reduced = model_comparisons.(response_var).(reduced_model_name).aic;
%             mse_reduced = model_comparisons.(response_var).(reduced_model_name).mse;
% 
%             fprintf('Reduced Model (dropped %s) - AIC: %f, MSE: %f\n', dropped_predictor, aic_reduced, mse_reduced);
% 
%             % Select the model with the lowest AIC
%             if aic_reduced < best_aic
%                 best_aic = aic_reduced;
%                 best_model = reduced_model_name;
%             end
% 
%             % Also check MSE (you can prioritize AIC or MSE as needed)
%             if mse_reduced < best_mse
%                 best_mse = mse_reduced;
%                 best_model = reduced_model_name;
%             end
%         end
%     end
% 
%     % Display the best model for this response variable
%     fprintf('Best Model for %s: %s (AIC: %f, MSE: %f)\n\n', response_var, best_model, best_aic, best_mse);
% end
% 
% 
% %% Scrap
% 
% 
% % 
% % 
% % % Specify the year to plot
% % i = 2020;
% % 
% % % Group 1: First set of predictor variables
% % group1_vars = {'pm2_5_sc_pm2_5_daily_slv_tt','temperature', 'BEUTI', 'wspd_along', 'outflow'};
% % group1_titles = {'PM2.5 SLV','Sea Surface Temperature', 'BEUTI', 'Wind Speed Along', 'Outflow'};
% % group1_ylabels = {'PM2.5 [μg/m³]','Temperature [deg C]', 'BEUTI', 'Wind Speed [m/s]', 'Outflow [m³/s]'};
% % 
% % % Group 2: Second set of predictor variables
% % group2_vars = {'surface_downwelling_photosynthetic_photon_flux_in_air', 'Nitrate', ...
% %                'Phosphate', 'Silicate', 'Ammonium'};
% % group2_titles = {'Surface PAR', 'Nitrate', 'Phosphate', 'Silicate', 'Ammonium'};
% % group2_ylabels = {'PAR [μmol photons m⁻² s⁻¹]', 'Nitrate [μmol/L]', 'Nitrite [μmol/L]', ...
% %                   'Phosphate [μmol/L]', 'Silicate [μmol/L]', 'Ammonium [μmol/L]'};
% % 
% % %% Plot Group 1 time series: 'temperature', 'BEUTI', 'mei', 'wspd_along', 'outflow'
% % close all;
% % for idx = 1:length(group1_vars)
% %     predictor = group1_vars{idx};
% %     predictor_title = group1_titles{idx};
% %     predictor_ylabel = group1_ylabels{idx};
% % 
% %     row = strcmp(color_table.Variable, predictor);
% %     color_plt = color_table.HighColor{row};  % High (positive) color
% % 
% %     % Filter for the selected year
% %     year_filter = physical_forcings_biology_table.datetime.Year == i;
% % 
% %     % Check if the current predictor exists in the data
% %     if ismember(predictor, physical_forcings_biology_table.Properties.VariableNames)
% %         subplot(5,1, idx);  % Create a subplot (5 rows, 1 column)
% % 
% %         % Plot the predictor with a 14-day moving average
% %         plot(physical_forcings_biology_table.datetime(year_filter), ...
% %              fillmissing(physical_forcings_biology_table.(predictor)(year_filter), 'linear'),'Color',color_plt,'LineWidth',2);
% %         ylabel(predictor_ylabel);
% %         title(predictor_title);
% %         hold on;
% % 
% %         % Add lines for 6 day lag and August Fire Smoke
% %         xline(czu_date, 'Color', 'k', 'LineWidth', linewdt, 'LineStyle', ':'); 
% %         xline(lag_date, 'Color', '#cf572b', 'LineWidth', linewdt, 'LineStyle', ':'); 
% % 
% %         xline(shade_start, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% %         fill([shade_start, shade_start, shade_end, shade_end], ...
% %              [min(ylim) * 1.2, max(ylim) * 1.2, max(ylim) * 1.2, min(ylim) * 1.2], ...
% %              [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4);  % Adjust 'FaceAlpha' for transparency
% %         xline(shade_end, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% % 
% %         xlim([datetime(2020,1,1,'TimeZone','UTC') datetime(2020,12,31,'TimeZone','UTC')]); % Set x-axis limits using datetime
% % 
% %         % Adjust limits to fit data dynamically (auto-scaling)
% %         y_data = movmean(physical_forcings_biology_table.(predictor)(year_filter), 1);
% %         ylim([min(y_data)*1.2 max(y_data)*1.2]);  % Add padding to the data
% % 
% %         hold off;
% %     else
% %         fprintf('Variable %s not found in the dataset.\n', predictor);
% %     end
% % end
% % set(gcf, 'Position', [0 0 1200 1200]);
% % 
% % % Save the figure if required
% % saving = 1;
% % if saving == 1
% %     saveas(gcf, "C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\group1_timeseries.png");
% %     saveas(gcf, "C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\group1_timeseries.png");
% % end
% % 
% % 
% % 
% % %% Plot Group 2 time series: Modify for Group 2 variables
% % figure;
% % close all;
% % for idx = 1:length(group2_vars)
% %     predictor = group2_vars{idx};
% %     predictor_title = group2_titles{idx};
% %     predictor_ylabel = group2_ylabels{idx};
% % 
% %     row = strcmp(color_table.Variable, predictor);
% %     color_plt = color_table.HighColor{row};  % High (positive) color
% % 
% %     % Filter for the selected year
% %     year_filter = physical_forcings_biology_table.datetime.Year == i;
% % 
% %     % Check if the current predictor exists in the data
% %     if ismember(predictor, physical_forcings_biology_table.Properties.VariableNames)
% %         subplot(5,1, idx);  % Create a subplot (5 rows, 1 column)
% % 
% %         % Plot the predictor with a 14-day moving average
% %         plot(physical_forcings_biology_table.datetime(year_filter), ...
% %              fillmissing(physical_forcings_biology_table.(predictor)(year_filter), 'linear'),'Color',color_plt,'LineWidth',2);
% %         ylabel(predictor_ylabel);
% %         title(predictor_title);
% %         hold on;
% % 
% %         % Add lines for 6 day lag and August Fire Smoke
% %         xline(czu_date, 'Color', 'k', 'LineWidth', linewdt, 'LineStyle', ':'); 
% %         xline(lag_date, 'Color', '#cf572b', 'LineWidth', linewdt, 'LineStyle', ':'); 
% %         xline(august_date, 'Color', 'k', 'LineWidth', linewdt, 'LineStyle', ':');
% % 
% %         xline(shade_start, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% %         fill([shade_start, shade_start, shade_end, shade_end], ...
% %              [min(ylim) * 1.2, max(ylim) * 1.2, max(ylim) * 1.2, min(ylim) * 1.2], ...
% %              [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4);  % Adjust 'FaceAlpha' for transparency
% %         xline(shade_end, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% % 
% %         xlim([datetime(2020,1,1,'TimeZone','UTC') datetime(2020,12,31,'TimeZone','UTC')]); % Set x-axis limits using datetime
% % 
% %         % Adjust limits to fit data dynamically (auto-scaling)
% %         y_data = movmean(physical_forcings_biology_table.(predictor)(year_filter), 1);
% %         ylim([min(y_data)*1.2 max(y_data)*1.2]);  % Add padding to the data
% % 
% %         hold off;
% %     else
% %         fprintf('Variable %s not found in the dataset.\n', predictor);
% %     end
% % end
% % set(gcf, 'Position', [0 0 1200 1200]);
% % 
% % % Save the figure if required
% % saving = 1;
% % if saving == 1
% %     saveas(gcf, "C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\group2_timeseries.png");
% %     saveas(gcf, "C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\group2_timeseries.png");
% % end
% % 
% % 
% % 
% % %% Next Anomalies. Group 1
% % 
% % 
% % 
% % close all
% % 
% % % Loop through the predictor variables in group1_vars
% % for idx = 1:length(group1_vars)
% %     predictor = group1_vars{idx};
% %     predictor_title = group1_titles{idx};
% %     predictor_ylabel = group1_ylabels{idx};
% % 
% % 
% %     % Find the color pair for this variable from the color_table
% %     row = strcmp(color_table.Variable, predictor);
% %     hcolor = color_table.HighColor{row};  % High (positive) color
% %     lcolor = color_table.LowColor{row};   % Low (negative) color
% % 
% %     % Check if the current predictor exists in the data
% %     if ismember(predictor, physical_forcings_biology_table.Properties.VariableNames)
% %         subplot(length(group1_vars), 1, idx);  % Create a subplot (5 rows, 1 column)
% % 
% %         % Calculate and plot the anomaly for the current predictor
% %         shade_anomaly_month(physical_forcings_biology_table.datetime, ...
% %                             physical_forcings_biology_table.(predictor), ...
% %                             hcolor, lcolor);
% % 
% %         % Add labels, titles, and other plot elements
% %         ylabel(predictor_ylabel,'FontSize',16);
% %         title(predictor_title);
% %         hold on;
% % 
% %         xlim([datenum(datetime(2020,1,1)) datenum(datetime(2020,12,31))])
% % 
% %         %Case specific limits
% %         if strcmp(predictor, 'pm2_5_sc_pm2_5_daily_slv_tt')
% %             ylim([-20 400])
% %         elseif strcmp(predictor, 'outflow')
% %             ylim([-300 300])
% %         elseif strcmp(predictor, 'temperature')
% %             ylim([-7 7])
% %         end
% % 
% %         datetick('x','keeplimits')
% % 
% % 
% %         % Add key lines and shaded areas
% %         xline(lag_date_datenum, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% %         xline(czu_date_datenum, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');
% %         xline(august_date_datenum, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');
% %         xline(start_date_datenum, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
% %         fill([shade_start_datenum, shade_start_datenum, shade_end_datenum, shade_end_datenum], ...
% %              [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25 min(ylim) * 1.25], ...
% %              [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4); 
% %         xline(end_date_datenum, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% % 
% %         hold off;
% %     else
% %         fprintf('Variable %s not found in the dataset.\n', predictor);
% %     end
% % end
% % % Set figure position
% % set(gcf, 'Position', [0 0 1200 1200]);
% % 
% % % Save the figure if required
% % saving = 1;
% % if saving == 1
% %     saveas(gcf, "C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\group1_anomalies.png");
% %     saveas(gcf, "C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\group1_anomalies.pdf");
% % end
% % 
% % 
% % %% Group2 Anomalies
% % clf
% % % Loop through the predictor variables in group2_vars
% % for idx = 1:length(group2_vars)
% %     predictor = group2_vars{idx};
% %     predictor_title = group2_titles{idx};
% %     predictor_ylabel = group2_ylabels{idx};
% % 
% % 
% %     % Find the color pair for this variable from the color_table
% %     row = strcmp(color_table.Variable, predictor);
% %     hcolor = color_table.HighColor{row};  % High
% % %  (positive) color
% %     lcolor = color_table.LowColor{row};   % Low (negative) color
% % 
% %     % Check if the current predictor exists in the data
% %     if ismember(predictor, physical_forcings_biology_table.Properties.VariableNames)
% %         subplot(length(group2_vars), 1, idx);  % Create a subplot (5 rows, 1 column)
% % 
% %         % Calculate and plot the anomaly for the current predictor
% %         shade_anomaly_month(physical_forcings_biology_table.datetime, ...
% %                             physical_forcings_biology_table.(predictor), ...
% %                             hcolor, lcolor);
% % 
% %         % Add labels, titles, and other plot elements
% %         ylabel(predictor_ylabel,'FontSize',16);
% %         title(predictor_title);
% %         hold on;
% % 
% %         %Case specific limits
% %         if strcmp(predictor, 'surface_downwelling_photosynthetic_photon_flux_in_air')
% %             ylim([-800 800])
% %         elseif strcmp(predictor, 'outflow')
% %             ylim([-300 300])
% %         end
% % 
% %         % Add key lines and shaded areas
% %         xline(lag_date_datenum, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% %         xline(czu_date_datenum, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');
% %         xline(august_date_datenum, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');
% %         xline(start_date_datenum, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
% %         fill([shade_start_datenum, shade_start_datenum, shade_end_datenum, shade_end_datenum], ...
% %              [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25 min(ylim) * 1.25], ...
% %              [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4); 
% %         xline(end_date_datenum, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% %         xlim([datenum(datetime(2020,1,1)) datenum(datetime(2020,12,31))])
% %         datetick('x','keeplimits')
% % 
% %         hold off;
% %     else
% %         fprintf('Variable %s not found in the dataset.\n', predictor);
% %     end
% % end
% % % Set figure position
% % set(gcf, 'Position', [0 0 1200 1200]);
% % 
% % % Save the figure if required
% % saving = 1;
% % if saving == 1
% %     saveas(gcf, "C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\group2_anomalies.png");
% %     saveas(gcf, "C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\group2_anomalies.pdf");
% % end
% % 
% 
% % 
% % %% Visualize each time series for the fire window
% % 
% % close all
% % i=2020;
% % 
% % %SST
% % plot(mlm_sst_tt.datetime(mlm_sst_tt.datetime.Year==i),movmean(mlm_sst_tt.sea_water_temperature(mlm_sst_tt.datetime.Year==i),14))
% % datetick('x')
% % ylabel(['Sea Surface',newline,'Temperature [deg C]']);
% % title('Moss Landing Marine Lab Temperature');
% % hold on
% % 
% % % Add shaded area to the plot
% % % Add lines for 8 day lag and August Fire Smoke
% % xline(lag_date, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% % xline(august_date, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':');
% % xline(start_date, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
% % fill([shade_start, shade_start, shade_end, shade_end], ...
% %     [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25 min(ylim) * 1.25], ...
% %     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4); 
% % xline(end_date, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% % 
% % %% BEUTI
% % close all
% % plot(BEUTI_SC.datetime(BEUTI_SC.datetime.Year==i),movmean(BEUTI_SC.BEUTI(BEUTI_SC.datetime.Year==i),14))
% % datetick('x','keeplimits')
% % ylabel(['BEUTI']);
% % title('BEUTI');
% % hold on
% % 
% % % Add shaded area to the plot
% % % Add lines for 8 day lag and August Fire Smoke
% % xline(lag_date, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% % xline(august_date, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':');
% % xline(start_date, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
% % fill([shade_start, shade_start, shade_end, shade_end], ...
% %     [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25 min(ylim) * 1.25], ...
% %     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4); 
% % xline(end_date, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% % 
% % 
% % % n_plot=n_plot+1;
% % % subplot(n_subplot,1,n_plot)
% % % shade_anomaly(datenum(HABs_SantaCruzWharf.datetime(HABs_SantaCruzWharf.datetime.Year==i)),HABs_SantaCruzWharf.anomaly_Nitrate(HABs_SantaCruzWharf.datetime.Year==i),'#f0c76e','#5bd4c8')
% % % datetick('x')
% % % ylabel(['Nitrate']);
% % % title('Nitrate');
% % % hold on
% % 
% % %% Plot values
% % i=2020;
% % close all; clf
% % figure(2)
% % n_subplot=4;
% % n_plot=0;
% % n_plot=n_plot+1;
% % subplot(n_subplot,1,n_plot)
% % plot(mlm_sst_tt.datetime(mlm_sst_tt.datetime.Year==i),movmean(mlm_sst_tt.sea_water_temperature(mlm_sst_tt.datetime.Year==i),14))
% % ylim([12 18])
% % datetick('x')
% % ylabel(['Sea Surface',newline,'Temperature [deg C]']);
% % title('Moss Landing Marine Lab Temperature');
% % hold on
% % 
% % % Add shaded area to the plot
% % % Add lines for 8 day lag and August Fire Smoke
% % xline(lag_date, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% % xline(august_date, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':');
% % xline(start_date, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
% % fill([shade_start, shade_start, shade_end, shade_end], ...
% %     [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25 min(ylim) * 1.25], ...
% %     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4); 
% % xline(end_date, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% % 
% % 
% % n_plot=n_plot+1;
% % subplot(n_subplot,1,n_plot)
% % plot(habs_tt.datetime(habs_tt.datetime.Year==i),habs_tt.Nitrate(habs_tt.datetime.Year==i),'Color','#5bd4c8')
% % datetick('x')
% % ylabel(['Nitrate']);
% % title('Nitrate');
% % hold on
% % 
% % % Add shaded area to the plot
% % % Add lines for 8 day lag and August Fire Smoke
% % xline(lag_date, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% % xline(august_date, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':');
% % xline(start_date, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
% % fill([shade_start, shade_start, shade_end, shade_end], ...
% %     [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25 min(ylim) * 1.25], ...
% %     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4); 
% % xline(end_date, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% % 
% % 
% % n_plot=n_plot+1;
% % subplot(n_subplot,1,n_plot)
% % plot(habs_tt.datetime(habs_tt.datetime.Year==i),habs_tt.Phosphate(habs_tt.datetime.Year==i),'Color','#afdef0')
% % datetick('x')
% % ylabel(['Phosphate']);
% % title('HABS Phosphate');
% % hold on
% % 
% % % Add shaded area to the plot
% % % Add lines for 8 day lag and August Fire Smoke
% % xline(lag_date, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% % xline(august_date, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':');
% % xline(start_date, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
% % fill([shade_start, shade_start, shade_end, shade_end], ...
% %     [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25 min(ylim) * 1.25], ...
% %     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4); 
% % xline(end_date, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% % 
% % n_plot=n_plot+1;
% % subplot(n_subplot,1,n_plot)
% % plot(habs_tt.datetime(habs_tt.datetime.Year==i),habs_tt.Silicate(habs_tt.datetime.Year==i),'Color','#c7c5c1')
% % datetick('x')
% % ylabel(['Silicate']);
% % title('HABS Silicate');
% % hold on
% % 
% % % Add shaded area to the plot and add lines for 8 day lag and August Fire Smoke
% % xline(lag_date, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% % xline(august_date, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':');
% % xline(start_date, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
% % fill([shade_start, shade_start, shade_end, shade_end], ...
% %     [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25 min(ylim) * 1.25], ...
% %     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4); 
% % xline(end_date, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% % 
% % 
% % % Set figure position
% % set(gcf, 'Position', [0 0 1600 800]);
% % 
% % %% Plot anomalies
% % i=2020;
% % close all; clf
% % figure(2)
% % n_subplot=4;
% % n_plot=0;
% % n_plot=n_plot+1;
% % subplot(n_subplot,1,n_plot)
% % shade_anomaly(datenum(mlm_sst_tt.datetime(mlm_sst_tt.datetime.Year==i)),mlm_sst_tt.anomaly(mlm_sst_tt.datetime.Year==i))
% % datetick('x')
% % ylabel(['Sea Surface',newline,'Temperature [deg C]']);
% % title('Moss Landing Marine Lab Temperature');
% % hold on
% % 
% % % Add shaded area to the plot
% % % Add lines for 8 day lag and August Fire Smoke
% % xline(lag_date_datenum, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% % xline(august_date_datenum, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':');
% % xline(start_date_datenum, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
% % fill([shade_start_datenum, shade_start_datenum, shade_end_datenum, shade_end_datenum], ...
% %     [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25 min(ylim) * 1.25], ...
% %     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4); 
% % xline(end_date_datenum, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% % 
% % 
% % n_plot=n_plot+1;
% % subplot(n_subplot,1,n_plot)
% % shade_anomaly(datenum(HABs_SantaCruzWharf.datetime(HABs_SantaCruzWharf.datetime.Year==i)),HABs_SantaCruzWharf.anomaly_Nitrate(HABs_SantaCruzWharf.datetime.Year==i),'#f0c76e','#5bd4c8')
% % datetick('x')
% % ylabel(['Nitrate']);
% % title('Nitrate');
% % hold on
% % 
% % % Add shaded area to the plot
% % % Add lines for 8 day lag and August Fire Smoke
% % xline(lag_date_datenum, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% % xline(august_date_datenum, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':');
% % xline(start_date_datenum, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
% % fill([shade_start_datenum, shade_start_datenum, shade_end_datenum, shade_end_datenum], ...
% %     [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25 min(ylim) * 1.25], ...
% %     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4); 
% % xline(end_date_datenum, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% % 
% 
% n_plot=n_plot+1;
% subplot(n_subplot,1,n_plot)
% shade_anomaly(datenum(HABs_SantaCruzWharf.datetime(HABs_SantaCruzWharf.datetime.Year==i)),HABs_SantaCruzWharf.anomaly_Phosphate(HABs_SantaCruzWharf.datetime.Year==i),'#ede32b','#afdef0')
% datetick('x')
% ylabel(['Phosphate']);
% title('HABS Phosphate');
% hold on
% 
% % Add shaded area to the plot
% % Add lines for 8 day lag and August Fire Smoke
% xline(lag_date_datenum, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% xline(august_date_datenum, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':');
% xline(start_date_datenum, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
% fill([shade_start_datenum, shade_start_datenum, shade_end_datenum, shade_end_datenum], ...
%     [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25 min(ylim) * 1.25], ...
%     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4); 
% xline(end_date_datenum, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% 
% n_plot=n_plot+1;
% subplot(n_subplot,1,n_plot)
% shade_anomaly(datenum(HABs_SantaCruzWharf.datetime(HABs_SantaCruzWharf.datetime.Year==i)),HABs_SantaCruzWharf.anomaly_Silicate(HABs_SantaCruzWharf.datetime.Year==i),'#0ad125','#c7c5c1')
% datetick('x')
% ylabel(['Silicate']);
% title('HABS Silicate');
% hold on
% 
% % Add shaded area to the plot and add lines for 8 day lag and August Fire Smoke
% xline(lag_date_datenum, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% xline(august_date_datenum, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':');
% xline(start_date_datenum, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
% fill([shade_start_datenum, shade_start_datenum, shade_end_datenum, shade_end_datenum], ...
%     [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25 min(ylim) * 1.25], ...
%     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4); 
% xline(end_date_datenum, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% 
% 
% % Set figure position
% set(gcf, 'Position', [0 0 1600 800]);
% 
% saving=1;
% if saving==1
%     saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\physical_drivers_"+num2str(n_subplot)+"_anomalies_sst_nuts.png"]);
%     saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\physical_drivers_"+num2str(n_subplot)+"_anomalies_sst_nuts.png"]);
% end
% saving=0;
% 
% 
% %% Physical drivers anomalies 2
% 
% n_subplot=4;
% n_plot=0;
% 
% 
% 
% clf
% figure(1)
% i=2020;
% 
% n_plot=n_plot+1;
% subplot(n_subplot,1,n_plot)
% shade_anomaly(datenum(Monterey_lvl_1_5.datetime(Monterey_lvl_1_5.datetime.Year==i)),Monterey_lvl_1_5.anomaly(Monterey_lvl_1_5.datetime.Year==i))
% datetick('x')
% ylim([-0.5 6])
% ylabel(['AOD 500nm']);
% title('Monterey Aeronet Aerosol Optical Depth');
% hold on
% 
% % Add shaded area to the plot and add lines for 8 day lag and August Fire Smoke
% xline(lag_date_datenum, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% xline(august_date_datenum, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':');
% xline(start_date_datenum, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
% fill([shade_start_datenum, shade_start_datenum, shade_end_datenum, shade_end_datenum], ...
%     [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25 min(ylim) * 1.25], ...
%     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4); 
% xline(end_date_datenum, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% 
% 
% n_plot=n_plot+1;
% subplot(n_subplot,1,n_plot)
% shade_anomaly(datenum(BEUTI_SC.datetime(BEUTI_SC.datetime.Year==i)),BEUTI_SC.anomaly(BEUTI_SC.datetime.Year==i),'#347a0d','#b9e3a1')
% datetick('x')
% ylim([-20 40])
% ylabel(['BEUTI at 37N']);
% title('BEUTI');
% hold on
% 
% % Add shaded area to the plot and add lines for 8 day lag and August Fire Smoke
% xline(lag_date_datenum, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% xline(august_date_datenum, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':');
% xline(start_date_datenum, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
% fill([shade_start_datenum, shade_start_datenum, shade_end_datenum, shade_end_datenum], ...
%     [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25 min(ylim) * 1.25], ...
%     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4); 
% xline(end_date_datenum, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% 
% 
% n_plot=n_plot+1;
% subplot(n_subplot,1,n_plot)
% shade_anomaly(datenum(ndbc46042_tt.datetime(ndbc46042_tt.datetime.Year==i)),ndbc46042_tt.anomaly(ndbc46042_tt.datetime.Year==i),'#eba42a','#71d3de')
% datetick('x')
% ylim([-8 8])
% ylabel(['Wind Speed [m/s]']);
% title('NDBC Alongshore Wind');
% hold on
% 
% % Add shaded area to the plot and add lines for 8 day lag and August Fire Smoke
% xline(lag_date_datenum, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% xline(august_date_datenum, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':');
% xline(start_date_datenum, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
% fill([shade_start_datenum, shade_start_datenum, shade_end_datenum, shade_end_datenum], ...
%     [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25 min(ylim) * 1.25], ...
%     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4); 
% xline(end_date_datenum, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% 
% 
% 
% n_plot=n_plot+1;
% subplot(n_subplot,1,n_plot)
% shade_anomaly(datenum(sanlorenzodailyoutflow.datetime(sanlorenzodailyoutflow.datetime.Year==i)),sanlorenzodailyoutflow.anomaly(sanlorenzodailyoutflow.datetime.Year==i),'#1e0db8','#c1bbfc')
% datetick('x')
% ylim([-250 300])
% ylabel(['Outflow [CFS]']);
% title('San Lorenzo River Outflow');
% hold on
% 
% % Add shaded area to the plot and add lines for 8 day lag and August Fire Smoke
% xline(lag_date_datenum, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% xline(august_date_datenum, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':');
% xline(start_date_datenum, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
% fill([shade_start_datenum, shade_start_datenum, shade_end_datenum, shade_end_datenum], ...
%     [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25 min(ylim) * 1.25], ...
%     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4); 
% xline(end_date_datenum, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% 
% % Set figure position
% set(gcf, 'Position', [0 0 1600 800]);
% 
% 
% saving=1;
% if saving==1
%     saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\physical_drivers_"+num2str(n_subplot)+"_anomalies.png"]);
%     saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\physical_drivers_"+num2str(n_subplot)+"_anomalies.pdf"]);
% end
% saving=0;
% 
% 
% %% MLML Surface PAR
% load('mlml_par_tt.mat')
% 
% peak_fire=datetime(2020,8,20);
% lag_days = 8;
% start_date = datetime(2020, 8, 18);
% end_date = datetime(2020, 9, 22);
% lag_date = datetime(2020, 8, 20 + lag_days);
% august_date = datetime(2020, 9, 5);
% peak_fire=datetime(2020,8,20);
% 
% %By day
% start_date_day = day(start_date, 'dayofyear');
% end_date_day = day(end_date, 'dayofyear');
% lag_date_day = day(lag_date, 'dayofyear');
% peak_fire_day=day(peak_fire,'dayofyear');
% august_date_day = day(august_date, 'dayofyear');
% 
% %By datenum
% start_date_datenum=datenum(start_date);
% end_date_datenum=datenum(end_date);
% lag_date_datenum=datenum(lag_date);
% august_date_datenum=datenum(august_date);
% 
% %For shading
% shade_start_day = start_date_day;
% shade_end_day = end_date_day;
% 
% 
% mlmlpar_tt=mlmlpar_tt(mlmlpar_tt.surface_downwelling_photosynthetic_photon_flux_in_air_qc_agg>1,:);
% mlmlpar_tt.dayofyear=day(mlmlpar_tt.datetime,'dayofyear');
% 
% close all
% subplot(2,1,1)
% plot(mlmlpar_tt.datetime(mlmlpar_tt.datetime.Year==2020),...
%     mlmlpar_tt.surface_downwelling_photosynthetic_photon_flux_in_air(mlmlpar_tt.datetime.Year==2020),'LineWidth',2)
% hold on
% % Add shaded area to the plot and add lines for 8 day lag and August Fire Smoke
% xline(peak_fire_day, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':'); 
% xline(lag_date_day, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% xline(august_date_day, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':');
% xline(start_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
% fill([shade_start_day, shade_start_day, shade_end_day, shade_end_day], ...
%     [min(ylim) * 1.5, max(ylim) * 1.5, max(ylim) * 1.5 min(ylim) * 1.5], ...
%     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4); 
% xline(end_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% 
% 
% plotClimatology(mlmlpar_tt,'surface_downwelling_photosynthetic_photon_flux_in_air')
% ylim([0 800])
% xlim([start_date_day-50 end_date_day+50])
% 
% 
% %PM 2.5
% subplot(2,1,2)
% sc_pm2_5_daily_sc_tt.dayofyear=day(sc_pm2_5_daily_sc_tt.datetime,'dayofyear');
% sc_pm2_5_daily_slv_tt.dayofyear=day(sc_pm2_5_daily_slv_tt.datetime,'dayofyear');
% 
% plot(sc_pm2_5_daily_sc_tt.dayofyear(sc_pm2_5_daily_sc_tt.datetime.Year==2020),...
%     sc_pm2_5_daily_sc_tt.pm2_5(sc_pm2_5_daily_sc_tt.datetime.Year==2020),'LineWidth',2)
% hold on
% plot(sc_pm2_5_daily_slv_tt.dayofyear(sc_pm2_5_daily_slv_tt.datetime.Year==2020),...
%     sc_pm2_5_daily_slv_tt.pm2_5(sc_pm2_5_daily_slv_tt.datetime.Year==2020),'LineWidth',2)
% % Add shaded area to the plot and add lines for 8 day lag and August Fire Smoke
% xline(peak_fire_day, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':'); 
% xline(lag_date_day, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% xline(august_date_day, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':');
% xline(start_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
% fill([shade_start_day, shade_start_day, shade_end_day, shade_end_day], ...
%     [min(ylim) * 1.5, max(ylim) * 1.5, max(ylim) * 1.5 min(ylim) * 1.5], ...
%     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4); 
% xline(end_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
% 
% 
% plotClimatology(sc_pm2_5_daily_tt,'pm2_5')
% ylim([0 800])
% xlim([start_date_day-50 end_date_day+50])
