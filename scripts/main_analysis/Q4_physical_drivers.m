%Physical Drivers/Alternative Forcings
addpath()

%% Load the data
load('Q:\Dante\data\MB_Wildfire_Obs\processed_data\joined_physical_drivers\physical_forcings_biology_table_alltime.mat')

addpath('Q:\Dante\Wildfire_Obs\functions')

% Calculate N:P Ratio (Nitrate/Phosphate) and add to table
physical_forcings_biology_table.N_P_ratio = physical_forcings_biology_table.Nitrate ./ physical_forcings_biology_table.Phosphate;

predictor_vars = {'temperature', 'BEUTI', 'wspd_along', 'AOD_500nm', ...
                  'surface_downwelling_photosynthetic_photon_flux_in_air', 'Nitrate', ...
                   'Phosphate', 'Silicate', 'Ammonium', 'outflow', 'pm2_5_sc_pm2_5_daily_slv_tt'};


%% Dates for plotting



% Dates
lag_days = 6;
start_date = datetime(2020, 8, 16);
start_date.TimeZone = 'UTC';
end_date = datetime(2020, 9, 22);
end_date.TimeZone = 'UTC';
lag_date = datetime(2020, 8, 21 + lag_days);
lag_date.TimeZone = 'UTC';

czu_date = datetime(2020, 8, 21);
czu_date.TimeZone = 'UTC';

august_date = datetime(2020, 9, 11);
august_date.TimeZone = 'UTC';
shade_start = datetime(2020,8,16);
shade_start.TimeZone = 'UTC';
shade_end = datetime(2020,9,22);
shade_end.TimeZone = 'UTC';

%By day
start_date_day = day(start_date, 'dayofyear');
end_date_day = day(end_date, 'dayofyear');
lag_date_day = day(lag_date, 'dayofyear');
czu_date_day=day(czu_date,'dayofyear');
august_date_day = day(august_date, 'dayofyear');

%By datenum
start_date_datenum = datenum(start_date);
end_date_datenum = datenum(end_date);
czu_date_datenum=datenum(czu_date);
lag_date_datenum = datenum(lag_date);
august_date_datenum = datenum(august_date);

%For shading
shade_start_datenum = start_date_datenum;
shade_end_datenum = end_date_datenum;

% Define the variables and corresponding high/low colors
variables = {'temperature', 'BEUTI', 'wspd_along','wspd_across','AOD_500nm', 'outflow', ...
             'surface_downwelling_photosynthetic_photon_flux_in_air', 'Nitrate', ...
             'Phosphate', 'Silicate', 'Ammonium', 'pm2_5_sc_pm2_5_daily_slv_tt','N_P_ratio'};

high_colors = {'#de4e6a', '#347a0d', '#eba42a', '#5bb4b3', '#5831b8', ...
               '#bc5090', '#5831b8', '#ff6f61', '#ede32b', '#0ad125', '#ffb6c1','#d62728','#d62728'};

low_colors  = {'#5698e8', '#b9e3a1', '#71d3de', '#f2a950', '#a89ff3', ...
               '#ffa600', '#a89ff3', '#6baed6', '#afdef0', '#c7c5c1', '#4682b4','#2ca02c','#2ca02c'};

% Create a table with the variables and corresponding colors
color_table = table(variables', high_colors', low_colors', ...
                    'VariableNames', {'Variable', 'HighColor', 'LowColor'});

%Plotting Specs
linewdt=2;


%% Group1: Time Series and Anomalies Plotting 

% Define the variables for Group1
group1_vars = {'pm2_5_sc_pm2_5_daily_slv_tt', 'temperature', 'wspd_along','wspd_across', ...
               'surface_downwelling_photosynthetic_photon_flux_in_air', 'Silicate'};
group1_titles = {'PM2.5 SLV', 'Temperature', 'Alongshore Wind Speed','Acrossshore Wind Speed', 'Surface PAR', 'Silicate'};
group1_ylabels = {'PM2.5 [μg ^{-3}]', 'Temperature [deg C]','Wind Speed [m s^{-1}]', 'Wind Speed [m s^{-1}]',"PAR"+newline+"[μmol photons m⁻² s⁻¹]", 'Silicate [μmol L^{-1}]'};

group1_ylabels = [ ...
    "[μg m^{-3}]", ...
    "[°C]", ...
    "[m s^{-1}]" + newline + "N ←→ S", ...
    "[m s^{-1}]" + newline + "W ←→ E", ...
    "[μmol photons" + newline + "m⁻² s⁻¹]", ...
    "[μmol ^{-1}]" ...
];
letters = 'abcdefghijklmnopqrstuvwxyz'; % Letters for labeling





%% Plot Group1 Time Series with Helvetica Font and Adjusted Dimensions for Publication

% Set publication settings
ftsz = 8;  % Font size for publication
ftname = 'Helvetica';  % Font name
linewdt = 2;  % Line width for plotting

figure;

for idx = 1:length(group1_vars)
    predictor = group1_vars{idx};
    predictor_title = group1_titles{idx};
    predictor_ylabel = group1_ylabels{idx};

    % Find the color for this variable
    row = strcmp(color_table.Variable, predictor);
    color_plt = color_table.HighColor{row};  % High (positive) color

    % Filter for the selected year (2020)
    year_filter = physical_forcings_biology_table.datetime.Year == 2020;

    % Check if the current predictor exists in the data
    if ismember(predictor, physical_forcings_biology_table.Properties.VariableNames)
        subplot(length(group1_vars), 1, idx);  % Create a subplot for each variable

        % Plot the predictor with a 14-day moving average
        plot(physical_forcings_biology_table.datetime(year_filter), ...
             fillmissing(physical_forcings_biology_table.(predictor)(year_filter), 'linear'), ...
             'Color', color_plt, 'LineWidth', linewdt);
        ylabel(predictor_ylabel, 'FontSize', ftsz, 'FontName', ftname);
        title(predictor_title, 'FontSize', ftsz, 'FontName', ftname);
        hold on;

        % Add lines for key dates and shaded area
        xline(czu_date, 'Color', 'k', 'LineWidth', linewdt, 'LineStyle', ':');
        xline(lag_date, 'Color', '#cf572b', 'LineWidth', linewdt, 'LineStyle', ':');
        xline(shade_start, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
        fill([shade_start, shade_start, shade_end, shade_end], [min(ylim) * 1.2, max(ylim) * 1.2, max(ylim) * 1.2, min(ylim) * 1.2], ...
             [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4);  % Shaded region
        xline(shade_end, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');

        xlim([datetime(2020, 1, 1, 'TimeZone', 'UTC'), datetime(2020, 12, 31, 'TimeZone', 'UTC')]);  % Set x-axis limits for the full year

        % Adjust y-axis limits dynamically based on data
        y_data = movmean(physical_forcings_biology_table.(predictor)(year_filter), 1);
        ylim([min(y_data) * 1.2, max(y_data) * 1.2]);  % Add padding to the data range

        hold off;

        ylabel(predictor_ylabel, 'FontSize', ftsz);
        xlabel("")
        text(0.02, 0.95, predictor_title, ...
            'Units', 'normalized', ...
            'FontSize', ftsz, ...
            'FontName', ftname, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top');        
        hold on;
        title("")


    else
        fprintf('Variable %s not found in the dataset.\n', predictor);
    end
end

% Adjust figure size for double-column width and appropriate height
set(gcf, 'Units', 'inches', 'Position', [0, 0, 5, 6]);  % 7.16 inches width and 10 inches height for subplots
set(gcf, 'PaperPositionMode', 'auto');  % Ensure the figure fits the paper size

% Apply consistent font settings across the figure
set(findall(gcf, '-property', 'FontName'), 'FontName', ftname, 'FontSize', ftsz);

%% Optional: Save the figure with high resolution for publication
saving = 1;
if saving == 1
    print(gcf, '-dpng', '-r1200', 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\group1_timeseries.png');
    print(gcf, '-dpdf', '-r1200', 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\group1_timeseries.pdf');
end

%%
figure;
for idx = 1:length(group1_vars)
    predictor = group1_vars{idx};
    predictor_title = group1_titles{idx};
    predictor_ylabel = group1_ylabels{idx};

    % Find the color pair for this variable from the color_table
    row = strcmp(color_table.Variable, predictor);
    hcolor = color_table.HighColor{row};  % High (positive) color
    lcolor = color_table.LowColor{row};   % Low (negative) color

    % Check if the current predictor exists in the data
    if ismember(predictor, physical_forcings_biology_table.Properties.VariableNames)
        subplot(length(group1_vars), 1, idx);  % Create a subplot (length(group1_vars) rows, 1 column)

        % Calculate and plot the anomaly for the current predictor
        shade_anomaly_month(physical_forcings_biology_table.datetime, ...
                            physical_forcings_biology_table.(predictor), ...
                            hcolor, lcolor);
        
        % Case-specific limits
        if strcmp(predictor, 'pm2_5_sc_pm2_5_daily_slv_tt')
            ylim([-20 400])
        elseif strcmp(predictor, 'outflow')
            ylim([-300 300])
        elseif strcmp(predictor, 'temperature')
            ylim([-7 7])
        elseif strcmp(predictor, 'surface_downwelling_photosynthetic_photon_flux_in_air')
            ylim([-800 800])
        end

        % Add labels, titles, and other plot elements
        ylabel(predictor_ylabel, 'FontSize', ftsz);
        xlabel("")
        text(0.02, 0.95, predictor_title, ...
            'Units', 'normalized', ...
            'FontSize', ftsz, ...
            'FontName', ftname, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top');        
        hold on;
        title("")

        % Add key lines and shaded areas
        xline(lag_date_datenum, 'Color', '#cf572b', 'LineWidth', 1, 'LineStyle', ':');
        xline(czu_date_datenum, 'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
        xline(start_date_datenum, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
        fill([shade_start_datenum, shade_start_datenum, shade_end_datenum, shade_end_datenum], ...
             [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25, min(ylim) * 1.25], ...
             [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        xline(end_date_datenum, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
        xlim([datenum(datetime(2020,1,1)) datenum(datetime(2020,12,31))]);
        datetick('x', 'keeplimits');
        set(gca, 'FontSize', ftsz);

        % Add subplot label in the upper-left corner
        text(-0.05, 1.25, letters(idx), 'Units', 'normalized', 'FontSize', ftsz + 2, ...
             'FontName', ftname, 'FontWeight', 'bold');

        hold off;
    else
        fprintf('Variable %s not found in the dataset.\n', predictor);
    end
end

% Adjust figure size for double-column width and appropriate height
set(gcf, 'Units', 'inches', 'Position', [0, 0, 5, 6]);  % 7.16 inches width and 10 inches height for subplots
set(gcf, 'PaperPositionMode', 'auto');  % Ensure the figure fits the paper size

% Apply consistent font settings across the figure
set(findall(gcf, '-property', 'FontName'), 'FontName', ftname, 'FontSize', ftsz);


%% Save Group1 Anomalies Plot if required
if saving == 1
    print(gcf, '-dpng', '-r1200', "C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\group1_anomalies.png");
    print(gcf, '-dpng', '-r1200', "C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\group1_anomalies.pdf");
end


%% Group: Time Series and Anomalies Plotting for Group 2 Variables
group2_vars = {'Nitrate', 'Phosphate', 'Ammonium', 'BEUTI', 'outflow', 'N_P_ratio'};
group2_titles = {'Nitrate', 'Phosphate', 'Ammonium', 'BEUTI', 'Outflow', 'N:P Ratio'};
group2_ylabels = {'Nitrate [μmol L^{-1}]', 'Phosphate [μmol L^{-1}]', 'Ammonium [μmol L^{-1}]', 'BEUTI', 'Outflow [m³ s^{-1}]', 'N:P Ratio'};
group2_ylabels = {'[μmol L^{-1}]', '[μmol L^{-1}]', '[μmol L^{-1}]', '', '[m³ s^{-1}]', ''};


% group2_vars = {'Nitrate', 'Phosphate', 'Ammonium', 'BEUTI', 'outflow'};
% group2_titles = {'Nitrate', 'Phosphate', 'Ammonium', 'BEUTI', 'Outflow'};
% group2_ylabels = {'Nitrate [μmolL^{-1}]', 'Phosphate [μmolL^{-1}]', 'Ammonium [μmolL^{-1}]', 'BEUTI', 'Outflow [m³s^{-1}]'};


%% Plot Group 2 Time Series
figure;
for idx = 1:length(group2_vars)
    predictor = group2_vars{idx};
    predictor_title = group2_titles{idx};
    predictor_ylabel = group2_ylabels{idx};

    % Find the color for this variable
    row = strcmp(color_table.Variable, predictor);
    color_plt = color_table.HighColor{row};  % High (positive) color

    % Filter for the selected year
    year_filter = physical_forcings_biology_table.datetime.Year == 2020;

    % Check if the current predictor exists in the data
    if ismember(predictor, physical_forcings_biology_table.Properties.VariableNames)
        subplot(length(group2_vars), 1, idx);  % Create a subplot for each variable

        % Plot the predictor with a 14-day moving average
        plot(physical_forcings_biology_table.datetime(year_filter), ...
             fillmissing(physical_forcings_biology_table.(predictor)(year_filter), 'linear'), 'Color', color_plt, 'LineWidth', 1);
        % Add labels, titles, and other plot elements
        ylabel(predictor_ylabel, 'FontSize', ftsz);
        xlabel("")
        text(0.02, 0.95, predictor_title, ...
            'Units', 'normalized', ...
            'FontSize', ftsz, ...
            'FontName', ftname, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top');        
        hold on;
        title("")

        hold on;

        % Add Redfield ratio line for N:P Ratio subplot
        if strcmp(predictor, 'N_P_ratio')
            redfield_value = 16;
            yline(redfield_value, 'k--', 'LineWidth', 1, 'DisplayName', 'Redfield Ratio (16)');
        end

        % Add lines for key dates and shading
        xline(czu_date, 'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
        xline(lag_date, 'Color', '#cf572b', 'LineWidth', 1.5, 'LineStyle', ':');
        xline(shade_start, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
        fill([shade_start, shade_start, shade_end, shade_end], [min(ylim) * 1.2, max(ylim) * 1.2, max(ylim) * 1.2, min(ylim) * 1.2], ...
             [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4);  % Shaded region
        xline(shade_end, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');

        xlim([datetime(2020,1,1,'TimeZone','UTC') datetime(2020,12,31,'TimeZone','UTC')]);  % Set x-axis limits using datetime

        % Adjust limits to fit data dynamically
        y_data = movmean(physical_forcings_biology_table.(predictor)(year_filter), 1);
        ylim([min(y_data)*1.2 max(y_data)*1.2]);  % Add padding to the data

        hold off;
    else
        fprintf('Variable %s not found in the dataset.\n', predictor);
    end
end

% Adjust figure size for double-column width and appropriate height
set(gcf, 'Units', 'inches', 'Position', [0, 0, 5, 6]);  % 7.16 inches width and 10 inches height for subplots
set(gcf, 'PaperPositionMode', 'auto');  % Ensure the figure fits the paper size

% Apply consistent font settings across the figure
set(findall(gcf, '-property', 'FontName'), 'FontName', ftname, 'FontSize', ftsz);

%% Save Group 2 Time Series Plot if required
saving = 1;
if saving == 1
    print(gcf, '-dpng', '-r1200', "C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\group2_timeseries.png");
    print(gcf, '-dpng', '-r1200', "C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\group2_timeseries.pdf");
end

%% Group 2 Anomalies Plot
% Plot Group 2 Anomalies with Helvetica Font and Publication-Ready Format

% Set publication settings
ftsz = 8;  % Font size for publication
ftname = 'Helvetica';  % Font name
linewdt = 2;  % Line width for plotting

figure;

letters = 'ghijklmnopqrstuvwxyz'; % Starting with 'f'

for idx = 1:length(group2_vars)
    predictor = group2_vars{idx};
    predictor_title = group2_titles{idx};
    predictor_ylabel = group2_ylabels{idx};

    % Find the color pair for this variable from the color_table
    row = strcmp(color_table.Variable, predictor);
    hcolor = color_table.HighColor{row};  % High (positive) color
    lcolor = color_table.LowColor{row};   % Low (negative) color

    % Check if the current predictor exists in the data
    if ismember(predictor, physical_forcings_biology_table.Properties.VariableNames)
        subplot(length(group2_vars), 1, idx);  % Create a subplot for each variable

        % Calculate and plot the anomaly for the current predictor
        shade_anomaly_month(physical_forcings_biology_table.datetime, ...
                            physical_forcings_biology_table.(predictor), ...
                            hcolor, lcolor);
        
        % Add labels, titles, and other plot elements
        ylabel(predictor_ylabel, 'FontSize', ftsz);
        xlabel("")
        text(0.02, 0.95, predictor_title, ...
            'Units', 'normalized', ...
            'FontSize', ftsz, ...
            'FontName', ftname, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top');        
        hold on;
        title("")

        % Set case-specific y-limits for certain predictors
        if strcmp(predictor, 'Nitrate')
            ylim([-10 10])
        elseif strcmp(predictor, 'Phosphate')
            ylim([-1.5 1.5])
        elseif strcmp(predictor, 'Ammonium')
            ylim([-8 8])
        elseif strcmp(predictor, 'outflow')
            ylim([-300 300])
        elseif strcmp(predictor, 'N_P_ratio')
            ylim([-40 40])  % Set reasonable limits for N:P ratio
        end
        % Add key lines and shaded areas
        xline(lag_date_datenum, 'Color', '#cf572b', 'LineWidth', 1, 'LineStyle', ':');
        xline(czu_date_datenum, 'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
        xline(start_date_datenum, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
        fill([shade_start_datenum, shade_start_datenum, shade_end_datenum, shade_end_datenum], ...
             [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25, min(ylim) * 1.25], ...
             [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        xline(end_date_datenum, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
        xlim([datenum(datetime(2020,1,1)) datenum(datetime(2020,12,31))]);
        datetick('x', 'keeplimits');
        set(gca, 'FontSize', ftsz);
        
        % Set font size and style for axes
        set(gca, 'FontSize', ftsz, 'FontName', ftname);

        % Add subplot label in the upper-left corner
        text(-0.05, 1.25, letters(idx), 'Units', 'normalized', 'FontSize', ftsz + 2, ...
             'FontName', ftname, 'FontWeight', 'bold');
        hold off;
    else
        fprintf('Variable %s not found in the dataset.\n', predictor);
    end
end

% Adjust figure size for double-column width and appropriate height
set(gcf, 'Units', 'inches', 'Position', [0, 0, 5, 6]);  % 7.16 inches width and 10 inches height for subplots
set(gcf, 'PaperPositionMode', 'auto');  % Ensure the figure fits the paper size

%% Save the figure if required
saving = 1;
if saving == 1
    print(gcf, '-dpng', '-r1200', 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\group2_anomalies.png');
    print(gcf, '-dpdf', '-r1200', 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\group2_anomalies.pdf');
end




%% Combine Anomaly Plots

%% Combined Plot: Group 1 and Group 2 Anomalies in a 2-Column Layout

% Set publication settings
ftsz = 10;  % Font size for publication
ftname = 'Helvetica';  % Font name
linewdt = 2;  % Line width for plotting

% Create figure for the combined plot
figure;

% Number of variables in each group
n_group1_vars = length(group1_vars);
n_group2_vars = length(group2_vars);

% Plot Group 1 in the first column (left side)
for idx = 1:n_group1_vars
    predictor = group1_vars{idx};
    predictor_title = group1_titles{idx};
    predictor_ylabel = group1_ylabels{idx};

    % Find the color pair for this variable from the color_table
    row = strcmp(color_table.Variable, predictor);
    hcolor = color_table.HighColor{row};  % High (positive) color
    lcolor = color_table.LowColor{row};   % Low (negative) color

    % Check if the current predictor exists in the data
    if ismember(predictor, physical_forcings_biology_table.Properties.VariableNames)
        subplot(max(n_group1_vars, n_group2_vars), 2, (idx - 1) * 2 + 1);  % Create subplot in the left column

        % Calculate and plot the anomaly for the current predictor
        shade_anomaly_month(physical_forcings_biology_table.datetime, ...
                            physical_forcings_biology_table.(predictor), ...
                            hcolor, lcolor);
        
        % Add labels, titles, and other plot elements
        ylabel(predictor_ylabel, 'FontSize', ftsz, 'FontName', ftname);
        title(predictor_title, 'FontSize', ftsz, 'FontName', ftname);
        hold on;

        % Add lines for key dates and shaded areas
        xline(lag_date_datenum, 'Color', '#cf572b', 'LineWidth', linewdt, 'LineStyle', ':');
        xline(czu_date_datenum, 'Color', 'k', 'LineWidth', linewdt, 'LineStyle', ':');
        xline(start_date_datenum, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
        fill([shade_start_datenum, shade_start_datenum, shade_end_datenum, shade_end_datenum], ...
             [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25, min(ylim) * 1.25], ...
             [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4);  % Shaded region
        xline(end_date_datenum, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');

        % Set x-axis limits and format date ticks
        xlim([datenum(datetime(2020,1,1)) datenum(datetime(2020,12,31))]);
        datetick('x', 'keeplimits');

        % Set font size and style for axes
        set(gca, 'FontSize', ftsz, 'FontName', ftname);

        hold off;
    else
        fprintf('Variable %s not found in the dataset.\n', predictor);
    end
end

% Plot Group 2 in the second column (right side)
for idx = 1:n_group2_vars
    predictor = group2_vars{idx};
    predictor_title = group2_titles{idx};
    predictor_ylabel = group2_ylabels{idx};

    % Find the color pair for this variable from the color_table
    row = strcmp(color_table.Variable, predictor);
    hcolor = color_table.HighColor{row};  % High (positive) color
    lcolor = color_table.LowColor{row};   % Low (negative) color

    % Check if the current predictor exists in the data
    if ismember(predictor, physical_forcings_biology_table.Properties.VariableNames)
        subplot(max(n_group1_vars, n_group2_vars), 2, idx * 2);  % Create subplot in the right column

        % Calculate and plot the anomaly for the current predictor
        shade_anomaly_month(physical_forcings_biology_table.datetime, ...
                            physical_forcings_biology_table.(predictor), ...
                            hcolor, lcolor);
        
        % Add labels, titles, and other plot elements
        ylabel(predictor_ylabel, 'FontSize', ftsz, 'FontName', ftname);
        title(predictor_title, 'FontSize', ftsz, 'FontName', ftname);
        hold on;

        % Set case-specific y-limits for certain predictors
        if strcmp(predictor, 'Nitrate')
            ylim([-10 10])
        elseif strcmp(predictor, 'Phosphate')
            ylim([-1.5 1.5])
        elseif strcmp(predictor, 'Ammonium')
            ylim([-8 8])
        elseif strcmp(predictor, 'outflow')
            ylim([-300 300])
        end

        % Add lines for key dates and shaded areas
        xline(lag_date_datenum, 'Color', '#cf572b', 'LineWidth', linewdt, 'LineStyle', ':');
        xline(czu_date_datenum, 'Color', 'k', 'LineWidth', linewdt, 'LineStyle', ':');
        xline(start_date_datenum, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
        fill([shade_start_datenum, shade_start_datenum, shade_end_datenum, shade_end_datenum], ...
             [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25, min(ylim) * 1.25], ...
             [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4);  % Shaded region
        xline(end_date_datenum, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
        
        % Set x-axis limits and format date ticks
        xlim([datenum(datetime(2020,1,1)) datenum(datetime(2020,12,31))]);
        datetick('x', 'keeplimits');

        % Set font size and style for axes
        set(gca, 'FontSize', ftsz, 'FontName', ftname);

        hold off;
    else
        fprintf('Variable %s not found in the dataset.\n', predictor);
    end
end

% Adjust figure size for double-column width and appropriate height
set(gcf, 'Units', 'inches', 'Position', [0, 0, 7.16, 10]);  % 7.16 inches width and 10 inches height for two columns
set(gcf, 'PaperPositionMode', 'auto');  % Ensure the figure fits the paper size

% Apply consistent font settings across the figure
set(findall(gcf, '-property', 'FontName'), 'FontName', ftname, 'FontSize', ftsz);

% Save the combined figure if required
saving = 1;
if saving == 1
    print(gcf, '-dpng', '-r1200', 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\combined_anomalies_plot.png');
    print(gcf, '-dpdf', '-r1200', 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\combined_anomalies_plot.pdf');
end





%% 2018

% Dates
lag_days = 6;
start_date = datetime(2020, 8, 16);
start_date.TimeZone = 'UTC';
end_date = datetime(2020, 9, 22);
end_date.TimeZone = 'UTC';
lag_date = datetime(2020, 8, 21 + lag_days);
lag_date.TimeZone = 'UTC';

czu_date = datetime(2020, 8, 21);
czu_date.TimeZone = 'UTC';

august_date = datetime(2020, 9, 11);
august_date.TimeZone = 'UTC';
shade_start = datetime(2020,8,16);
shade_start.TimeZone = 'UTC';
shade_end = datetime(2020,9,22);
shade_end.TimeZone = 'UTC';

%By day
start_date_day = day(start_date, 'dayofyear');
end_date_day = day(end_date, 'dayofyear');
lag_date_day = day(lag_date, 'dayofyear');
czu_date_day=day(czu_date,'dayofyear');
august_date_day = day(august_date, 'dayofyear');

%By datenum
start_date_datenum = datenum(start_date);
end_date_datenum = datenum(end_date);
czu_date_datenum=datenum(czu_date);
lag_date_datenum = datenum(lag_date);
august_date_datenum = datenum(august_date);

%% Group1: Time Series and Anomalies Plotting 
ftsz=8
% Define the variables for Group1
group1_vars = {'pm2_5_sc_pm2_5_daily_slv_tt', 'temperature', 'wspd_along','wspd_across', ...
               'surface_downwelling_photosynthetic_photon_flux_in_air', 'Silicate'};
group1_titles = {'PM2.5 SLV', 'Temperature', 'Alongshore Wind Speed','Acrossshore Wind Speed', 'Surface PAR', 'Silicate'};
group1_ylabels = {'PM2.5 [μg ^{-3}]', 'Temperature [deg C]','Wind Speed [m s^{-1}]', 'Wind Speed [m s^{-1}]',"PAR"+newline+"[μmol photons m⁻² s⁻¹]", 'Silicate [μmol L^{-1}]'};

group1_ylabels = [ ...
    "[μg m^{-3}]", ...
    "[°C]", ...
    "[m s^{-1}]" + newline + "N ←→ S", ...
    "[m s^{-1}]" + newline + "W ←→ E", ...
    "[μmol photons" + newline + "m⁻² s⁻¹]", ...
    "[μmol ^{-1}]" ...
];
letters = 'abcdefghijklmnopqrstuvwxyz'; % Letters for labeling



%%
figure;
for idx = 1:length(group1_vars)
    predictor = group1_vars{idx};
    predictor_title = group1_titles{idx};
    predictor_ylabel = group1_ylabels{idx};

    % Find the color pair for this variable from the color_table
    row = strcmp(color_table.Variable, predictor);
    hcolor = color_table.HighColor{row};  % High (positive) color
    lcolor = color_table.LowColor{row};   % Low (negative) color

    % Check if the current predictor exists in the data
    if ismember(predictor, physical_forcings_biology_table.Properties.VariableNames)
        subplot(length(group1_vars), 1, idx);  % Create a subplot (length(group1_vars) rows, 1 column)

        % Calculate and plot the anomaly for the current predictor
        shade_anomaly_month(physical_forcings_biology_table.datetime, ...
                            physical_forcings_biology_table.(predictor), ...
                            hcolor, lcolor);
        
        % Case-specific limits
        if strcmp(predictor, 'pm2_5_sc_pm2_5_daily_slv_tt')
            ylim([-20 400])
        elseif strcmp(predictor, 'outflow')
            ylim([-300 300])
        elseif strcmp(predictor, 'temperature')
            ylim([-7 7])
        elseif strcmp(predictor, 'surface_downwelling_photosynthetic_photon_flux_in_air')
            ylim([-800 800])
        end

        % Add labels, titles, and other plot elements
        ylabel(predictor_ylabel, 'FontSize', ftsz);
        xlabel("")
        text(0.02, 0.95, predictor_title, ...
            'Units', 'normalized', ...
            'FontSize', ftsz, ...
            'FontName', ftname, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top');        
        hold on;
        title("")

        % Add key lines and shaded areas
        xline(lag_date_datenum, 'Color', '#cf572b', 'LineWidth', 1, 'LineStyle', ':');
        xline(czu_date_datenum, 'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
        xline(start_date_datenum, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
        fill([shade_start_datenum, shade_start_datenum, shade_end_datenum, shade_end_datenum], ...
             [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25, min(ylim) * 1.25], ...
             [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        xline(end_date_datenum, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
        xlim([datenum(datetime(2018,1,1)) datenum(datetime(2018,12,31))]);
        datetick('x', 'keeplimits');
        set(gca, 'FontSize', ftsz);

        % Add subplot label in the upper-left corner
        text(-0.05, 1.25, letters(idx), 'Units', 'normalized', 'FontSize', ftsz + 2, ...
             'FontName', ftname, 'FontWeight', 'bold');

        hold off;
    else
        fprintf('Variable %s not found in the dataset.\n', predictor);
    end
end

% Adjust figure size for double-column width and appropriate height
set(gcf, 'Units', 'inches', 'Position', [0, 0, 5, 6]);  % 7.16 inches width and 10 inches height for subplots
set(gcf, 'PaperPositionMode', 'auto');  % Ensure the figure fits the paper size

% Apply consistent font settings across the figure
set(findall(gcf, '-property', 'FontName'), 'FontName', ftname, 'FontSize', ftsz);


%% Save Group1 Anomalies Plot if required
if saving == 1
    print(gcf, '-dpng', '-r1200', "C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\group1_anomalies_2018.png");
    print(gcf, '-dpng', '-r1200', "C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\group1_anomalies_2018.pdf");
end


%% Group: Time Series and Anomalies Plotting for Group 2 Variables
group2_vars = {'Nitrate', 'Phosphate', 'Ammonium', 'BEUTI', 'outflow', 'N_P_ratio'};
group2_titles = {'Nitrate', 'Phosphate', 'Ammonium', 'BEUTI', 'Outflow', 'N:P Ratio'};
group2_ylabels = {'Nitrate [μmol L^{-1}]', 'Phosphate [μmol L^{-1}]', 'Ammonium [μmol L^{-1}]', 'BEUTI', 'Outflow [m³ s^{-1}]', 'N:P Ratio'};
group2_ylabels = {'[μmol L^{-1}]', '[μmol L^{-1}]', '[μmol L^{-1}]', '', '[m³ s^{-1}]', ''};


% group2_vars = {'Nitrate', 'Phosphate', 'Ammonium', 'BEUTI', 'outflow'};
% group2_titles = {'Nitrate', 'Phosphate', 'Ammonium', 'BEUTI', 'Outflow'};
% group2_ylabels = {'Nitrate [μmolL^{-1}]', 'Phosphate [μmolL^{-1}]', 'Ammonium [μmolL^{-1}]', 'BEUTI', 'Outflow [m³s^{-1}]'};



%% Group 2 Anomalies Plot
% Plot Group 2 Anomalies with Helvetica Font and Publication-Ready Format

% Set publication settings
ftsz = 8;  % Font size for publication
ftname = 'Helvetica';  % Font name
linewdt = 2;  % Line width for plotting

figure;

letters = 'ghijklmnopqrstuvwxyz'; % Starting with 'f'

for idx = 1:length(group2_vars)
    predictor = group2_vars{idx};
    predictor_title = group2_titles{idx};
    predictor_ylabel = group2_ylabels{idx};

    % Find the color pair for this variable from the color_table
    row = strcmp(color_table.Variable, predictor);
    hcolor = color_table.HighColor{row};  % High (positive) color
    lcolor = color_table.LowColor{row};   % Low (negative) color

    % Check if the current predictor exists in the data
    if ismember(predictor, physical_forcings_biology_table.Properties.VariableNames)
        subplot(length(group2_vars), 1, idx);  % Create a subplot for each variable

        % Calculate and plot the anomaly for the current predictor
        shade_anomaly_month(physical_forcings_biology_table.datetime, ...
                            physical_forcings_biology_table.(predictor), ...
                            hcolor, lcolor);
        
        % Add labels, titles, and other plot elements
        ylabel(predictor_ylabel, 'FontSize', ftsz);
        xlabel("")
        text(0.02, 0.95, predictor_title, ...
            'Units', 'normalized', ...
            'FontSize', ftsz, ...
            'FontName', ftname, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top');        
        hold on;
        title("")

        % Set case-specific y-limits for certain predictors
        if strcmp(predictor, 'Nitrate')
            ylim([-10 10])
        elseif strcmp(predictor, 'Phosphate')
            ylim([-1.5 1.5])
        elseif strcmp(predictor, 'Ammonium')
            ylim([-8 8])
        elseif strcmp(predictor, 'outflow')
            ylim([-300 300])
        elseif strcmp(predictor, 'N_P_ratio')
            ylim([-40 40])  % Set reasonable limits for N:P ratio
        end
        % Add key lines and shaded areas
        xline(lag_date_datenum, 'Color', '#cf572b', 'LineWidth', 1, 'LineStyle', ':');
        xline(czu_date_datenum, 'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
        xline(start_date_datenum, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
        fill([shade_start_datenum, shade_start_datenum, shade_end_datenum, shade_end_datenum], ...
             [min(ylim) * 1.25, max(ylim) * 1.25, max(ylim) * 1.25, min(ylim) * 1.25], ...
             [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        xline(end_date_datenum, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
        xlim([datenum(datetime(2018,1,1)) datenum(datetime(2018,12,31))]);
        datetick('x', 'keeplimits');
        set(gca, 'FontSize', ftsz);
        
        % Set font size and style for axes
        set(gca, 'FontSize', ftsz, 'FontName', ftname);

        % Add subplot label in the upper-left corner
        text(-0.05, 1.25, letters(idx), 'Units', 'normalized', 'FontSize', ftsz + 2, ...
             'FontName', ftname, 'FontWeight', 'bold');
        hold off;
    else
        fprintf('Variable %s not found in the dataset.\n', predictor);
    end
end

% Adjust figure size for double-column width and appropriate height
set(gcf, 'Units', 'inches', 'Position', [0, 0, 5, 6]);  % 7.16 inches width and 10 inches height for subplots
set(gcf, 'PaperPositionMode', 'auto');  % Ensure the figure fits the paper size

%% Save the figure if required
saving = 1;
if saving == 1
    print(gcf, '-dpng', '-r1200', 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\group2_anomalies_2018.png');
    print(gcf, '-dpdf', '-r1200', 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\anomalies\group2_anomalies_2018.pdf');
end









%% XCORR Each variable with these taxa

load('Q:\Dante\data\MB_Wildfire_Obs\processed_data\joined_physical_drivers\physical_forcings_biology_table_alltime.mat')

%% Define predictor and taxa variables
predictor_vars = {'pm2_5_sc_pm2_5_daily_slv_tt', 'temperature', 'wspd_along','wspd_across', ...
               'surface_downwelling_photosynthetic_photon_flux_in_air', 'Silicate', ...
                  'Nitrate', 'Phosphate', 'Ammonium', 'BEUTI', 'outflow'};
predictor_labels= {'PM2.5 SLV', 'Temperature', 'Along wind speed','Across wind wpeed', 'Surface PAR', 'Silicate',...
    'Nitrate', 'Phosphate', 'Ammonium', 'BEUTI', 'Outflow'};
response_taxa_vars = physical_forcings_biology_table.Properties.VariableNames([3:53,59]);

% Parameters
maxLag = 30;          % Maximum lag
minLag = -10;         % Minimum lag for the window of interest
numtests = 5000;      % Number of surrogate tests for significance testing
alpha = 0.05;         % Significance level for dynamic threshold

% Initialize containers for results for all taxa and predictors
numPredictors = length(predictor_vars);
numTaxa = length(response_taxa_vars);
all_max_correlations = nan(numTaxa, numPredictors);  % Max cross-correlation values
all_max_lags = nan(numTaxa, numPredictors);          % Corresponding lags for max CC
all_p_values = nan(numTaxa, numPredictors);          % Raw p-values
all_adjusted_p_values = nan(numTaxa, numPredictors); % Adjusted p-values after FDR

% Loop through each predictor and each taxon
for i = 1:numPredictors
    predictorData = physical_forcings_biology_table.(predictor_vars{i});  % Extract predictor data
    
    for j = 1:numTaxa
        % Extract taxon data
        taxaData = (physical_forcings_biology_table.(response_taxa_vars{j}));

        % Remove NaNs for xcorr calculation
        valid_idx = ~isnan(predictorData) & ~isnan(taxaData);
        predictorDataClean = predictorData(valid_idx);
        taxaDataClean = taxaData(valid_idx);
        
        % Compute cross-correlation and filter for lag window of interest
        [cc, lags] = xcorr(taxaDataClean, predictorDataClean, maxLag, 'coeff');
        lag_filter = (lags >= minLag) & (lags <= maxLag);
        filtered_cc = cc(lag_filter);
        filtered_lags = lags(lag_filter);
        
        % Find the maximum cross-correlation and corresponding lag
        [max_cc, max_idx] = max(abs(filtered_cc));
        max_lag = filtered_lags(max_idx);
        all_max_correlations(j, i) = max_cc;
        all_max_lags(j, i) = max_lag;
        
        % Permutation test using Fourier surrogate method to establish null distribution
        max_rand_ccs = zeros(numtests, 1);
        for test = 1:numtests
            rand_series1 = generate_fourier_surrogate(taxaDataClean);  % Generate Fourier surrogate
            
            % Compute cross-correlation for surrogate data
            [xc_rand, ~] = xcorr(rand_series1, predictorDataClean, maxLag, 'coeff');
            rand_filtered_cc = xc_rand(lag_filter);
            
            % Store the maximum correlation for surrogate data
            max_rand_ccs(test) = max(abs(rand_filtered_cc));
        end

        % Dynamic threshold and p-value based on null distribution
        dynamic_threshold = prctile(max_rand_ccs, 95);
        p_value = mean(max_rand_ccs >= max_cc);
        all_p_values(j, i) = max(p_value, 1 / numtests);  % Avoid exact zero p-values
    end
end

% Benjamini-Hochberg FDR correction for p-values across all taxa and predictors
all_p_values_vector = all_p_values(:);  % Flatten matrix to vector for FDR processing
[sorted_p_values, sort_idx] = sort(all_p_values_vector, 'ascend');
num_tests = length(all_p_values_vector);

% Calculate FDR-adjusted p-values
adjusted_p_values_vector = sorted_p_values .* num_tests ./ (1:num_tests)';
adjusted_p_values_vector(adjusted_p_values_vector > 1) = 1;

% Re-map the adjusted p-values back to their original positions
unsorted_adjusted_p_values = nan(size(all_p_values_vector));
unsorted_adjusted_p_values(sort_idx) = adjusted_p_values_vector;
all_adjusted_p_values = reshape(unsorted_adjusted_p_values, numTaxa, numPredictors);


%% Plot the heatmap for local maximum cross-correlations
% Define font properties and figure dimensions for publication
ftsz = 8;              % Font size for axis labels and ticks
ftname = 'Helvetica';  % Font name for a professional look

italic_taxa = {
    'Akashiwo', 'Asterionellopsis','Boreadinium', 'Ceratium','Chaetoceros', ...
    'Cochlodinium', 'Corethron', 'Cyl Nitz', ...
    'Det Cer Lau', 'Dictyocha', 'Dinophysis', ...
    'Ditylum', 'Entomoneis', 'Eucampia', ...
    'Gymnodinium','Gyrodinium', 'Hemiaulus', 'Leptocylindrus', ...
    'Peridinium','Phaeocystis', 'Polykrikos','Protocentrum','Protoperidinium', 'Pseudo nitzschia', ...
    'Pyramimonas','Skeletonema', 'Thalassionema', 'Thalassiosira', ...
    'Tiarina', 'Tontonia', 'Torodinium', ...
    'Tropidoneis', 'Vicicitus'
};

% Plot the heatmap for local maximum cross-correlations
figure;
imagesc(all_max_correlations);
colormap(jet);
xlabel('Predictor Variables');
ylabel('');

% Customize X and Y tick labels
xticks(1:numPredictors);
xticklabels(predictor_labels);  % Use custom labels for predictors
yticks(1:numTaxa);

% Taxa label formatting
taxa_clean = strrep(taxa, '_', ' ');
set(gca, 'YTick', 1:numel(taxa_clean));
set(gca, 'YTickLabel', taxa_clean);

% Italicize only taxa in italic_taxa
set(gca, 'TickLabelInterpreter', 'tex');  % Switch from LaTeX back to default interpreter
yticklabels = taxa_clean;

for i = 1:numel(taxa_clean)
    base_name = strrep(taxa_clean{i}, '*', '');
    
    if ismember(base_name, italic_taxa)
        yticklabels{i} = ['\it ', taxa_clean{i}];  % TeX-style italics
    end
end

% Update X and Y axis font properties
ax = gca;  % Get current axis
ax.XAxis.TickLabelRotation = 45;   % Rotate x-axis labels for readability
ax.XAxis.FontSize = ftsz;          % Set x-axis font size
ax.XAxis.FontName = ftname;        % Set x-axis font name
set(gca, 'YTickLabel', yticklabels);
ax.YAxis.FontSize = ftsz;          % Set y-axis font size
ax.YAxis.FontName = ftname;        % Set y-axis font name

cb = colorbar('Location', 'eastoutside');
cb.FontSize = ftsz;
cb.FontName = ftname;
cb.Label.Rotation = 270;         % Keep horizontal
cb.Label.String = '\it Cross correlation';  % Bold italic C
cb.Label.FontSize = ftsz;
cb.Label.FontName = ftname;

% Adjust figure size for publication (double-column width for Limnology and Oceanography)
set(gcf, 'Units', 'inches', 'Position', [0, 0, 5, 6]);  % Width: 7.16 inches, height to fit
set(gcf, 'PaperPositionMode', 'auto');  % Ensure the figure fits the paper size when printed

% Overlay the lag with an asterisk for significant cross-correlations
hold on;
for i = 1:numPredictors
    for j = 1:numTaxa
        if all_adjusted_p_values(j, i) < alpha
            % Display lag with an asterisk for significant correlations
            text(i, j, sprintf('*Lag:%d', all_max_lags(j, i)), ...
                'HorizontalAlignment', 'center', 'Color', 'k', 'FontSize', ftsz-1, 'FontWeight','bold');
        end
    end
end
hold off;


% Save the figure if required
saving = 1;
if saving == 1
    print(gcf,'-dpng', '-r1200', "C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\xcorr\xcorr_heatmap_physical_vars_with_significant.png");
    print(gcf,'-dpdf', '-r1200', "C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\xcorr\xcorr_heatmap_physical_vars_with_significant.pdf");
end

%% Compare SST measurements for 2020
% Load the data
data = readtable('physical_forcings_biology_table_alltime.csv');

%% Convert 'datetime' to a MATLAB datetime format
data.datetime = datetime(data.datetime, 'InputFormat', 'dd-MMM-yyyy');
data.weekofyear = week(data.datetime, "weekofyear");

% Filter data for the year 2020
data_2020 = data(year(data.datetime) == 2020, {'datetime', 'sst', 'sea_water_temperature', 'temperature', 'weekofyear'});

% Retime for weekly data
data_weekly = retime(table2timetable(data), 'weekly', 'median');

% Calculate temperature difference for all years
data_weekly.temp_diff = data_weekly.sst - data_weekly.sea_water_temperature;

% Set publication settings
ftsz = 8;  % Font size for publication
ftname = 'Helvetica';  % Font name
linewdt = 2;  % Line width for plotting

% Adjust figure size for double-column width and appropriate height
figure('Units', 'inches', 'Position', [0, 0, 5, 6]);  % Increased width for outside legends
set(gcf, 'PaperPositionMode', 'auto');

% Use tiled layout for better control over subplot spacing
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% First subplot: Original temperature plot for 2020
nexttile;
plot(data_2020.datetime, data_2020.sst, 'Color', '#e466e8', 'LineWidth', linewdt);  % SST - muted orange
hold on;
plot(data_2020.datetime, data_2020.sea_water_temperature, 'Color', [0.64, 0.08, 0.18], 'LineWidth', linewdt);  % Sea Water Temp - dark red
plot(data_2020.datetime, data_2020.temperature, 'Color', [0.93, 0.69, 0.13], 'LineWidth', linewdt);  % Temp - soft orange

% Customize the first subplot appearance
xlabel('Date', 'FontSize', ftsz, 'FontName', ftname);
ylabel('Temperature (°C)', 'FontSize', ftsz, 'FontName', ftname);
set(gca, 'FontSize', ftsz, 'FontName', ftname);
title('Temperature', 'FontSize', ftsz, 'FontName', ftname);
legend({'Satellite SST', 'MLML (16m)', 'SCMW'}, 'FontSize', ftsz, 'Location', 'eastoutside');
grid on;

% Set x-ticks to the first day of each month in 2020
months_2020 = datetime(2020, 1:12, 1);
set(gca, 'XTick', months_2020);
set(gca, 'XTickLabel', datestr(months_2020, 'mmm'));

hold off;

% Second subplot: Temperature difference with climatology and 2020 data
nexttile;
% Step 1: Calculate climatology data for temperature difference
climatology_out = climatology_se_week(data_weekly, 'temp_diff');

% Step 2: Apply smoothing to the climatology mean and SE, excluding NaN values
smooth_weeks = 1;  % Define the smoothing window size
mean_smooth = movmean(climatology_out.Mean, smooth_weeks, 'omitnan');              % Smoothed mean excluding NaN
se_smooth = abs(movmean(climatology_out.StandardError, smooth_weeks, 'omitnan'));       % Smoothed standard error excluding NaN

% Calculate the upper and lower bounds of the 95% confidence interval
ci95_smooth_upper = mean_smooth + 1.96 * se_smooth;                     % Upper bound of 95% confidence interval
ci95_smooth_lower = mean_smooth - 1.96 * se_smooth;                     % Lower bound of 95% confidence interval

% Step 3: Extract weeks for plotting
weeks = climatology_out.weekOfYear;

% Remove NaN values for plotting
validIndices = ~isnan(mean_smooth) & ~isnan(ci95_smooth_upper) & ~isnan(ci95_smooth_lower);
weeks = weeks(validIndices);
mean_smooth = mean_smooth(validIndices);
ci95_smooth_upper = ci95_smooth_upper(validIndices);
ci95_smooth_lower = ci95_smooth_lower(validIndices);

% Step 4: Plot the smoothed climatology mean
plot(weeks, mean_smooth, '-', 'LineWidth', 1, 'Color', [0.59, 0.61, 0.65]);  % RGB for #979da6
hold on;

% Step 5: Plot the confidence interval using area for shading
% Plot the area between the upper and lower confidence interval bounds
fill([weeks; flipud(weeks)], [ci95_smooth_upper; flipud(ci95_smooth_lower)], ...
     [0.59, 0.61, 0.65], 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

hold on
% Plot the 2020 temperature difference data on top of the climatology
plot(data_weekly.weekofyear(data_weekly.datetime.Year == 2020), ...
     data_weekly.sst(data_weekly.datetime.Year == 2020) - data_weekly.sea_water_temperature(data_weekly.datetime.Year == 2020), ...
     'Color', '#1f77b4', 'LineWidth', linewdt);  % 2020 temperature difference in blue

% Add a horizontal line at 0 for reference
yline(0, 'k--', 'LineWidth', 1);

% Customize the second subplot appearance
xlabel('Date', 'FontSize', ftsz, 'FontName', ftname);
ylabel('Temperature Difference (°C)', 'FontSize', ftsz, 'FontName', ftname);
set(gca, 'FontSize', ftsz, 'FontName', ftname);
title('Temperature Difference (Satellite SST-MLML Temperature)', 'FontSize', ftsz, 'FontName', ftname);
legend({'Climatology', '2020'}, 'FontSize', ftsz, 'Location', 'eastoutside');
grid on;
% Set x-ticks for each month as approximate week numbers
month_weeks = [1, 5, 9, 14, 18, 23, 27, 31, 36, 40, 45, 49];  % Approximate start of each month in weeks
month_labels = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
set(gca, 'XTick', month_weeks);
set(gca, 'XTickLabel', month_labels);  % Label ticks with month abbreviations

% Ensure the x-axis shows the full range of weeks
xlim([min(weeks), max(weeks)]);

hold off;

% Save the combined figure
print(gcf, '-dpng', '-r1200', 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\temperature_compare\temperature_compare_mixing.png');
print(gcf, '-dpdf', '-r1200', 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\physical_drivers\temperature_compare\temperature_compare_mixing.pdf');
