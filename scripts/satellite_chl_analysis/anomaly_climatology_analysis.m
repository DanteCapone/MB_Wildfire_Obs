% MATLAB script to read and process anomaly and climatology data from CSV
% files for MB and upwelling shadow
addpath('Q:\Dante\data\satellite_chl_CCE')

%% Explore Data by timescale and averaging region
areas = {'EOF2', 'MBarea'};
timescales = {'5d', '15d', 'Mo'};
filePrefix = 'chl';  % Assuming all files have the prefix 'chl'

% Prompt user to select area and timescale
areaSelection = listdlg('PromptString', 'Select area:', ...
                        'SelectionMode', 'single', ...
                        'ListString', areas);
timescaleSelection = listdlg('PromptString', 'Select timescale:', ...
                             'SelectionMode', 'single', ...
                             'ListString', timescales);

% Generate the appropriate filename based on user selection
area = areas{areaSelection};
timescale = timescales{timescaleSelection};
filename = sprintf('%s%sAnomaly_%s.csv', filePrefix, timescale, area);

% Read the selected CSV file
data = readtable(filename);

% Generate datetime columns for start and end dates
data.StartDate = datetime(data.SYear, 1, 1) + days(data.SDay - 1);
data.EndDate = datetime(data.EYear, 1, 1) + days(data.EDay - 1);

% Plot mean time series with standard deviation as shaded area

% Create a figure
clf
figure(1);
hold on;

% Plot the mean as a line plot
plot(startDates, meanValues, '-b*', 'LineWidth', 2);

% Add the black reference line at y = 1
yline(1, '-k', 'LineWidth', 1.5);

% Create the shaded area for the confidence intervals (pastel green)
fill([startDates; flipud(startDates)], ...
     [meanValues - stdDevValues; flipud(meanValues + stdDevValues)], ...
     [0.6, 0.9, 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.3);  % Pastel green

% Shade areas above and below y = 1
above_one = meanValues > 1;
below_one = meanValues < 1;

% Set y-limits dynamically based on mean and standard deviation
ylim([min(meanValues - stdDevValues)*1.2, max(meanValues + stdDevValues)*1.2]);

% Add labels and title
xlabel('Start Date');
ylabel('Mean Anomaly');
title('Time Series of Mean Anomaly with Confidence Intervals and Anomalies');
grid on;

hold off;

dual_monitor=0;
if dual_monitor == 1
    set(gcf, 'Position', [0 1100 1400 1000])
else
    set(gcf, 'Position', [0 100 1400 1000])
end



%% Visualizing the different timescales and averaging regions for the climatologies

close all
areas = {'EOF2shadow', 'MBarea'};
timescales = {'5d', '15d', 'Mo'};
filePrefix = 'chl';  % Assuming all files have the prefix 'chl'


% Iterate over the two areas (EOF2 and MBarea)
for areaIdx = 1:length(areas)
    area = areas{areaIdx};

    % Create a new figure for the current area
    figure('Name', sprintf('Time Series for %s', area), 'NumberTitle', 'off');
    
    % Iterate over the three timescales (5day, 15day, mo)
    for timescaleIdx = 1:length(timescales)
        timescale = timescales{timescaleIdx};
        filename = sprintf('%s%s_%s.csv', filePrefix, timescale, area);

        % Read the CSV file for the current area and timescale
        data = readtable(filename);

        % Generate datetime columns for start and end dates
        data.StartDate = datetime(data.SYear, 1, 1) + days(data.SDay - 1);
        data.EndDate = datetime(data.EYear, 1, 1) + days(data.EDay - 1);

        % Calculate the time series mean and standard deviation
        meanValues = data.Mean;
        stdDevValues = data.StDev;
        startDates = data.StartDate;

        % Create a subplot for the current timescale
        subplot(3, 1, timescaleIdx);
        hold on;

        % Plot the mean time series as a blue line
        plot(startDates, meanValues, '-k', 'LineWidth', 2);

        % Add the confidence intervals in pastel green
        fill([startDates; flipud(startDates)], ...
             [meanValues - stdDevValues; flipud(meanValues + stdDevValues)], ...
             [0.6, 0.9, 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

        % Set y-limits dynamically based on mean and standard deviation
        ylim([min(meanValues - stdDevValues) * 1.2, max(meanValues + stdDevValues) * 1.2]);
        start_date = datetime(2020, 1, 1);
        end_date = datetime(2020, 12, 31);
        xlim([start_date end_date])


        % Add labels, grid, and title
        xlabel('Start Date');
        ylabel('Mean Value');
        title(sprintf('Normal Data Time Series for %s - %s', area, timescale));
        grid on;
        
        hold off;
    end
    
    % Adjust the layout for better visualization
    sgtitle(sprintf('Normal Data Time Series for %s', area));

    dual_monitor=0;
    if dual_monitor == 1
        set(gcf, 'Position', [0 1100 1400 1000])
    else
        set(gcf, 'Position', [0 100 1400 1000])
    end

end



%% 

% Define the areas and timescales
areas = {'EOF2shadow', 'MBarea'};
timescales = {'5d', '15d', 'Mo'};
filePrefix = 'chl';  % Assuming all files have the prefix 'chl'


% Iterate over the two areas (EOF2 and MBarea)
for areaIdx = 1:length(areas)
    area = areas{areaIdx};
    
    % Iterate over the three timescales (5day, 15day, mo)
    for timescaleIdx = 1:length(timescales)
        timescale = timescales{timescaleIdx};
        filename = sprintf('%s%s_%s.csv', filePrefix, timescale, area);

        % Read the CSV file for the current area and timescale
        data = readtable(filename);

        % Generate datetime columns for start and end dates
        data.StartDate = datetime(data.SYear, 1, 1) + days(data.SDay - 1);
        
        % Determine the number of periods in a year for the given timescale
        if strcmp(timescale, '5d')
            numPeriods = 365 / 5;  % 73 periods per year
        elseif strcmp(timescale, '15d')
            numPeriods = round(365 / 15); % 24 periods per year
        else
            numPeriods = 12;       % 12 months per year
        end

        % Initialize variables to store climatology data
        climatology_mean = zeros(numPeriods, 1);
        climatology_max = zeros(numPeriods, 1);
        climatology_min = zeros(numPeriods, 1);
        climatology_median = zeros(numPeriods, 1);
        climatology_std = zeros(numPeriods, 1);
        climatology_ci95 = zeros(numPeriods, 1);  % 95% confidence interval

        % Iterate over each period in a typical year
        for periodIdx = 1:numPeriods
            % Extract data for this specific period across all years
            if strcmp(timescale, '5d')
                periodData = data.Mean(data.SDay >= (periodIdx-1)*5 + 1 & data.SDay <= periodIdx*5);
            elseif strcmp(timescale, '15d')
                periodData = data.Mean(data.SDay >= (periodIdx-1)*15 + 1 & data.SDay <= periodIdx*15);
            else
                periodData = data.Mean(month(data.StartDate) == periodIdx);
            end

            % Compute statistics for the current period
            climatology_mean(periodIdx) = mean(periodData, 'omitnan');
            climatology_max(periodIdx) = max(periodData);
            climatology_min(periodIdx) = min(periodData);
            climatology_median(periodIdx) = median(periodData, 'omitnan');
            climatology_std(periodIdx) = std(periodData, 'omitnan');
            
            % 95% confidence interval (assuming normal distribution)
            n = length(periodData);
            climatology_ci95(periodIdx) = 1.96 * climatology_std(periodIdx) / sqrt(n);
        end

        % Create a table for the climatology data
        climatologyTable = table((1:numPeriods)', climatology_mean, climatology_max, climatology_min, ...
                                 climatology_median, climatology_std, climatology_ci95, ...
                                 'VariableNames', {'Period', 'Mean', 'Max', 'Min', 'Median', 'Std', 'CI95'});

        % Save the climatology as CSV
        csvFilename = sprintf('%s_climatology_%s_%s.csv', filePrefix, timescale, area);
        writetable(climatologyTable, csvFilename);
        
        % Save the climatology as MAT file
        matFilename = sprintf('%s_climatology_%s_%s.mat', filePrefix, timescale, area);
        save(matFilename, 'climatologyTable');
    end
end


%% Plot climatologies for MB and shadow region side by side

% Initialize storage for climatology data
climatologies = struct();

% Define the areas and timescales
areas = {'EOF2shadow', 'MBarea'};
timescales_for_struct = {'d5', 'd15', 'Mo'};

% Load climatology data from the saved files for both EOF2 and MBarea
for areaIdx = 1:length(areas)
    area = areas{areaIdx};
    
    % Load climatology data for all timescales (5day, 15day, mo)
    for timescaleIdx = 1:length(timescales)
        timescale_name = timescales_for_struct{timescaleIdx};
        timescale = timescales{timescaleIdx};
        csvFilename = sprintf('%s_climatology_%s_%s.csv', filePrefix, timescale, area);
        climatologies.(area).(timescale_name) = readtable(csvFilename);
    end
end

MB_satellite_climatologies_struct=climatologies;
save('Q:\Dante\data\MB_Wildfire_Obs\processed_data\satellite\new_climatologies_anomalies\MB_satellite_climatologies_struct.mat','MB_satellite_climatologies_struct')
% Create a new figure for plotting climatologies
figure('Name', 'Climatology Comparison: EOF2 vs MBarea', 'NumberTitle', 'off');

% Define month labels
monthLabels = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};

% Iterate over the three timescales (5-day, 15-day, monthly)
for timescaleIdx = 1:length(timescales)
    timescale = timescales{timescaleIdx};
    timescale_name = timescales_for_struct{timescaleIdx};

    % Create a subplot for the current timescale
    subplot(3, 1, timescaleIdx);
    hold on;

    % Plot climatology for EOF2
    plot(climatologies.EOF2shadow.(timescale_name).Period, ...
         climatologies.EOF2shadow.(timescale_name).Mean, '-b', 'LineWidth', 1.5, 'DisplayName', 'EOF2');
     
    % Plot climatology for MBarea
    plot(climatologies.MBarea.(timescale_name).Period, ...
         climatologies.MBarea.(timescale_name).Mean, '-r', 'LineWidth', 1.5, 'DisplayName', 'MBarea');

    % Add shaded error bars for standard deviation (confidence intervals are optional)
    fill([climatologies.EOF2shadow.(timescale_name).Period; flipud(climatologies.EOF2shadow.(timescale_name).Period)], ...
         [climatologies.EOF2shadow.(timescale_name).Mean - climatologies.EOF2shadow.(timescale_name).Std; ...
          flipud(climatologies.EOF2shadow.(timescale_name).Mean + climatologies.EOF2shadow.(timescale_name).Std)], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    fill([climatologies.MBarea.(timescale_name).Period; flipud(climatologies.MBarea.(timescale_name).Period)], ...
         [climatologies.MBarea.(timescale_name).Mean - climatologies.MBarea.(timescale_name).Std; ...
          flipud(climatologies.MBarea.(timescale_name).Mean + climatologies.MBarea.(timescale_name).Std)], ...
         'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Add limits
    % Determine the number of periods in a year for the given timescale
   if strcmp(timescale_name, 'd5')
        xlim([1 73]);  % 5-day intervals over a year
        xticks([1, 7, 13, 19, 25, 31, 37, 43, 49, 55, 61, 67]);  % Approximate positions for the first of each month
        xticklabels(monthLabels);  % Label them with month names

    elseif strcmp(timescale_name, 'd15')
        xlim([1 24]);  % 15-day intervals over a year
        xticks([1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23]);  % Approximate positions for the first of each month
        xticklabels(monthLabels);  % Label them with month names

    else
        xlim([1 12]);  % Monthly intervals
        xticks(1:12);  % Custom ticks for each month
        xticklabels(monthLabels);  % Month labels
    end

    % Customize plot appearance
    xlabel('Period');
    ylabel('Mean Value');
    title(sprintf('Climatology - %s', timescale_name));
    % legend('show');
    grid on;
    
    hold off;
end

% Add a super title for the whole figure
sgtitle('Climatology Comparison for EOF2 and MBarea');

dual_monitor=1;
if dual_monitor == 1
    set(gcf, 'Position', [0 1400 1400 800])
else
    set(gcf, 'Position', [0 100 1400 1000])
end


%% Finally make plot of 2020
close all
% Define constants for the plot
linewdt = 2;
month_labs = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
xvals = 1:73;  % 73 periods for 5-day climatology

lag_days = 6;
start_date = datetime(2020, 8, 16);
end_date = datetime(2020, 9, 22);
lag_date = datetime(2020, 8, 20 + lag_days);
august_date = datetime(2020, 9, 5);

%By day
start_date_day = day(start_date, 'dayofyear');
end_date_day = day(end_date, 'dayofyear');
lag_date_day = day(lag_date, 'dayofyear');
august_date_day = day(august_date, 'dayofyear');

%By week
start_date_week = week(start_date, 'weekofyear');
end_date_week = week(end_date, 'weekofyear');
lag_date_week = week(lag_date, 'weekofyear');
august_date_week = week(august_date, 'weekofyear');

% By 5-day window
% Each 5-day window starts at day 1, 6, 11, 16, ..., so we divide by 5 and round up
start_date_5day = ceil(start_date_day / 5);
end_date_5day = ceil(end_date_day / 5);
lag_date_5day = ceil(lag_date_day / 5);
august_date_5day = ceil(august_date_day / 5);

%Years
start_year=2016;
end_year=2023;


%For shading
shade_start_day = start_date_day;
shade_end_day = end_date_day;

shade_start_5day = start_date_5day;
shade_end_5day = end_date_5day;

shade_start_week = start_date_week;
shade_end_week = end_date_week;




% Load the 5-day climatology data (EOF2shadow)
% Assuming climatologies.EOF2shadow.d5 contains Mean and Std columns
% You should have climatology data for Mean and Std deviation

% Plot climatology with error bounds
figure;

% Plot the shaded error region (standard deviation)
fill([xvals, fliplr(xvals)], ...
     [climatologies.EOF2shadow.d5.Mean' + climatologies.EOF2shadow.d5.Std', ...
      fliplr(climatologies.EOF2shadow.d5.Mean' - climatologies.EOF2shadow.d5.Std')], ...
      'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
hold on;

% Plot the climatology line
plot(xvals, climatologies.EOF2shadow.d5.Mean, '-', 'LineWidth', linewdt, 'Color', '#979da6');
hold on
plot(xvals, climatologies.MBarea.d5.Mean, '-', 'LineWidth', linewdt, 'Color', '#979da6');

% Customize the x-axis for months
xlabel('Month');
ylabel('Average Chlorophyll-a [mg/m^3]');
xlim([0 73]);
ylim([0 25]);
set(gca, 'XTick', 1:73/12:73, 'XTickLabel', month_labs, 'XTickLabelRotation', 0);
grid on;

% Read the normal 2020 data from 'chl5d_EOF2shadow.csv'
% Make sure the file has columns like SYear, SDay, Mean for each time point
filename = 'chl5d_EOF2shadow.csv';  % Path to the CSV file
data_2020 = readtable(filename);

% Filter data for the year 2020
data_2020 = data_2020(data_2020.SYear == 2020, :);

% Map the day of the year (SDay) to 5-day periods (1-73)
dayOfYear_2020 = data_2020.SDay;
periodIndex_2020 = ceil(dayOfYear_2020 / 5);  % Convert day of the year to 5-day periods

% Plot 2020 normal data on the right y-axis
yyaxis right;

% Plot the 2020 data
plot(periodIndex_2020, data_2020.Mean, '-', 'Color', '#66b398', 'DisplayName', '2020', 'LineWidth', linewdt);

% Customize the right y-axis
ylabel('Chlorophyll-a (2020) [mg/m^3]');
legend(["Climatology", "2020"], 'AutoUpdate', 'Off');


% Add lines for 8 day lag and August Fire Smoke
xline(lag_date_5day, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
xline(august_date_5day, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':');




% Add shaded area to the plot
xline(start_date_5day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
fill([shade_start_5day, shade_start_5day, shade_end_5day, shade_end_5day], ...
    [min(ylim) * 4.5, max(ylim) * 1.5, max(ylim) * 1.5 min(ylim) * 4.5], ...
    [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4); 
xline(end_date_5day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');

% Final plot customization
ylim([0 25])

set(gca, 'XTick', 1:73/12:73, 'XTickLabel', month_labs, 'XTickLabelRotation', 0);
grid on;

% Adjust font and other aesthetics
ax = gca;
ax.FontSize = 16;
ax.FontName = 'Times';

title('Satellite Chlorophyll-a Data for EOF2 (2020 vs Climatology)');

dual_monitor=0;
if dual_monitor == 1
    set(gcf, 'Position', [0 1400 1400 600])
else
    set(gcf, 'Position', [0 100 1400 1000])
end
