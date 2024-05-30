% Analysis: Basic Exploration of Satellite Chlorophyll Data
% Date: 4.24.2024
% Data Link: https://spg-satdata.ucsd.edu/MBay/MBay.htm

addSubfoldersToPath('Q:\Dante\')
savepath

%% Part 2: Analysis
load('merra_MB_bc_avg_daily_plt.mat')
load('satellite_chl_MB_5day_anomalies.mat')
%% Climatologies

base_path = 'Q:\Dante\data\MB_Wildfire_Obs\satellite_chl';

%% Compute climatologies
start_year=2008;
end_year=2023;
years=start_year:1:end_year;
base_path = 'Q:\Dante\data\MB_Wildfire_Obs\satellite_chl';
data_type = 'day';  % Change as needed: 'daily', '5day', '15day', 'month'

%Implement function to calculate climatology (Takes a while)
results = computeChlorophyllMetrics(base_path, data_type, 2008, 2023);
results.years=years;

%% 

%% Plot climatologies
% Now plot the output in a figure with two subplots
load('merra_MB_bc_avg_daily_plt.mat')


% For axis labels
month_labs = {'Jan', 'Feb', 'Ma', 'Apr', 'May', 'Jun', ...
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};


% Calculate the months for x-axis labels
% Each 5-day period's month is calculated assuming the period starts on Jan 1
% month_5weeks = arrayfun(@(x) datestr(datetime(start_year,1,1) + days(x*5-2.5), 'mmm'), 1:73, 'UniformOutput', false);
clf
figure(1);
% First subplot for climatology with shaded confidence intervals
subplot(2,1,1);
xvals = 1:length(results.climatology);
fill([xvals, fliplr(xvals)], [results.climatology' + results.error', fliplr(results.climatology' - results.error')], [0.9 0.9 1], 'EdgeColor', 'none'); hold on;
plot(xvals, results.climatology, 'g-','LineWidth',2);
title('Monterey Bay Chlorophyll-a Climatology (5-day)');
xlabel('Month');
ylabel('Chlorophyll-a [mg/m^3] Climatology');
% xlim([0 73])
% set(gca, 'XTick', 1:73/12:73, 'XTickLabel', month_labs, 'XTickLabelRotation', 45);  % Use the same x-axis labels
grid on;

% Second subplot for anomalies with a zero line
yyaxis right
years_sel=[2020];
legendInfo = cell(end_year - start_year + 1, 1);  % Initialize cell array for legend labels
index = 1;  % Initialize index for storing labels in legendInfo
counter=0;

for i = 1:size(results.anomalies, 2)  % Loop through each year
    currentYear = start_year + i - 1;  % Calculate the current year
    if ismember(currentYear, years_sel)
        counter=counter+1;
        yearLabel = sprintf('%d', currentYear);  % Create a label for the current year
        % plot(1:size(results.anomalies, 1), results.anomalies(:, i), '-', 'Color', colors(counter, :), 'DisplayName', yearLabel, 'LineWidth', 2);
        shade_anomaly(1:size(results.anomalies, 1), results.anomalies(:, i))
        hold on
        ylabel('Chlorophyll-a Anomaly 2020 [mg/m^3]');
        legendInfo{index} = yearLabel;  % Store the year label in the appropriate position
        index = index + 1;  % Increment the index for the next valid year
    end
end
ylim([-2 10])
xlim([0 366])


% Update legend to only include selected years
plot(xlim, [0 0], 'k-','LineWidth',2,'Color', [0,0,0, 0.5],'LineStyle','-');  % Adding a zero line
% legend(legendInfo(1:index-1), 'Location', 'northwest');  % Adjust legend to show only included years
title('Monterey Bay Chlorophyll-a Anomaly');
% xlabel('Month');
ylabel('Chlorophyll-a [mg/m^3]');
% xlim([0 73])
% set(gca, 'XTick', 1:73/12:73, 'XTickLabel', month_labs, 'XTickLabelRotation', 45);  % Use the same x-axis labels
set(gca, 'XTick', 1:30:366, 'XTickLabel', month_labs, 'XTickLabelRotation', 45);  % Use the same x-axis labels

grid on;
% legend(legendInfo, 'Location', 'northeastoutside');  % Display the legend outside the plot



subplot(2,1,2); hold on;

%Change BC to g/m3 for plotting
merra_bc_plt.Median_BCCMASS=merra_bc_plt.Median_BCCMASS.*1000;
merra_bc_plt.anomaly=merra_bc_plt.anomaly.*1000;

yyaxis left;
plotClimatology(merra_bc_plt((merra_bc_plt.datetime.Year > 2002) & ...
                            (merra_bc_plt.datetime.Year ~= 2008) & ...
                            (merra_bc_plt.datetime.Year ~= 2016) & ...
                            (merra_bc_plt.datetime.Year ~= 2020) & ...
                            (merra_bc_plt.datetime.Year ~= 2021),:),'Median_BCCMASS')
ylabel(['BC Surface',newline(),'Mass Concentration [g/m^3] Climatology'])
ylim([0 1.5e-3])
alpha(0.1)
hold on 
yyaxis right;
% plot(merra_bc_plt.dayOfYear(merra_bc_plt.datetime.Year == 2020),merra_bc_plt.anomaly(merra_bc_plt.datetime.Year == 2020),'Color',colorz(3,:),'LineStyle','-','LineWidth',2);
shade_anomaly(merra_bc_plt.dayOfYear(merra_bc_plt.datetime.Year == 2020),merra_bc_plt.anomaly(merra_bc_plt.datetime.Year == 2020))
xlim([0 366])
ylabel(['BC Surface',newline(),'Mass Concentration 2020 Anomaly [g/m^3]'])
set(gca, 'XTick', 1:30:366, 'XTickLabel', month_labs, 'XTickLabelRotation', 45);  % Use the same x-axis labels


% set(gcf,'Position',[-400 1200 1600 600])
set(gcf,'Position',[0 0 1600 600])

saving=1;
if saving==1
    if data_type=='day'
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\satellite_chl\satellite_chl_climatology_5day_MB_fire_years_with_merra_daily.png"]);
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\satellite_chl\satellite_chl_climatology_5day_MB_fire_years_with_merra_daily.pdf"]);
    end
end
saving=0;



%% Plot

close all
figure(1)
for ii=1:width(results.anomalies)
    plot(results.anomalies(:,ii))
    xlim([0 366])
    hold on
end
datetick('x')
legend()
plot([0 366], [0 0],'LineStyle','--','Color','k')

set(gcf,'Position',[0 0 1200 800])


%% Write time series of Chl 
% Extracting the years and anomalies from the structure
yearss = results.years;
anomalies = results.anomalies;

% Initialize a table for satellite chlorophyll data with a datetime column
satellite_chl_MB = table([], [], [], [], 'VariableNames', {'Year', '5_Day_Window', 'Chlorophyll', 'datetime'});
satellite_chl_MB = table([], [], [], [], 'VariableNames', {'Year', 'Day', 'Chlorophyll', 'datetime'});

% Populate the table
for j = 1:length(yearss)
    current_year = yearss(j);

    if data_type=="5day"
        num_windows = size(anomalies, 1); % number of 5-day windows (73)
        for i = 1:num_windows
            % Calculate the start day of the 5-day window
            start_day = (i - 1) * 5 + 1; % first day of the current window
            date_str = datetime(current_year, 1, 1) + days(start_day - 1); % calculate date
            
            % Assume each row in anomalies corresponds to a 5-day window
            new_row = {current_year, i, anomalies(i, j), date_str};
            satellite_chl_MB = [satellite_chl_MB; new_row];
        end
    end
end

% save('satellite_chl_MB_5day_anomalies.mat','satellite_chl_MB')
% writetable(satellite_chl_MB,'satellite_chl_MB_5day_anomalies.csv')

% satellite_chl_MB_daily_anomalies=results;
% save('satellite_chl_MB_daily_anomalies.mat','satellite_chl_MB_daily_anomalies')

%% ============ Comparative analysis Satellite vs. Merra-2

%% Now correlate with Merra


satellite_all_raw=outerjoin(merra_bc_plt,satellite_chl_MB,'Keys','datetime','Type','left');

% Remove rows where Chlorophyll is NaN
validRows = ~isnan(satellite_all_raw.Chlorophyll);  % Logical index for rows with non-NaN Data2
satellite_all = satellite_all_raw(validRows, :);


%% Xcorr
start_year=2020;
end_year=2020;
satellite_all_sel=satellite_all(satellite_all.datetime_merra_bc_plt.Year>=start_year & satellite_all.datetime_merra_bc_plt.Year <= end_year,:);
x=detrend(satellite_all_sel.anomaly-mean(satellite_all_sel.anomaly));
y=detrend(satellite_all_sel.Chlorophyll-mean(satellite_all_sel.Chlorophyll));

% Initialize variables to store coefficients and lags

start_date_plot=satellite_all_sel.datetime_merra_bc_plt(1);
end_date_plot=satellite_all_sel.datetime_merra_bc_plt(end);
start_date_plot=datetime(2020,1,1);
end_date_plot=datetime(2020,12,31);


%Cross-correlate
[cc,lags] = xcorr(x,y,'coeff');


 % Compute first derivative of cc values
cc_derivative = diff(cc(:));
lags_midpoints = (lags(1:end-1) + lags(2:end)) / 2; % Compute midpoints of lags for derivative

% Find indices where derivative changes from positive to negative (relative maxima)
maxima_indices = find(cc_derivative(1:end-1) > 0 & cc_derivative(2:end) < 0) + 1;
maxima_lags = lags(maxima_indices);
    
maxima_lags


%% Plot the cross-correlation
plotting=0;
clf
figure(1)
% stem(lags, cc);
subplot(2,1,1)
plot(lags,cc,'LineStyle','-','LineWidth',2)
hold on

i=find(cc==max(cc));
plot(lags(i),cc(i),'*','MarkerSize',12)
title("Cross-correlation between AOD 500nm &"+newline+ "Satellite Chlorophyll-a");
xlabel('Lags (Days)');
ylabel('Normalized Cross-correlation');
% xlim([-365 365])

subplot(2,1,2)
yyaxis left
plot(merra_bc_plt.datetime,merra_bc_plt.anomaly,'r')
hold on 
ylabel(['BC Surface',newline(),'Mass Concentration Anomaly [g/m^3] '])

yyaxis right
plot(satellite_chl_MB.datetime,satellite_chl_MB.Chlorophyll,'g-')
ylabel('5-day Chlorophyll-a Anomaly [mg/m^3]');
xlim([start_date_plot end_date_plot])

set(gcf,'Position',[0 0 1200 800])


saving=0;
if saving==1
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\satellite_chl\Xcorr_satellite_chl_climatology_5day_MB_with_merra.png"]);
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\satellite_chl\Xcorr_satellite_chl_climatology_5day_MB_with_merra.pdf"]);
end


%% Compute Climatolog by grid cell
data_type = 'month';  % Change as needed: 'daily', '5day', '15day', 'month'

results_grid = computeChlorophyllMetricsByGrid(base_path, data_type, 2020, 2020);

%% Coordinate calculations
maxX = 307; % Maximum x index for chlorophyll data
maxY = 318; % Maximum y index for chlorophyll data

[x, y] = meshgrid(0:maxX, 0:maxY);
lon = -122.68085 + 0.003290 * x;
lat = 37.180664 - 0.002701 * y;

%% Define bounds for subsetting (example bounds, adjust as needed)
latBounds = [36.3219, 37.180664]; % latitude range
lonBounds = [-122.68085, -121.6675]; % longitude range

subsetMask = (lat >= latBounds(1) & lat <= latBounds(2)) & ...
             (lon >= lonBounds(1) & lon <= lonBounds(2));


% Create a figure and plot mean anomalies
clf
figure(1);
pcolor(lon(1,:), (lat(:, 1)),results_grid.climatology);  % Use imagesc to scale the data to color map
shading interp
col = colorbar;  % Adds a color bar to indicate the scale
ylabel(col, 'Chlorophyll Concentration (mg/m^3)');  % Label the colorbar
axis xy;  % Ensures the axes are correctly oriented
xlabel('Longitude');
ylabel('Latitude');
title('Monterey Bay 2020 Chlorophll-a Climatology');

set(gcf,'Position',[0 900 800 600])

saving=0;
if saving==1
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\Santa_Cruz\figures\satellite\satellite_chla_2020_climatology.png"]);
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\Santa_Cruz\figures\satellite\satellite_chla_2020_climatology.pdf"]);
end
saving=0;

%% ========= Functions ========= 

%% Function to parse the filename and determine the year and month
function [year, month] = parseFilename(filename)
    yearStr = filename(2:5);
    year = str2double(yearStr);
    doyStart = str2double(filename(6:8));
    doyEnd = str2double(filename(13:15));
    avgDoy = round(mean([doyStart doyEnd])); % Compute average day of the year

    % Define month day ranges
    monthDays = [0 31 59 90 120 151 181 212 243 273 304 334 366];
    months = {'January', 'February', 'March', 'April', 'May', 'June', ...
              'July', 'August', 'September', 'October', 'November', 'December'};

    % Find the corresponding month
    monthIdx = find(avgDoy > monthDays, 1, 'last');
    month = months{monthIdx};
end


function [dateStr] = parseFilenameDay(filename)
    % Extract the year and day of year from the filename
    yearStr = filename(2:5);
    year = str2double(yearStr);
    doyStr = filename(6:8);
    doy = str2double(doyStr);

    % Determine if it is a leap year
    if mod(year, 4) == 0
        if mod(year, 100) == 0
            if mod(year, 400) == 0
                leapYear = true;  % Leap year
            else
                leapYear = false; % Not a leap year
            end
        else
            leapYear = true;  % Leap year
        end
    else
        leapYear = false;  % Not a leap year
    end

    % Compute the date from the day of year
    startDate = datenum(year, 1, 1); % January 1 of the given year
    dateNum = startDate + doy - 1; % Convert DOY to MATLAB date number
    
    % Format the date string as 'Month, Day Year'
    dateStr = datestr(dateNum, 'mmmm, dd yyyy');
end


%% Averaging Functions

function monthlyChlAvg = plotChlorophyllTimeSeries(file_path, createNewFig,colors,yearIdx)
    files = dir(fullfile(file_path, '*.hdf'));
    fileNames = {files.name};
    monthlyChlAvg = nan(1, length(fileNames));  % Initialize with NaN to handle missing data

    for f = 1:length(fileNames)
        fullFilePath = fullfile(file_path, fileNames{f});
        info = hdfinfo(fullFilePath);
        dataSetName = info.SDS(1).Name;
        data = hdfread(fullFilePath, dataSetName); % Specify dataset name
        data = int16(data); 
        data(data < 0) = data(data < 0) + 256;
        data(data < 2 | data > 254) = NaN;
        validData = data > 1 & data < 255;
        chlConcentration = NaN(size(data));
        chlConcentration(validData) = 10.^(0.015 * double(data(validData)) - 2.0);
        monthlyChlAvg(f) = nanmean(chlConcentration(:));  % Calculate mean only for valid data
    end

    if createNewFig
        figure;
    end
    months = 1:length(monthlyChlAvg);
    
    % Extract year from file path using regexp
    yearMatch = regexp(file_path, '\d{4}', 'match');
    if ~isempty(yearMatch)
        year = str2double(yearMatch{1});  % Assumes the first match is the year
    else
        year = NaN;  % Default or error value if no year is found
    end
    
    plot(months, monthlyChlAvg, '-', 'LineWidth', 2, 'MarkerSize', 8, 'Color', colors(yearIdx, :), 'DisplayName', sprintf('Year %d', str2double(regexp(file_path, '\d{4}', 'match', 'once'))));
    hold on;
    xlabel('Month');
    ylabel('Average Chlorophyll-a [mg m^-^3]');
    title('Monthly Average Chlorophyll-a Concentration Over Multiple Years');
    grid on;
    legend('Location','bestoutside');
end


function dailyChlAvg = plotChlorophyllTimeSeriesDaily(file_path, createNewFig, colors, yearIdx)
    files = dir(fullfile(file_path, '*.hdf'));
    fileNames = {files.name};
    dailyChlAvg = nan(1, length(fileNames));  % Initialize with NaN to handle missing data

    for f = 1:length(fileNames)
        fullFilePath = fullfile(file_path, fileNames{f});
        info = hdfinfo(fullFilePath);
        dataSetName = info.SDS(1).Name;
        data = hdfread(fullFilePath, dataSetName); % Read specific dataset
        data = int16(data); 
        data(data < 0) = data(data < 0) + 256;
        data(data < 2 | data > 254) = NaN;
        validData = data > 1 & data < 255;
        chlConcentration = NaN(size(data));
        chlConcentration(validData) = 10.^(0.015 * double(data(validData)) - 2.0);
        dailyChlAvg(f) = nanmean(chlConcentration(:));  % Calculate mean only for valid data
    end

    if createNewFig
        figure;
    end
    daysOfYear = 1:length(dailyChlAvg);
    
    % Extract year from file path using regexp
    yearMatch = regexp(file_path, '\d{4}', 'match');
    year = str2double(yearMatch{1});  % Assume the first match is the year
    
    plot(daysOfYear, dailyChlAvg, '-', 'LineWidth', 2, 'MarkerSize', 8, 'Color', colors(yearIdx, :), 'DisplayName', sprintf('Year %d', year));
    hold on;
    xlabel('Day of Year');
    ylabel('Average Chlorophyll-a [mg m^-^3]');
    title('Daily Average Chlorophyll-a Concentration');
    grid on;
    legend show;
end

%% Climatology Functions

function output = computeChlorophyllMetrics(base_path, data_type, start_year, end_year)
    % Define the aggregation period based on data type
    switch data_type
        case 'day'
            datePattern = '(?i)C(\d{4})(\d{3})_chl_comp_mapped\.hdf';
            numPeriods = 366;  % Consider leap years
        case '5day'
            datePattern = 'C(\d{4})(\d{3})(\d{4})(\d{3})__comp\.hdf';
            numPeriods = 73;   % 5-day periods in a year
            periodLength = 5;
        case '15day'
            datePattern = 'C(\d{4})(\d{3})(\d{4})(\d{3})__comp\.hdf';
            numPeriods = 24;   % 15-day periods in a year
            periodLength = 15;
        case 'month'
            datePattern = 'C(\d{4})(\d{3})(\d{4})(\d{3})__comp\.hdf';
            numPeriods = 12;   % Months in a year
        otherwise
            error('Invalid data type specified.');
    end

    % Initialize storage for all data across all years
    dailyData = nan(numPeriods, end_year - start_year + 1); 

    for year = start_year:end_year
        yearIndex = year - start_year + 1;
        folder_path = fullfile(base_path, sprintf('%d', year), sprintf('C%d_chl_%s', year, data_type));
        files = dir(fullfile(folder_path, '*.hdf'));

        for f = 1:length(files)
            filename = files(f).name;
            tokens = regexp(filename, datePattern, 'tokens');
            if ~isempty(tokens)
                startDay = str2double(tokens{1}(2));
                if strcmp(data_type, 'day')
                    period = str2double(tokens{1}(2));
                elseif strcmp(data_type, 'month')
                    endDay = str2double(tokens{1}(4));
                    representativeDay = round(mean([startDay, endDay]));
                    period = dayOfYearToMonth(representativeDay, year);
                else
                    if exist('periodLength', 'var')
                        period = ceil(startDay / periodLength);
                    else
                        error('Period length is not set for data type %s', data_type);
                    end
                end
            else
                error('Filename does not match expected format');
            end


            fullFilePath = fullfile(folder_path, filename);
            info = hdfinfo(fullFilePath);
            dataSetName = info.SDS(1).Name;
            data = hdfread(fullFilePath, dataSetName);
            data = int16(data);
            data(data < 0) = data(data < 0) + 256;
            data(data < 2 | data > 254) = NaN;
            validData = data > 1 & data < 255;
            chlConcentration = NaN(size(data));
            chlConcentration(validData) = 10.^(0.015 * double(data(validData)) - 2.0);
            dailyData(period, year - start_year + 1) = nanmean(chlConcentration(:));
        end
    end

    % Calculate climatology and confidence intervals
    climatology = nanmedian(dailyData, 2);
    sem = nanstd(dailyData, 0, 2) / sqrt(sum(~isnan(dailyData), 2));
    ci95 = 1.96 * sem;
    anomalies = dailyData - climatology;

    % Prepare the output structure
    output.climatology = climatology;
    output.error = ci95;
    output.anomalies = anomalies;
end


% function output = computeChlorophyllMetrics(base_path, start_year, end_year)
%     dailyData = nan(366, end_year - start_year + 1);  % Assuming leap years are included
% 
%     for year = start_year:end_year
%         file_path = fullfile(base_path, sprintf('%d', year), sprintf('C%d_chl_day', year));
%         files = dir(fullfile(file_path, '*.hdf'));
% 
%         for f = 1:length(files)
%             filename = files(f).name;
%             tokens = regexp(filename, '(?i)C(\d{4})(\d{3})_chl_comp_mapped\.hdf', 'tokens');
%             if ~isempty(tokens)
%                 tokens = tokens{1}; % First and only match set
%                 dayOfYear = str2double(tokens{2});
%             else
%                 error('Filename does not match expected format');
%             end
% 
%             fullFilePath = fullfile(file_path, filename);
%             info = hdfinfo(fullFilePath);
%             dataSetName = info.SDS(1).Name;
%             data = hdfread(fullFilePath, dataSetName);
%             data = int16(data);
%             data(data < 0) = data(data < 0) + 256;
%             data(data < 2 | data > 254) = NaN;
%             validData = data > 1 & data < 255;
%             chlConcentration = NaN(size(data));
%             chlConcentration(validData) = 10.^(0.015 * double(data(validData)) - 2.0);
%             dailyData(dayOfYear, year - start_year + 1) = nanmean(chlConcentration(:));
%         end
%     end
% 
%     % Reshape dailyData into a single vector
%     dailyData = reshape(dailyData, [], 1);  % Reshape to 366*(end_year - start_year + 1) automatically
% 
%     climatology = nanmean(reshape(dailyData, 366, []), 2);
%     sem = nanstd(reshape(dailyData, 366, []), 0, 2) / sqrt(sum(~isnan(reshape(dailyData, 366, [])), 2));
%     ci95 = 1.96 * sem;
%     anomalies = reshape(dailyData, 366, []) - climatology;
% 
%     output.climatology = climatology;
%     output.error = ci95;
%     output.anomalies = anomalies;
% end

%% Helper Functions
function addSubfoldersToPath(rootPath)
    % Add the root path first
    addpath(rootPath);
    
    % Generate a list of all subdirectories
    allSubFolders = genpath(rootPath);
    
    % Parse the string of folders returned by genpath
    % genpath returns a path string that MATLAB uses which separates
    % folders with a semicolon (;) on Windows or a colon (:) on UNIX/Mac.
    % We need to split this string to get individual paths.
    folderList = strsplit(allSubFolders, pathsep);
    
    % Add each subdirectory to the MATLAB path
    for i = 1:length(folderList)
        addpath(folderList{i});
    end
    
    % Optionally, save the path if you want MATLAB to remember these settings
    % in future sessions
    savepath;
end

function output = computeChlorophyllMetricsByGrid(base_path, data_type, start_year, end_year)
    % Define the aggregation period based on data type
    switch data_type
        case 'day'
            datePattern = '(?i)C(\d{4})(\d{3})_chl_comp_mapped\.hdf';
            numPeriods = 366;  % Consider leap years
            periodLength = 1;  % Every day is a period by itself
        case '5day'
            datePattern = 'C(\d{4})(\d{3})(\d{4})(\d{3})__comp\.hdf';
            periodLength = 5;
            numPeriods = 73;   % 5-day periods in a year
        case '15day'
            datePattern = 'C(\d{4})(\d{3})(\d{4})(\d{3})__comp\.hdf';
            periodLength = 15;
            numPeriods = 24;   % 15-day periods in a year
        case 'month'
            datePattern = 'C(\d{4})(\d{3})(\d{4})(\d{3})__comp\.hdf';
            periodLength = 30;  % Approximate, mainly for consistency in handling
            numPeriods = 12;   % Months in a year
        otherwise
            error('Invalid data type specified.');
    end
    
    totalPeriods = numPeriods * (end_year - start_year + 1);
    % Initialize datetime array
    timeLabels = repmat(datetime('2000-01-01'), totalPeriods, 1);  % Initialize with a default date

    % Initialize data storage with dimensions for lat/lon
    gridData = [];

    for year = start_year:end_year
        folder_path = fullfile(base_path, sprintf('%d', year), sprintf('C%d_chl_%s', year, data_type));
        files = dir(fullfile(folder_path, '*.hdf'));
        yearIndex = year - start_year;

        for f = 1:length(files)
            filename = files(f).name;
            tokens = regexp(filename, datePattern, 'tokens');
            if isempty(tokens)
                error('Filename does not match expected format');
            end
            startDay = str2double(tokens{1}(2));
            period = determinePeriod(startDay, data_type, periodLength);
            periodIndex = yearIndex * numPeriods + period;

            % Assign datetime for each period start
            if yearIndex == 0 && f == 1  % Assume first file of first year sets the pattern
                timeLabels(periodIndex) = datetime(year, 1, 1) + days((period - 1) * periodLength);
            end

            fullFilePath = fullfile(folder_path, filename);
            data = readAndProcessChlData(fullFilePath);

            % Initialize or aggregate data
            if isempty(gridData)
                gridData = nan([size(data), totalPeriods]);
            end
            gridData(:, :, periodIndex) = nanmean(cat(4, gridData(:, :, periodIndex), data), 4);
        end
    end

    % Compute statistics
    output = computeStatistics(gridData, timeLabels);

end

%%
function period = determinePeriod(day, dataType, periodLength)
    switch dataType
        case 'day'
            period = day;
        otherwise
            period = ceil(day / periodLength);
    end
end

function data = readAndProcessChlData(filePath)
    info = hdfinfo(filePath);
    dataSetName = info.SDS(1).Name;
    data = hdfread(filePath, dataSetName);
    data = int16(data);
    data(data < 0) = data(data < 0) + 256;  % Adjust data range
    validData = (data > 1 & data < 255);
    data(~validData) = NaN;  % Set invalid data to NaN
    data(validData) = 10.^(0.015 * double(data(validData)) - 2.0);  % Convert to concentration
    return
end

function stats = computeStatistics(data, timeLabels)
    climatology = nanmean(data, 3);  % Median over all periods
    sem = nanstd(data, 0, 3) ./ sqrt(sum(~isnan(data), 3));  % Standard error of mean
    ci95 = 1.96 * sem;  % 95% confidence interval
    anomalies = data - climatology;

    % Compute the mean of anomalies over time
    meanAnomalies = nanmean(anomalies, 3);

    stats.climatology = climatology;
    stats.error = ci95;
    stats.anomalies = anomalies;
    stats.meanAnomalies = meanAnomalies;
    stats.timeLabels = timeLabels;
end