% Analysis: Basic Exploration of Satellite Chlorophyll Data
% Date: 4.24.2024
% Data Link: https://spg-satdata.ucsd.edu/MBay/MBay.htm


%% Using Monthly Data first
%% Define the file path for Chlorophyll data
file_path = 'E:\satellite_chl\2020\C2020_chl_month';

%% List all HDF files in the directory
files = dir(fullfile(file_path, '*.hdf'));
% Extract the names
fileNames = {files.name};

% Sort files by the numeric part of the name which represents the date
[~, idx] = sort(str2double(regexp(fileNames, '\d+', 'match', 'once')));
fileNames = fileNames(idx)';

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

%% Create an animation of monthly Chlorophyll data
hFig = figure;
axis tight manual;

filename = 'E:\satellite_chl\animations\chlDataAnimationmonthly_2020.gif';

for f = 1:length(fileNames)
    fullFilePath = fullfile(file_path, fileNames{f});
    info = hdfinfo(fullFilePath);
    dataSetName = info.SDS(1).Name;
    data = hdfread(fullFilePath, dataSetName); % Specify dataset name

    % Parse filename to get year and month
    [year, month] = parseFilename(fileNames{f});


    % Create coastline mask
    coastlineMask = data == -1;
    [coastlineY, coastlineX] = find(coastlineMask);

    % Correct for the signed byte issue
    data = int16(data);  % Convert data to a 16-bit integer to handle negative values correctly
    negativeIndices = data < 0;
    data(negativeIndices) = data(negativeIndices) + 256;

    % Filter invalid data
    data(data < 2 | data > 254) = NaN; % Excludes invalid pixel values
    data(data == 0) = 0.01; % Represents no data
    data(data == 1) = NaN; % Coastline, excluded from analysis
    data(data == 255) = 66.834; % Invalid data (land/clouds), excluded from analysis



  % Extract the boundaries of the coastline
    boundaries = bwboundaries(coastlineMask);

    % Convert pixel values to Chlorophyll concentration
    % validData = data > 0 & data < 255;
    chlConcentration = NaN(size(data));
    chlConcentration = 10.^(0.015 * double(data) - 2.0);

    for k = 1:size(chlConcentration, 3)
        imagesc(lon(1,:), (lat(:, 1)), squeeze(chlConcentration(:, :, k)));
        shading interp
        hold on; % Hold the current plot
        % Plot each boundary
        for b = 1:length(boundaries)
            boundary = boundaries{b};
            validIndices = boundary(:,1) <= size(lat,1) & boundary(:,2) <= size(lon,2);
            % Ensure that the indices are not out of the matrix dimension limits
            if all(validIndices)
                % Extract the corresponding latitude and longitude values
                boundaryLat = lat(sub2ind(size(lat), boundary(validIndices,1), boundary(validIndices,2)));
                boundaryLon = lon(sub2ind(size(lon), boundary(validIndices,1), boundary(validIndices,2)));
                plot(boundaryLon, boundaryLat, 'k', 'LineWidth', 1);
            end
        end
        hold off;
        colormap([1 1 1; jet]); % Add white color for NaN values    
        caxis([0.01 15]); % Adjust color axis based on expected range
        c = colorbar;
        xlabel(c, 'Chlorophyll-a [mg m^-^3]');        
        title(sprintf('Chlorophyll-a [mg m^-^3]: %s, %d', month, year));
        set(gca, 'YDir', 'normal');
        set(gcf,'Position',[300 50 600 500])

        [imind, cm] = rgb2ind(frame2im(getframe(hFig)), 256);

        if f == 1 && k == 1
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf);
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append');
        end

        pause(0.01);
    end
end


%% 15-day data

%% Define the file path for Chlorophyll data
file_path = 'E:\satellite_chl\2020\C2020_chl_15day';

%% List all HDF files in the directory
files = dir(fullfile(file_path, '*.hdf'));
% Extract the names
fileNames = {files.name};

% Sort files by the numeric part of the name which represents the date
[~, idx] = sort(str2double(regexp(fileNames, '\d+', 'match', 'once')));
fileNames = fileNames(idx)';

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

%% Create an animation of 15-day Chlorophyll data
hFig = figure;
axis tight manual;

filename = 'E:\satellite_chl\animations\chlDataAnimation15day_2020.gif';

for f = 1:length(fileNames)
    fullFilePath = fullfile(file_path, fileNames{f});
    info = hdfinfo(fullFilePath);
    dataSetName = info.SDS(1).Name;
    data = hdfread(fullFilePath, dataSetName); % Specify dataset name

    % Parse filename to get year and month
    datestr = parseFilenameDay(fileNames{f});


    % Create coastline mask
    coastlineMask = data == -1;
    [coastlineY, coastlineX] = find(coastlineMask);

    % Correct for the signed byte issue
    data = int16(data);  % Convert data to a 16-bit integer to handle negative values correctly
    negativeIndices = data < 0;
    data(negativeIndices) = data(negativeIndices) + 256;

    % Filter invalid data
    data(data < 2 | data > 254) = NaN; % Excludes invalid pixel values
    data(data == 0) = 0.01; % Represents no data
    data(data == 1) = NaN; % Coastline, excluded from analysis
    data(data == 255) = 66.834; % Invalid data (land/clouds), excluded from analysis



  % Extract the boundaries of the coastline
    boundaries = bwboundaries(coastlineMask);

    % Convert pixel values to Chlorophyll concentration
    % validData = data > 0 & data < 255;
    chlConcentration = NaN(size(data));
    chlConcentration = 10.^(0.015 * double(data) - 2.0);

    for k = 1:size(chlConcentration, 3)
        imagesc(lon(1,:), (lat(:, 1)), squeeze(chlConcentration(:, :, k)));
        shading interp
        hold on; % Hold the current plot
        % Plot each boundary
        for b = 1:length(boundaries)
            boundary = boundaries{b};
            validIndices = boundary(:,1) <= size(lat,1) & boundary(:,2) <= size(lon,2);
            % Ensure that the indices are not out of the matrix dimension limits
            if all(validIndices)
                % Extract the corresponding latitude and longitude values
                boundaryLat = lat(sub2ind(size(lat), boundary(validIndices,1), boundary(validIndices,2)));
                boundaryLon = lon(sub2ind(size(lon), boundary(validIndices,1), boundary(validIndices,2)));
                plot(boundaryLon, boundaryLat, 'k', 'LineWidth', 1);
            end
        end
        hold off;
        colormap([1 1 1; jet]); % Add white color for NaN values    
        caxis([0.01 15]); % Adjust color axis based on expected range
        c = colorbar;
        xlabel(c, 'Chlorophyll-a [mg m^-^3]');        

        xlabel('Longitude')
        ylabel('Lattitude')
        title(sprintf('Chlorophyll-a [mg m^-^3]: %s', datestr));
        set(gca, 'YDir', 'normal');
        set(gcf,'Position',[300 50 1200 1000])

        [imind, cm] = rgb2ind(frame2im(getframe(hFig)), 256);

        if f == 1 && k == 1
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf);
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append');
        end

        pause(0.001);
    end
end

%% Daily data

%% Define the file path for Chlorophyll data
file_path = 'E:\satellite_chl\2020\C2020_chl_day';

%% List all HDF files in the directory
files = dir(fullfile(file_path, '*.hdf'));
% Extract the names
fileNames = {files.name};

% Sort files by the numeric part of the name which represents the date
[~, idx] = sort(str2double(regexp(fileNames, '\d+', 'match', 'once')));
fileNames = fileNames(idx)';

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

%% Create an animation of monthly Chlorophyll data
hFig = figure;
axis tight manual;

filename = 'E:\satellite_chl\animations\chlDataAnimationdaily_2020.gif';

for f = 1:length(fileNames)
    fullFilePath = fullfile(file_path, fileNames{f});
    info = hdfinfo(fullFilePath);
    dataSetName = info.SDS(1).Name;
    data = hdfread(fullFilePath, dataSetName); % Specify dataset name

    % Parse filename to get year and month
    datestr = parseFilenameDay(fileNames{f});


    % Create coastline mask
    coastlineMask = data == -1;
    [coastlineY, coastlineX] = find(coastlineMask);

    % Correct for the signed byte issue
    data = int16(data);  % Convert data to a 16-bit integer to handle negative values correctly
    negativeIndices = data < 0;
    data(negativeIndices) = data(negativeIndices) + 256;

    % Filter invalid data
    data(data < 2 | data > 254) = NaN; % Excludes invalid pixel values
    data(data == 0) = 0.01; % Represents no data
    data(data == 1) = NaN; % Coastline, excluded from analysis
    data(data == 255) = 66.834; % Invalid data (land/clouds), excluded from analysis



  % Extract the boundaries of the coastline
    boundaries = bwboundaries(coastlineMask);

    % Convert pixel values to Chlorophyll concentration
    % validData = data > 0 & data < 255;
    chlConcentration = NaN(size(data));
    chlConcentration = 10.^(0.015 * double(data) - 2.0);

    for k = 1:size(chlConcentration, 3)
        imagesc(lon(1,:), (lat(:, 1)), squeeze(chlConcentration(:, :, k)));
        shading interp
        hold on; % Hold the current plot
        % Plot each boundary
        for b = 1:length(boundaries)
            boundary = boundaries{b};
            validIndices = boundary(:,1) <= size(lat,1) & boundary(:,2) <= size(lon,2);
            % Ensure that the indices are not out of the matrix dimension limits
            if all(validIndices)
                % Extract the corresponding latitude and longitude values
                boundaryLat = lat(sub2ind(size(lat), boundary(validIndices,1), boundary(validIndices,2)));
                boundaryLon = lon(sub2ind(size(lon), boundary(validIndices,1), boundary(validIndices,2)));
                plot(boundaryLon, boundaryLat, 'k', 'LineWidth', 1);
            end
        end
        hold off;
        colormap([1 1 1; jet]); % Add white color for NaN values    
        caxis([0.01 15]); % Adjust color axis based on expected range
        c = colorbar;
        xlabel(c, 'Chlorophyll-a [mg m^-^3]');        

        xlabel('Longitude')
        ylabel('Lattitude')
        title(sprintf('Chlorophyll-a [mg m^-^3]: %s', datestr));
        set(gca, 'YDir', 'normal');
        set(gcf,'Position',[300 50 600 500])

        [imind, cm] = rgb2ind(frame2im(getframe(hFig)), 256);

        if f == 1 && k == 1
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf);
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append');
        end

        pause(0.001);
    end
end



%% Averages across all space for all years
    colors = [
    0.00, 0.45, 0.74;  % MATLAB Blue
    0.85, 0.33, 0.10;  % MATLAB Red
    0.93, 0.69, 0.13;  % MATLAB Yellow
    0.49, 0.18, 0.56;  % MATLAB Purple
    0.47, 0.67, 0.19;  % MATLAB Green
    0.30, 0.75, 0.93;  % MATLAB Light Blue
    0.64, 0.08, 0.18;  % MATLAB Dark Red
    1.00, 0.41, 0.16;  % Orange
    0.08, 0.17, 0.55;  % Navy Blue
    0.29, 0.00, 0.51;  % Indigo
    0.00, 0.50, 0.00;  % Office Green
    0.74, 0.72, 0.42;  % Olive Green
];

base_path = 'E:\\satellite_chl\\';
dataset_name = 'Chlorophyll_Data';
start_year = 2012;
end_year = 2023;
createNewFig = true;

for year = start_year:end_year
    yearIdx=year-start_year+1;
    file_path = fullfile(base_path, sprintf('%d', year), sprintf('C%d_chl_month', year));
    plotChlorophyllTimeSeries(file_path, createNewFig,colors,yearIdx);
    createNewFig = false;  % Only create a new figure on the first call
end

set(gcf,'Position',[0 1000 1200 800])
hold off;  


%% Plot daily data 

base_path = 'E:\\satellite_chl\\';
dataset_name = 'Chlorophyll_Data';
start_year = 2020;
end_year = 2020;
createNewFig = true;

for year = start_year:end_year
    yearIdx=year-start_year+1;
    file_path = fullfile(base_path, sprintf('%d', year), sprintf('C%d_chl_day', year));
    plotChlorophyllTimeSeries(file_path, createNewFig,colors,yearIdx);
    createNewFig = false;  % Only create a new figure on the first call
end

set(gcf,'Position',[0 1000 1200 800])
hold off;  



%% Climatologies

base_path = 'E:\\satellite_chl\\';
results = computeChlorophyllMetrics(base_path, 2012, 20);

%%
start_year=2012;
end_year=2023;
base_path = 'E:\\satellite_chl\\';
data_type = '5day';  % Change as needed: 'daily', '5day', '15day', 'month'
results = computeChlorophyllMetrics(base_path, data_type, 2012, 2023);


%% Plot climatologies
% Now plot the output in a figure with two subplots

% Calculate the months for x-axis labels
% Each 5-day period's month is calculated assuming the period starts on Jan 1
% months = arrayfun(@(x) datestr(datetime(start_year,1,1) + days(x*5-2.5), 'mmm'), 1:73, 'UniformOutput', false);

figure(1);
% First subplot for climatology with shaded confidence intervals
subplot(2,1,1);
xvals = 1:length(results.climatology);
fill([xvals, fliplr(xvals)], [results.climatology' + results.error', fliplr(results.climatology' - results.error')], [0.9 0.9 1], 'EdgeColor', 'none'); hold on;
plot(xvals, results.climatology, 'g-','LineWidth',2);
title('Monterey Bay Chlorophyll-a Climatology ');
xlabel('Month');
ylabel('Chlorophyll Concentration');
xlim([0 73])
set(gca, 'XTick', 1:73/12:73, 'XTickLabel', months, 'XTickLabelRotation', 45);  % Set custom x-axis labels
grid on;

% Second subplot for anomalies with a zero line
subplot(2,1,2); hold on;
% legendInfo = cell(end_year - start_year + 1, 1);  % Initialize cell array for legend labels
% for i = 1:size(results.anomalies, 2)  % Loop through each year
i=9;
        yearLabel = sprintf('%d', start_year + i - 1);  % Create a label for the current year
        plot(1:size(results.anomalies, 1), results.anomalies(:, i), '-', 'Color', colors(i, :), 'DisplayName', yearLabel,'LineWidth',2);
        legendInfo{i} = yearLabel;  % Store the year label
% end
plot(xlim, [0 0], 'k-','LineWidth',2);  % Adding a zero line
title('Chlorophyll-a Anomaly');
xlabel('Month');
ylabel('Chlorophyll-a [mg/m^3]');
xlim([0 73])
set(gca, 'XTick', 1:73/12:73, 'XTickLabel', months, 'XTickLabelRotation', 45);  % Use the same x-axis labels
grid on;
% legend(legendInfo, 'Location', 'northeastoutside');  % Display the legend outside the plot

set(gcf,'Position',[0 50 1200 600])


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

set(gcf,'Position',[0 1000 1200 800])


%% ========= Functions

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
        case 'daily'
            datePattern = '(?i)C(\d{4})(\d{3})_Chl_comp_mapped\.hdf';
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
                if strcmp(data_type, 'daily')
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
