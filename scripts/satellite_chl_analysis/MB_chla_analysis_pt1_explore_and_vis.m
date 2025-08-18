% Analysis: Basic Exploration of Satellite Chlorophyll Data
% Date: 4.24.2024
% Data Link: https://spg-satdata.ucsd.edu/MBay/MBay.htm


%% Part 1: Testing
addpath('Q:\Dante\')

%% Using Monthly Data first
%% Define the file path for Chlorophyll data
file_path = 'Q:\Dante\data\MB_Wildfire_Obs\satellite_chl\2020\C2020_chl_5day';

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

filename = 'Q:\Dante\Wildfire_Obs\animations\chlDataAnimation5day_2020.gif';

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
file_path = 'Q:\Dante\data\MB_Wildfire_Obs\satellite_chl\2020\C2020_chl_day';

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

filename = 'Q:\Dante\Wildfire_Obs\animations\chlDataAnimationdaily_2020.gif';

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
% Original color palette
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

% Number of additional colors needed
additionalColorsNeeded = 22 - size(colors, 1);

% Generate additional colors by interpolating existing colors
newColors = zeros(additionalColorsNeeded, 3);
for i = 1:additionalColorsNeeded
    % Pick two random indices from the existing color palette
    idx1 = randi(size(colors, 1));
    idx2 = randi(size(colors, 1));
    while idx2 == idx1  % Ensure different indices for interpolation
        idx2 = randi(size(colors, 1));
    end
    
    % Interpolation factor
    t = rand();  % Random interpolation factor between 0 and 1
    
    % Linearly interpolate between two colors
    newColors(i, :) = (1 - t) * colors(idx1, :) + t * colors(idx2, :);
end

% Combine original and new colors
extendedColors = [colors; newColors];

% Optionally, you can shuffle the colors to randomize their order
% extendedColors = extendedColors(randperm(size(extendedColors, 1)), :);

% Display the extended color palette
disp(extendedColors);

%%

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
