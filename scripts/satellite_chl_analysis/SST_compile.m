% Define the directory where your CSV files are located
directory = 'Q:\Dante\data\MB_Wildfire_Obs\satellite\sst\sst';
fileList = dir(fullfile(directory, '*.csv'));

% Check if fileList is empty
if isempty(fileList)
    error('No CSV files found in the specified directory. Please check the path.');
end

% Initialize the main structure
dataStruct = struct();

% First pass: Load Mo data to establish monthly means and std deviations
monthlyData = struct();

for i = 1:length(fileList)
    % Get the file name
    fileName = fileList(i).name;
    
    % Extract region and time period from the file name
    if contains(fileName, 'EOF2shadow')
        region = 'EOF2shadow';
    elseif contains(fileName, 'MBarea')
        region = 'MBarea';
    else
        continue; % Skip files that don't match expected pattern
    end
    
    % Check if the file corresponds to monthly data
    if contains(fileName, 'Mo')
        timePeriod = 'Mo';
    else
        continue; % Skip files that are not monthly
    end
    
    % Read the CSV file
    filePath = fullfile(directory, fileName);
    data = readtable(filePath);
    
    % Calculate month from start day-of-year
    data.Month = month(datetime(data.SYear, 1, 1) + days(data.SDay - 1));
    
    % Initialize monthly statistics for the current region if not already set
    if ~isfield(monthlyData, region)
        monthlyData.(region) = struct();
    end
    
    % Calculate monthly mean and std deviation for each year
    for monthNum = 1:12
        monthData = data.Mean(data.Month == monthNum);
        monthlyData.(region).mean(monthNum) = mean(monthData, 'omitnan');
        monthlyData.(region).std(monthNum) = std(monthData, 'omitnan');
    end
end

% Second pass: Load all data and filter d5 and d15 data based on monthly stats
for i = 1:length(fileList)
    % Get the file name
    fileName = fileList(i).name;
    
    % Extract region and time period from the file name
    if contains(fileName, 'EOF2shadow')
        region = 'EOF2shadow';
    elseif contains(fileName, 'MBarea')
        region = 'MBarea';
    else
        continue; % Skip files that don't match expected pattern
    end
    
    % Assign a valid field name for time periods
    if contains(fileName, 't5d')
        timePeriod = 'd5';
    elseif contains(fileName, 't15d')
        timePeriod = 'd15';
    elseif contains(fileName, 'Mo')
        timePeriod = 'Mo';
    else
        continue; % Skip files that don't match expected pattern
    end

    % Read the CSV file
    filePath = fullfile(directory, fileName);
    data = readtable(filePath);
    
    % Calculate the month based on the start date for alignment with Mo data
    data.Month = month(datetime(data.SYear, 1, 1) + days(data.SDay - 1));
    
    % Filter d5 and d15 data based on monthly stats
    if ~strcmp(timePeriod, 'Mo')
        for j = 1:height(data)
            monthNum = data.Month(j);
            monthlyMean = monthlyData.(region).mean(monthNum);
            monthlyStd = monthlyData.(region).std(monthNum);
            
            % Mark values as NaN if they deviate by more than 3 std from the monthly mean
            if abs(data.Mean(j) - monthlyMean) > 5
                data.Mean(j) = NaN;
            end
        end
    end
    
    % Convert SYear and SDay to day of year
    startDate = datetime(data.SYear, 1, 1) + days(data.SDay - 1);
    endDate = datetime(data.EYear, 1, 1) + days(data.EDay - 1);
    
    % Initialize region structure if it doesn't exist
    if ~isfield(dataStruct, region)
        dataStruct.(region) = struct();
    end
    
    % Initialize time period structure if it doesn't exist
    if ~isfield(dataStruct.(region), timePeriod)
        dataStruct.(region).(timePeriod) = struct('data', [], 'startDate', [], 'endDate', []);
    end
    
    % Store the data in the structure
    dataStruct.(region).(timePeriod).data = data;
    dataStruct.(region).(timePeriod).startDate = startDate;
    dataStruct.(region).(timePeriod).endDate = endDate;
end

% Display the structure to verify contents
disp(dataStruct);


%% Save 5day time series for subsequent analysis
% Directory for saving the processed data
saveDirectory = 'Q:\Dante\data\MB_Wildfire_Obs\processed_data\satellite\sst';

% Ensure save directory exists
if ~isfolder(saveDirectory)
    mkdir(saveDirectory);
end

% Process only the EOF2shadow region for d5 data
region = 'EOF2shadow';

% Check if 'd5' data exists for EOF2shadow region
if isfield(dataStruct.(region), 'd5')
    % Extract d5 data
    d5Data = dataStruct.(region).('d5').data;
    timeData = dataStruct.(region).('d5').startDate;
    sstData = d5Data.Mean;

    % Sort data by time to ensure retime works correctly
    [timeData, sortIdx] = sort(timeData);
    sstData = sstData(sortIdx);

    % Create timetable for interpolation
    d5Timetable = timetable(timeData, sstData, 'VariableNames', {'sst'});

    % Interpolate to daily values
    dailyTimetable = retime(d5Timetable, 'daily', 'linear');

    % Prepare the output structure
    EOF2shadow_d5_daily_interpolated=timetable();
    EOF2shadow_d5_daily_interpolated.datetime = dailyTimetable.timeData;
    EOF2shadow_d5_daily_interpolated.datetime.TimeZone='UTC';
    EOF2shadow_d5_daily_interpolated.sst = dailyTimetable.sst;

    % Save the interpolated data as a .mat file
    saveFileName = fullfile(saveDirectory, [region '_d5_daily_interpolated.mat']);
    save(saveFileName, 'EOF2shadow_d5_daily_interpolated');
    
    disp(['Saved daily interpolated data for region ' region ' to ' saveFileName]);
else
    disp(['No d5 data found for region ' region]);
end



%% Climatology: Define the regions and time periods
% Define the regions and time periods
regions = {'EOF2shadow', 'MBarea'};
timePeriods = {'d5', 'd15', 'Mo'};
incrementsPerYear = [73, 24, 12];  % Number of increments for 5-day, 15-day, and monthly periods
incrementDays = [5, 15, 30];  % Approximate number of days per increment for each time period

% Loop over each region
for r = 1:length(regions)
    region = regions{r};
    
    % Loop over each time period
    for i = 1:length(timePeriods)
        period = timePeriods{i};
        numIncrements = incrementsPerYear(i);  % Number of increments
        daysPerIncrement = incrementDays(i);   % Approximate days per increment
        
        % Initialize arrays for climatology data and day-of-year mapping
        climatologyMean = nan(numIncrements, 1);
        climatologyMedian = nan(numIncrements, 1);
        climatologyStd = nan(numIncrements, 1);
        climatologyDOY = nan(numIncrements, 1);  % Store day-of-year for each increment
        
        % Debugging output to verify processing
        disp(['Processing climatology for region: ' region ', period: ' period]);
        
        % Loop over each increment within the year
        for increment = 1:numIncrements
            allDataForIncrement = [];
            
            % Loop over all rows in the data for this period and region
            for j = 1:height(dataStruct.(region).(period).data)
                % Check if the current row corresponds to the current increment
                if mod(j - 1, numIncrements) + 1 == increment
                    allDataForIncrement = [allDataForIncrement; dataStruct.(region).(period).data.Mean(j)];
                end
            end
            
            % Verify data collection for each increment
            disp(['Increment ' num2str(increment) ': Data count = ' num2str(length(allDataForIncrement))]);
            
            % Calculate climatology for the current increment
            climatologyMean(increment) = mean(allDataForIncrement, 'omitnan');
            climatologyMedian(increment) = median(allDataForIncrement, 'omitnan');
            climatologyStd(increment) = std(allDataForIncrement, 'omitnan');
            
            % Calculate the approximate DOY for this increment
            climatologyDOY(increment) = (increment - 1) * daysPerIncrement + (daysPerIncrement / 2);
        end
        
        % Store the climatology arrays and DOY mapping in the structure
        dataStruct.(region).(period).climatology.mean = climatologyMean;
        dataStruct.(region).(period).climatology.median = climatologyMedian;
        dataStruct.(region).(period).climatology.std = climatologyStd;
        dataStruct.(region).(period).climatology.DOY = climatologyDOY;  % Store the DOY mapping
    end
end



%% Plot all data across time for each region
% Define the regions and time periods
regions = {'EOF2shadow', 'MBarea'};
timePeriods = {'d5', 'd15', 'Mo'};

for r = 1:length(regions)
    region = regions{r};
    figure('Name', [region ' - All Time Periods (Filtered)'], 'NumberTitle', 'off');
    
    for t = 1:length(timePeriods)
        period = timePeriods{t};
        
        % Check if the region has data for this period
        if isfield(dataStruct.(region), period)
            % Extract the data
            timeData = dataStruct.(region).(period).startDate;
            meanData = dataStruct.(region).(period).data.Mean;
            
            % Filter out outliers based on 3 standard deviations from a moving mean window
            movMean = movmean(meanData, 3);  % Moving mean with window size of 3
            movStd = movstd(meanData, 3);    % Moving std deviation with window size of 3
            outlierMask = abs(meanData - movMean) > (3 * movStd);
            meanDataFiltered = meanData;  % Copy data
            meanDataFiltered(outlierMask) = NaN;  % Replace outliers with NaN
            
            % Plot in a subplot
            subplot(3, 1, t);
            hold on;
            
            % Plot the filtered mean data
            plot(timeData, meanDataFiltered, '-', 'DisplayName', 'Filtered Mean');
            
            % Set plot limits
            ylim([0 20]);
            
            % Add labels, title, and legend
            xlabel('Time');
            ylabel('Mean Value');
            title([region ' - ' period ' (Filtered)']);
            legend('show', 'AutoUpdate', 'off');
            
            hold off;
        end
    end
end

%% Plot rate of change (derivative) of the filtered data for each region and time period
close all
for r = 1:length(regions)
    region = regions{r};
    figure('Name', [region ' - Rate of Change'], 'NumberTitle', 'off');
    
    for t = 1:length(timePeriods)
        period = timePeriods{t};
        
        % Check if the region has data for this period
        if isfield(dataStruct.(region), period)
            % Extract the data
            timeData = dataStruct.(region).(period).startDate;
            meanData = dataStruct.(region).(period).data.Mean;
            
            % Filter out outliers as before
            movMean = movmean(meanData, 3);
            movStd = movstd(meanData, 3);
            outlierMask = abs(meanData - movMean) > (3 * movStd);
            meanDataFiltered = meanData;
            meanDataFiltered(outlierMask) = NaN;
            
            % Calculate the rate of change (derivative)
            deltaTime = diff(timeData);  % Time increments (assumes equal intervals)
            deltaMean = diff(meanDataFiltered);  % Change in mean values
            rateOfChange = deltaMean ./ days(deltaTime);  % Calculate rate of change per day
            
            % Plot in a subplot
            subplot(3, 1, t);
            hold on;
            
            % Plot the rate of change
            plot(timeData(1:end-1), rateOfChange, '-', 'DisplayName', 'Rate of Change');
            
            % Set plot limits
            ylim([-2 2]);  % Adjust based on expected rate of change range
            
            % Add labels, title, and legend
            xlabel('Time');
            ylabel('Rate of Change');
            title([region ' - ' period ' (Rate of Change)']);
            legend('show', 'AutoUpdate', 'off');
            
            hold off;
        end
    end
end


%% Plot 2020 data vs. climatology for each region
for r = 1:length(regions)
    region = regions{r};
    figure('Name', [region ' - 2020 vs Climatology'], 'NumberTitle', 'off');
    
    for t = 1:length(timePeriods)
        period = timePeriods{t};
        
        % Check if the region has data for this period
        if isfield(dataStruct.(region), period)
            % Extract 2020 data and day-of-year
            timeIndices2020 = find(dataStruct.(region).(period).data.SYear == 2020);
            timeData2020 = dataStruct.(region).(period).startDate(timeIndices2020);
            meanData2020 = dataStruct.(region).(period).data.Mean(timeIndices2020);
            dayOfYear2020 = day(timeData2020, 'dayofyear');
            
            % Extract climatology data and corresponding DOY
            climatologyDOY = dataStruct.(region).(period).climatology.DOY;
            climatologyMean = dataStruct.(region).(period).climatology.mean;
            climatologyStd = dataStruct.(region).(period).climatology.std;
            
            % Plot in a subplot
            subplot(2, 2, t);
            hold on;
            
            % Plot the 2020 mean data
            plot(dayOfYear2020, meanData2020, '-', 'DisplayName', '2020 Mean');
            
            % Plot the climatology mean and standard deviation as shaded area
            fill([climatologyDOY; flipud(climatologyDOY)], ...
                 [climatologyMean - climatologyStd; flipud(climatologyMean + climatologyStd)], ...
                 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Climatology Std');
            plot(climatologyDOY, climatologyMean, '-r', 'DisplayName', 'Climatology Mean');
            

            if t==3
                ylim([10 20])

            end


            % Add labels, title, and legend
            xlabel('Day of Year');
            ylabel('Mean Value');
            title([region ' - ' period ' (2020 vs Climatology)']);
            legend('show', 'AutoUpdate', 'off');
            
            hold off;
        end
    end
end
