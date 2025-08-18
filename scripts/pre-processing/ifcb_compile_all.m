% Initialize the final dataframe
ifcb_ucsc_all = [];

% Define the main directory containing subdirectories
main_directory = 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\Santa_Cruz\UCSC_Phyto-20240214T054110Z-001\UCSC_Phyto';

% Get a list of all subdirectories
subdirectories = dir(main_directory);
subdirectories = subdirectories([subdirectories.isdir]); % Remove non-directory entries
subdirectories = subdirectories(~ismember({subdirectories.name}, {'.', '..'})); % Remove '.' and '..' entries

%  Loop through each subdirectory
for subdir_idx = 1:numel(subdirectories)
    subdir_name = subdirectories(subdir_idx).name;
    
    % Get a list of all CSV files in the subdirectory
    files = dir(fullfile(main_directory, subdir_name, '*.csv'));
    
    % Loop through each CSV file
    for file_idx = 1:numel(files)
        file_name = files(file_idx).name;
        file_path = fullfile(main_directory, subdir_name, file_name);
        
        % Print the file name
        disp(['Processing file: ', file_name]);
        
        % Read the CSV file
        data = readtable(file_path);
        
        % Calculate sample volume
        data.sample_volume = ((data.runTime - data.inhibitTime) * 0.25) / 60;
        
        % Divide the first 52 columns by sample volume
        % data{:, 1:52} = data{:, 1:52};
        % 
        % Exclude datetime and filename variables from the computation
        cols_to_include = ~ismember(1:size(data, 2), [find(strcmp(data.Properties.VariableNames, 'dateTime')), find(strcmp(data.Properties.VariableNames, 'fileName'))]);
        
        % Calculate the max of each column, excluding datetime and filename variables
        cells_L_daily = max(data{:, cols_to_include}, 1);
        
        % Create a new row with max values
        new_row = array2table(cells_L_daily, 'VariableNames', data.Properties.VariableNames(cols_to_include));
        
        % Add dateTime column with year, month, and day from the first entry of the dateTime column
        first_dateTime = datetime(data.dateTime(1));
        new_row.dateTime = repmat(first_dateTime, size(new_row, 1), 1);
        
        % Add the new row to the final dataframe
        ifcb_ucsc_all = [ifcb_ucsc_all; new_row];
    end
end

%%  Save the table as a .mat file
save("C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\Santa_Cruz\ifcb_processed\ifcb_ucsc_all.mat", 'ifcb_ucsc_all');


%% Plot
figure(1)
subplot(2,1,1)
plot(ifcb_ucsc_all.dateTime,ifcb_ucsc_all.Cryptophyte,'.')
subplot(2,1,2)
plot(ifcb_ucsc_all.dateTime,ifcb_ucsc_all.Centric,'r.')



%% Make another table that's just daily average
% Initialize the final dataframe
ifcb_ucsc_all_mean = [];

% Define the main directory containing subdirectories
main_directory = 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\Santa_Cruz\UCSC_Phyto-20240214T054110Z-001\UCSC_Phyto';

% Get a list of all subdirectories
subdirectories = dir(main_directory);
subdirectories = subdirectories([subdirectories.isdir]); % Remove non-directory entries
subdirectories = subdirectories(~ismember({subdirectories.name}, {'.', '..'})); % Remove '.' and '..' entries

%  Loop through each subdirectory
for subdir_idx = 1:numel(subdirectories)
    subdir_name = subdirectories(subdir_idx).name;
    
    % Get a list of all CSV files in the subdirectory
    files = dir(fullfile(main_directory, subdir_name, '*.csv'));
    
    % Loop through each CSV file
    for file_idx = 1:numel(files)
        file_name = files(file_idx).name;
        file_path = fullfile(main_directory, subdir_name, file_name);
        
        % Print the file name
        disp(['Processing file: ', file_name]);
        
        % Read the CSV file
        data = readtable(file_path);
        
        % Calculate sample volume
        data.sample_volume = ((data.runTime - data.inhibitTime) * 0.25) / 60;
        
        % Divide the first 52 columns by sample volume
        % data{:, 1:52} = data{:, 1:52};
        % 
        % Exclude datetime and filename variables from the computation
        cols_to_include = 1:52;
        
        % Calculate the max of each column, excluding datetime and filename variables
        cells_L_daily = sum(data{:, 1:52}, 1)./sum(data.sample_volume).*1000;
        
        % Create a new row with max values
        new_row = array2table(cells_L_daily, 'VariableNames', data.Properties.VariableNames(cols_to_include));
        
        % Add dateTime column with year, month, and day from the first entry of the dateTime column
        first_dateTime = datetime(data.dateTime(1));
        new_row.dateTime = repmat(first_dateTime, size(new_row, 1), 1);
        
        % Add the new row to the final dataframe
        ifcb_ucsc_all_mean = [ifcb_ucsc_all_mean; new_row];
    end
end

%%  Save the table as a .mat file
save("C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\Santa_Cruz\ifcb_processed\ifcb_ucsc_all_mean.mat", 'ifcb_ucsc_all_mean');
writetable(ifcb_ucsc_all_mean,"C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\Santa_Cruz\ifcb_processed\ifcb_ucsc_all_mean.csv")


%% Compile biovolume

% Define the base directory
baseDir = 'Q:\Dante\data\MB_Wildfire_Obs\ifcb\scw-features-v4-20240624T172452Z-001\scw-features-v4\';

% Initialize an empty table to store summary statistics
summaryTable = table();

% Get list of all CSV files in the directory and its subdirectories
filePattern = fullfile(baseDir, '**\*.csv');
csvFiles = dir(filePattern);

% Conversion factor from pixels^3 to um^3
conversionFactor = (1 / 2.7)^3;

%% Loop through each file and process it
for k = 1:10
    % Get the file name and path
    baseFileName = csvFiles(k).name;
    fullFileName = fullfile(csvFiles(k).folder, baseFileName);
    
    % Read the CSV file
    data = readtable(fullFileName);
    
    % Convert Biovolume from pixels^3 to um^3
    data.Biovolume = data.Biovolume * conversionFactor;
    
    % Calculate summary statistics
    meanBiovolume = mean(data.Biovolume);
    medianBiovolume = median(data.Biovolume);
    maxBiovolume = max(data.Biovolume);
    minBiovolume = min(data.Biovolume);
    stdBiovolume = std(data.Biovolume);
    numRows = height(data);
    
    % Extract the date and time from the file name
    dateTimeStr = regexp(baseFileName, 'D(\d{8})T(\d{6})', 'tokens');
    dateStr = dateTimeStr{1}{1};
    timeStr = dateTimeStr{1}{2};
    
    dateNum = datenum(dateStr, 'yyyymmdd');
    dateFormatted = datestr(dateNum, 'mm/dd/yyyy');
    
    % Append the results to the summary table
    newRow = {dateFormatted, timeStr, meanBiovolume, medianBiovolume, maxBiovolume, minBiovolume, stdBiovolume, numRows};
    summaryTable = [summaryTable; newRow];
end

% Add column names to the summary table
summaryTable.Properties.VariableNames = {'Date', 'Time', 'MeanBiovolume', 'MedianBiovolume', 'MaxBiovolume', 'MinBiovolume', 'StdBiovolume', 'NumRows'};

% Display the summary table
disp(summaryTable);

%% Optionally, write the summary table to a CSV file
writetable(summaryTable, 'ifcb_biovolume_summary_statistics.csv');
save('ifcb_biovolume_summary_statistics.mat',"summaryTable")




