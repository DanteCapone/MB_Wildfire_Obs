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