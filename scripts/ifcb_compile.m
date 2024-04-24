% Define the folder path
% base_path = 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\Santa_Cruz\SCW IFCB Classified Sep-Aug 2020-20231011T193037Z-001\SCW IFCB Classified Sep-Aug 2020\';
base_path = 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\Santa_Cruz\UCSC_Phyto-20240214T054110Z-001\UCSC_Phyto';

% Get a list of all subfolders in the base path
subfolders = dir(base_path);
subfolders = subfolders([subfolders.isdir]);  % filter out files, keep only directories
subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'}));  % remove current and parent directory

% Initialize arrays to store categories and results
all_categories = {};
combined_results = table();

% Traverse each subfolder
for folder_idx = 1:length(subfolders)
    folder_path = fullfile(base_path, subfolders(folder_idx).name);
    mat_files = dir(fullfile(folder_path, '*.mat'));

    % Step 1: Find unique categories across all files
    for i = 1:length(mat_files)
        data = load(fullfile(folder_path, mat_files(i).name), 'TBclass_above_threshold');
        cell_array = data.TBclass_above_threshold;
        unique_categories = unique(cell_array);
        all_categories = union(all_categories, unique_categories);
    end

    % Step 2: Construct the table for this folder path
    for i = 1:length(mat_files)
        file_name = mat_files(i).name;

        % Extract datetime
        pattern = 'D(\d{4})(\d{2})(\d{2})T(\d{2})(\d{2})(\d{2})';
        tokens = regexp(file_name, pattern, 'tokens');
        date_str = tokens{1};
        datetime_str = sprintf('%s/%s/%s-%s:%s:%s', date_str{2}, date_str{3}, date_str{1}, date_str{4}, date_str{5}, date_str{6});
        dt = datetime(datetime_str, 'InputFormat', 'MM/dd/yyyy-HH:mm:ss');

        data = load(fullfile(folder_path, file_name), 'TBclass_above_threshold');
        cell_array = data.TBclass_above_threshold;

        category_counts = zeros(1, length(all_categories));
        for j = 1:length(all_categories)
            category_counts(j) = sum(strcmp(cell_array, all_categories{j}));
        end

        tmp_table = array2table(category_counts, 'VariableNames', matlab.lang.makeValidName(all_categories));
        tmp_table.filename = {file_name};
        tmp_table.datetime = dt;

        combined_results = [combined_results; tmp_table];
    end
end

% Group by day and compute daily average
combined_results.Day = dateshift(combined_results.datetime, 'start', 'day');
daily_avg = varfun(@mean, combined_results, 'GroupingVariables', 'Day', 'InputVariables', matlab.lang.makeValidName(all_categories));

% Remove the 'mean_' prefix from the variable names
for i = 1:length(all_categories)
    oldVarName = sprintf('mean_%s', matlab.lang.makeValidName(all_categories{i}));
    daily_avg.Properties.VariableNames{oldVarName} = matlab.lang.makeValidName(all_categories{i});
end

% Display the final daily average table
disp(daily_avg);


%% Save
ifcb_aug_oct_2020=daily_avg;
ifcb_aug_oct_2020.datetime=ifcb_aug_oct_2020.Day;
save("C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\Santa_Cruz\IFCB_compiled_8_10_2020.mat","ifcb_aug_oct_2020")