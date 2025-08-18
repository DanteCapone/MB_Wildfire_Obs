% Load the dataset
load('Q:\Dante\data\MB_Wildfire_Obs\processed_data\ifcb_processed\ifcb_bray.mat');


% Step 1: Extract the 'datetime' column to use for filtering
datetime_values = ifcbbray.datetime; % Assuming 'datetime' is a datetime column

% Step 2: Extract taxa counts from columns 2 to 52
taxa_counts = ifcbbray{:, 2:52}; % Community composition data (taxa counts)

% Step 3: Identify samples that have another sample exactly 6 days apart
valid_indices = false(size(datetime_values)); % Preallocate a logical array

for i = 1:length(datetime_values)

    % Look for another sample exactly 6 days before or after
    match_before = find(datetime_values == datetime_values(i) - days(6), 1);
    match_after = find(datetime_values == datetime_values(i) + days(6), 1);
    
    if ~isempty(match_before) || ~isempty(match_after)
        valid_indices(i) = true; % Mark this sample as valid if it has a match
    end
end

% Step 4: Filter the taxa counts and datetime values based on valid indices
filtered_taxa_counts = taxa_counts(valid_indices, :);
filtered_datetime_values = datetime_values(valid_indices);

% Step 5: Compute Bray-Curtis dissimilarity on the filtered data
function D = braycurtis_batch(X)
    n = size(X, 1); % Number of samples
    D = zeros(n*(n-1)/2, 1); % Preallocate distance matrix in pdist format
    k = 1;
    for i = 1:n-1
        for j = i+1:n
            numerator = sum(abs(X(i,:) - X(j,:)));
            denominator = sum(X(i,:) + X(j,:));
            D(k) = numerator / denominator;
            k = k + 1;
        end
    end
end

bray_dist = braycurtis_batch(filtered_taxa_counts);

% Step 6: Perform NMDS using 'mdscale' on Bray-Curtis dissimilarities
nmds_dim = 2; % For 2D NMDS

% Convert the condensed distance vector into a full distance matrix
square_brays_dist = squareform(bray_dist); % Converts vector to square matrix

% Perform NMDS
[nmds_coords, stress] = mdscale(square_brays_dist, nmds_dim, 'Criterion', 'stress');

%%
% Step 1: Perform k-means clustering with different numbers of clusters on NMDS coordinates
max_clusters = 10; % Maximum number of clusters to test
silhouette_avg = zeros(max_clusters-1, 1); % Store average silhouette scores

for k = 2:max_clusters
    % Apply k-means clustering
    cluster_labels = kmeans(nmds_coords, k, 'Replicates', 10, 'Distance', 'sqeuclidean');
    
    % Calculate silhouette scores for current number of clusters
    silhouette_avg(k-1) = mean(silhouette(nmds_coords, cluster_labels)); % Silhouette score
end


%% Step 2: Plot the silhouette scores to find the optimal number of clusters
figure;
plot(2:max_clusters, silhouette_avg, '-o');
title('Silhouette Scores for Different Numbers of Clusters');
xlabel('Number of Clusters');
ylabel('Average Silhouette Score');

% Step 3: Use the elbow method to find the optimal number of clusters
sum_of_squared_distances = zeros(max_clusters-1, 1);
for k = 2:max_clusters
    % Apply k-means clustering
    [~, ~, sumd] = kmeans(nmds_coords, k, 'Replicates', 10, 'Distance', 'sqeuclidean');
    sum_of_squared_distances(k-1) = sum(sumd); % Total within-cluster sum of squares
end

% Plot SSE (Elbow Method)
figure;
plot(2:max_clusters, sum_of_squared_distances, '-o');
title('Elbow Method for Optimal Number of Clusters');
xlabel('Number of Clusters');
ylabel('Sum of Squared Distances');

% Step 4: Based on the silhouette and elbow method plots, select the optimal number of clusters
optimal_clusters = 4; % Example: select optimal number of clusters (adjust based on results)

% Step 5: Apply k-means with the optimal number of clusters
[cluster_labels, cluster_centers] = kmeans(nmds_coords, optimal_clusters, 'Replicates', 10, 'Distance', 'sqeuclidean');

%% Step 6: Plot the NMDS results and color by cluster
close all
% Define updated date range for labeling points of interest
label_start_date = datetime(2020, 8, 21);
label_end_date = datetime(2020, 8, 2);

% Define 30-day buffer period
buffer_start_date = label_start_date - days(30);
buffer_end_date = label_end_date + days(30);

% Proximity threshold for NMDS space
proximity_threshold = 0.1;

% Step 1: Identify points in the given date range and 30-day buffer
labeled_indices = find(filtered_datetime_values >= label_start_date & filtered_datetime_values <= label_end_date);
buffer_indices = find(filtered_datetime_values >= buffer_start_date & filtered_datetime_values <= buffer_end_date);

% Step 2: Perform k-means clustering
optimal_clusters = 4; % Number of clusters from elbow method
[cluster_labels, cluster_centers] = kmeans(nmds_coords, optimal_clusters, 'Replicates', 10, 'Distance', 'sqeuclidean');

% Step 3: Plot NMDS points - color all points by cluster using consistent colors
figure;
all_indices = 1:length(filtered_datetime_values);
colors_cluster = lines(optimal_clusters); % Use consistent colormap for clusters

% Plot all points colored by cluster using consistent cluster colors
scatter_handles_clusters = gscatter(nmds_coords(:, 1), nmds_coords(:, 2), cluster_labels, colors_cluster, '.', 20);
hold on;

% Step 4: Identify proximal dates outside the buffer period
outside_buffer_indices = setdiff(all_indices, buffer_indices); % Exclude buffer points

% Find proximal points (within proximity to labeled dates)
proximal_indices = [];
for i = 1:length(outside_buffer_indices)
    distances_to_labeled = pdist2(nmds_coords(outside_buffer_indices(i), :), nmds_coords(labeled_indices, :));
    if min(distances_to_labeled) < proximity_threshold
        proximal_indices = [proximal_indices, outside_buffer_indices(i)];
    end
end

% Step 5: Extract the month-year from the dates for the criteria meeting points
month_year_labels = cellstr(datestr(filtered_datetime_values, 'mmm-yyyy')); % Extract month-year in 'MMM-YYYY' format

% Initialize separate handles and labels for the clusters and month-year
scatter_handles_criteria = [];
legend_labels_clusters = cell(1, optimal_clusters); % For clusters
legend_labels_criteria = {}; % For month-year points
used_colors = zeros(0, 3); % Track used colors for legend creation
used_month_years = {}; % Track used month-years to avoid duplicates in legend
valid_month_years = {}; % Track month-years that meet criteria

% Add cluster labels to the legend
for k = 1:optimal_clusters
    legend_labels_clusters{k} = ['Cluster ', num2str(k)];
end

% Step 6: Identify valid month-years that meet the criteria
for i = 1:length(filtered_datetime_values) - 7 + 1
    % Extract 7 days of data
    date_window = filtered_datetime_values(i:i + 7 - 1);
    cluster_window = cluster_labels(i:i + 7 - 1);
    
    % Skip the window if there is a gap of more than 7 days
    if max(diff(date_window)) > days(7)
        continue; % Skip this window due to large gap
    end
    
    % Check if 4 or more days are in the same cluster
    most_common_cluster = mode(cluster_window);
    count_in_cluster = sum(cluster_window == most_common_cluster);
    
    % If 4 or more days are in the same cluster
    if count_in_cluster >= 4
        % Find the points in this window that belong to the most common cluster
        cluster_indices_in_window = find(cluster_window == most_common_cluster);
        date_indices_in_cluster = i + cluster_indices_in_window - 1; % Convert to global indices
        
        % Check if these points meet the proximity condition
        distances_to_labeled_cluster_points = pdist2(nmds_coords(date_indices_in_cluster, :), nmds_coords(labeled_indices, :));
        min_distances = min(distances_to_labeled_cluster_points, [], 2);
        
        % If 4 or more points in the cluster meet the proximity threshold
        if sum(min_distances < proximity_threshold) >= 4
            % Add the month-year of these points to the list of valid month-years
            month_year_label = month_year_labels{date_indices_in_cluster(1)}; % Get the month-year of the first point
            if ~ismember(month_year_label, valid_month_years)
                valid_month_years{end+1} = month_year_label; % Add the month-year to the valid list
            end
        end
    end
end

% Now generate colors only for valid month-years
num_valid_month_years = length(valid_month_years);
colors_month_year = [
    0.13, 0.55, 0.49;   % Soft but dark blue-green
    0.79, 0.18, 0.57;   % Magenta
    0, 0, 0];            % Black

% Step 7: Plot and color only criteria meeting points by month-year
for i = 1:length(filtered_datetime_values) - 7 + 1
    % Extract 7 days of data
    date_window = filtered_datetime_values(i:i + 7 - 1);
    cluster_window = cluster_labels(i:i + 7 - 1);
    
    % Skip the window if there is a gap of more than 7 days
    if max(diff(date_window)) > days(7)
        continue; % Skip this window due to large gap
    end
    
    % Check if 4 or more days are in the same cluster
    most_common_cluster = mode(cluster_window);
    count_in_cluster = sum(cluster_window == most_common_cluster);
    
    % If 4 or more days are in the same cluster
    if count_in_cluster >= 4
        % Find the points in this window that belong to the most common cluster
        cluster_indices_in_window = find(cluster_window == most_common_cluster);
        date_indices_in_cluster = i + cluster_indices_in_window - 1; % Convert to global indices
        
        % Check if these points meet the proximity condition
        distances_to_labeled_cluster_points = pdist2(nmds_coords(date_indices_in_cluster, :), nmds_coords(labeled_indices, :));
        min_distances = min(distances_to_labeled_cluster_points, [], 2);
        
        % If 4 or more points in the cluster meet the proximity threshold
        if sum(min_distances < proximity_threshold) >= 4
            % Assign the color based on the month-year of the criteria meeting points
            month_year_label = month_year_labels{date_indices_in_cluster(1)}; % Get the month-year of the first point in the cluster
            month_year_idx = find(strcmp(valid_month_years, month_year_label)); % Find the corresponding index for color
            
            % Only create a scatter handle for the first occurrence of a new month-year
            if ~ismember(month_year_label, used_month_years)
                % Color the first criteria meeting points by their month-year and add to the legend
                ptsz=140;
                scatter_handle = scatter(nmds_coords(date_indices_in_cluster(1), 1), nmds_coords(date_indices_in_cluster(1), 2), ...
                                         ptsz, colors_month_year(month_year_idx, :), 'filled','diamond');
                scatter_handles_criteria = [scatter_handles_criteria, scatter_handle];
                
                % Add the legend entry for the first occurrence
                legend_labels_criteria{end+1} = month_year_label; % Add month-year to the legend
                used_month_years{end+1} = month_year_label; % Track the used month-year
            end
            
            % Plot the rest of the points for this month-year without adding to the legend
            for j = 2:length(date_indices_in_cluster)
                date_idx = date_indices_in_cluster(j);
                
                % Color the subsequent criteria meeting points by their month-year
                scatter(nmds_coords(date_idx, 1), nmds_coords(date_idx, 2), ptsz, colors_month_year(month_year_idx, :), 'filled','diamond');
            end
        end
    end
end

% Step 8: Label points of interest in black
black_color = [0 0 0]; % Color for points of interest
for i = labeled_indices
    scatter(nmds_coords(i, 1), nmds_coords(i, 2), 60, black_color, 'filled');
    
    % Label the points of interest
    text(nmds_coords(i, 1), nmds_coords(i, 2), datestr(filtered_datetime_values(i), 'dd-mmm-yyyy'), ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'FontSize', 16, 'Rotation', 45, 'Color', black_color);
end

% Step 9: Draw ellipses around clusters using proper cluster colors
for i = 1:optimal_clusters
    % Get the coordinates for points in the current cluster
    cluster_points = nmds_coords(cluster_labels == i, :);
    
    % Plot an ellipse around the cluster points using the cluster color
    plot_ellipse(cluster_points(:, 1), cluster_points(:, 2), cluster_centers(i, :), colors_cluster(i, :)); % Use cluster color for ellipse
end

% Step 10: Manually add legend for clusters and month-year points
full_legend_labels = [legend_labels_clusters, legend_labels_criteria]; % Combine cluster and month-year labels
scatter_handles = [scatter_handles_clusters', scatter_handles_criteria]; % Combine handles for clusters and month-year points
legend(scatter_handles, full_legend_labels, 'Location', 'bestoutside','AutoUpdate','Off'); % Add the legend after plotting everything

% Custom function to plot ellipses around clusters using cluster color
function plot_ellipse(x, y, center, color)
    theta = linspace(0, 2*pi, 100); % Array of angles
    x_mean = mean(x); % Mean of x-coordinates
    y_mean = mean(y); % Mean of y-coordinates
    x_radius = std(x) * 2; % Use 2 standard deviations for the radius
    y_radius = std(y) * 2;
    ellipse_x = x_mean + x_radius * cos(theta); % X points of the ellipse
    ellipse_y = y_mean + y_radius * sin(theta); % Y points of the ellipse
    
    % Plot the ellipse using the specified cluster color
    plot(ellipse_x, ellipse_y, 'Color', color, 'LineWidth', 3);
    hold on;
end

% Final plot adjustments
ylim([-0.4 0.5])
title('NMDS of 6-day Bray-Curtis Dissimilarity with similar anomalous shifts highlighted');
xlabel('NMDS Dimension 1');
ylabel('NMDS Dimension 2');
hold off;
set(gca,'FontSize',16)

set(gcf, 'Position', [0 0 1600 1200]);


saving=0;
if saving==1
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\bray_curtis\nmds_v0.png"]);
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\bray_curtis\nmds_v0.pdf"]);
end
saving=0;






%% II. Using only significant taxa


% Define significant taxa (as provided)
significant_taxa = {'Asterionellopsis', 'Centric', 'Hemiaulus', 'Leptocylindrus', 'Skeletonema', 'Thalassionema', 'Thalassiosira'};

% Step 1: Extract the 'datetime' column for later use
datetime_values = ifcbbray.datetime; % Assuming 'datetime' is a datetime column

% Step 2: Focus on columns 2 through 52 for taxa counts
taxa_data = ifcbbray(:, 2:52); % Only select columns 2 to 52, excluding datetime

% Step 3: Identify significant taxa within these columns
significant_indices = find(ismember(taxa_data.Properties.VariableNames, significant_taxa));

% Create a logical array to identify significant taxa in columns 2 to 52
is_significant = false(1, width(taxa_data));
is_significant(significant_indices) = true;

% Step 4: Sum the abundances of non-significant taxa in columns 2 to 52
other_taxa_abundance = sum(table2array(taxa_data(:, ~is_significant)), 2);

% Step 5: Create a new table with significant taxa and 'Other Taxa' column
new_taxa_abundances = [taxa_data(:, is_significant), array2table(other_taxa_abundance)]; % Do not include datetime yet

% Update variable names to include 'Other Taxa'
new_variable_names = [taxa_data.Properties.VariableNames(is_significant), 'Other Taxa'];
new_taxa_abundances.Properties.VariableNames = new_variable_names;

% Step 6: Extract the taxa counts (significant + Other Taxa)
taxa_counts = new_taxa_abundances{:,:}; % Extract the table as an array for further analysis

% Step 7: Identify samples that have another sample exactly 6 days apart
valid_indices = false(size(datetime_values)); % Preallocate a logical array

for i = 1:length(datetime_values)
    % Look for another sample exactly 6 days before or after
    match_before = find(datetime_values == datetime_values(i) - days(6), 1);
    match_after = find(datetime_values == datetime_values(i) + days(6), 1);
    
    if ~isempty(match_before) || ~isempty(match_after)
        valid_indices(i) = true; % Mark this sample as valid if it has a match
    end
end

% Step 8: Filter the taxa counts and datetime values based on valid indices
filtered_taxa_counts = taxa_counts(valid_indices, :);
filtered_datetime_values = datetime_values(valid_indices);


bray_dist = braycurtis_batch(filtered_taxa_counts);

% Step 10: Perform NMDS using 'mdscale' on Bray-Curtis dissimilarities
nmds_dim = 3; % For 2D NMDS

% Convert the condensed distance vector into a full distance matrix
square_brays_dist = squareform(bray_dist); % Converts vector to square matrix

% Perform NMDS
[nmds_coords, stress] = mdscale(square_brays_dist, nmds_dim, 'Criterion', 'stress');

%% Step 1: Perform k-means clustering with different numbers of clusters on NMDS coordinates
max_clusters = 10; % Maximum number of clusters to test
silhouette_avg = zeros(max_clusters-1, 1); % Store average silhouette scores

for k = 2:max_clusters
    % Apply k-means clustering
    cluster_labels = kmeans(nmds_coords, k, 'Replicates', 10, 'Distance', 'sqeuclidean');
    
    % Calculate silhouette scores for current number of clusters
    silhouette_avg(k-1) = mean(silhouette(nmds_coords, cluster_labels)); % Silhouette score
end


%% Step 2: Plot the silhouette scores to find the optimal number of clusters
figure;
plot(2:max_clusters, silhouette_avg, '-o');
title('Silhouette Scores for Different Numbers of Clusters');
xlabel('Number of Clusters');
ylabel('Average Silhouette Score');

% Step 3: Use the elbow method to find the optimal number of clusters
sum_of_squared_distances = zeros(max_clusters-1, 1);
for k = 2:max_clusters
    % Apply k-means clustering
    [~, ~, sumd] = kmeans(nmds_coords, k, 'Replicates', 10, 'Distance', 'sqeuclidean');
    sum_of_squared_distances(k-1) = sum(sumd); % Total within-cluster sum of squares
end

% Plot SSE (Elbow Method)
figure;
plot(2:max_clusters, sum_of_squared_distances, '-o');
title('Elbow Method for Optimal Number of Clusters');
xlabel('Number of Clusters');
ylabel('Sum of Squared Distances');

% Step 4: Based on the silhouette and elbow method plots, select the optimal number of clusters
optimal_clusters = 4; % Example: select optimal number of clusters (adjust based on results)

% Step 5: Apply k-means with the optimal number of clusters
[cluster_labels, cluster_centers] = kmeans(nmds_coords, optimal_clusters, 'Replicates', 10, 'Distance', 'sqeuclidean');

%% Step 6: Plot the NMDS results and color by cluster
close all
% Define updated date range for labeling points of interest
label_start_date = datetime(2020, 8, 21);
label_end_date = datetime(2020, 8, 2);

% Define 30-day buffer period
buffer_start_date = label_start_date - days(30);
buffer_end_date = label_end_date + days(30);

% Proximity threshold for NMDS space
proximity_threshold = 0.08;

% Step 1: Identify points in the given date range and 30-day buffer
labeled_indices = find(filtered_datetime_values >= label_start_date & filtered_datetime_values <= label_end_date);
buffer_indices = find(filtered_datetime_values >= buffer_start_date & filtered_datetime_values <= buffer_end_date);

% Step 2: Perform k-means clustering
optimal_clusters = 4; % Number of clusters from elbow method
[cluster_labels, cluster_centers] = kmeans(nmds_coords, optimal_clusters, 'Replicates', 10, 'Distance', 'sqeuclidean');

% Step 3: Plot NMDS points - color all points by cluster using consistent colors
figure;
all_indices = 1:length(filtered_datetime_values);
colors_cluster = lines(optimal_clusters); % Use consistent colormap for clusters

% Plot all points colored by cluster using consistent cluster colors
scatter_handles_clusters = gscatter(nmds_coords(:, 1), nmds_coords(:, 2), cluster_labels, colors_cluster, '.', 20);
hold on;

% Step 4: Identify proximal dates outside the buffer period
outside_buffer_indices = setdiff(all_indices, buffer_indices); % Exclude buffer points

% Find proximal points (within proximity to labeled dates)
proximal_indices = [];
for i = 1:length(outside_buffer_indices)
    distances_to_labeled = pdist2(nmds_coords(outside_buffer_indices(i), :), nmds_coords(labeled_indices, :));
    if min(distances_to_labeled) < proximity_threshold
        proximal_indices = [proximal_indices, outside_buffer_indices(i)];
    end
end

% Step 5: Extract the month-year from the dates for the criteria meeting points
month_year_labels = cellstr(datestr(filtered_datetime_values, 'mmm-yyyy')); % Extract month-year in 'MMM-YYYY' format

% Initialize separate handles and labels for the clusters and month-year
scatter_handles_criteria = [];
legend_labels_clusters = cell(1, optimal_clusters); % For clusters
legend_labels_criteria = {}; % For month-year points
used_colors = zeros(0, 3); % Track used colors for legend creation
used_month_years = {}; % Track used month-years to avoid duplicates in legend
valid_month_years = {}; % Track month-years that meet criteria

% Add cluster labels to the legend
for k = 1:optimal_clusters
    legend_labels_clusters{k} = ['Cluster ', num2str(k)];
end

% Step 6: Identify valid month-years that meet the criteria
for i = 1:length(filtered_datetime_values) - 7 + 1
    % Extract 7 days of data
    date_window = filtered_datetime_values(i:i + 7 - 1);
    cluster_window = cluster_labels(i:i + 7 - 1);
    
    % Skip the window if there is a gap of more than 7 days
    if max(diff(date_window)) > days(7)
        continue; % Skip this window due to large gap
    end
    
    % Check if 4 or more days are in the same cluster
    most_common_cluster = mode(cluster_window);
    count_in_cluster = sum(cluster_window == most_common_cluster);
    
    % If 4 or more days are in the same cluster
    if count_in_cluster >= 4
        % Find the points in this window that belong to the most common cluster
        cluster_indices_in_window = find(cluster_window == most_common_cluster);
        date_indices_in_cluster = i + cluster_indices_in_window - 1; % Convert to global indices
        
        % Check if these points meet the proximity condition
        distances_to_labeled_cluster_points = pdist2(nmds_coords(date_indices_in_cluster, :), nmds_coords(labeled_indices, :));
        min_distances = min(distances_to_labeled_cluster_points, [], 2);
        
        % If 4 or more points in the cluster meet the proximity threshold
        if sum(min_distances < proximity_threshold) >= 4
            % Add the month-year of these points to the list of valid month-years
            month_year_label = month_year_labels{date_indices_in_cluster(1)}; % Get the month-year of the first point
            if ~ismember(month_year_label, valid_month_years)
                valid_month_years{end+1} = month_year_label; % Add the month-year to the valid list
            end
        end
    end
end

% Now generate colors only for valid month-years
num_valid_month_years = length(valid_month_years);
colors_month_year = [
    0.13, 0.55, 0.49;   % Soft but dark blue-green (original)
    0.79, 0.18, 0.57;   % Magenta (original)
    0.00, 0.00, 0.00;   % Black (original)
    0.98, 0.68, 0.13;   % Orange-gold for high contrast
    0.18, 0.33, 0.66;   % Blue with medium contrast
    0.85, 0.33, 0.10;   % Dark orange for warm contrast
    0.57, 0.44, 0.86;   % Soft lavender for light contrast
    0.39, 0.76, 0.65;   % Aqua-green for soft contrast
    0.70, 0.11, 0.11;   % Dark red for strong contrast
    0.31, 0.50, 0.74;   % Muted blue for a neutral contrast
];
          % Black

% Step 7: Plot and color only criteria meeting points by month-year
for i = 1:length(filtered_datetime_values) - 7 + 1
    % Extract 7 days of data
    date_window = filtered_datetime_values(i:i + 7 - 1);
    cluster_window = cluster_labels(i:i + 7 - 1);
    
    % Skip the window if there is a gap of more than 7 days
    if max(diff(date_window)) > days(7)
        continue; % Skip this window due to large gap
    end
    
    % Check if 4 or more days are in the same cluster
    most_common_cluster = mode(cluster_window);
    count_in_cluster = sum(cluster_window == most_common_cluster);
    
    % If 4 or more days are in the same cluster
    if count_in_cluster >= 4
        % Find the points in this window that belong to the most common cluster
        cluster_indices_in_window = find(cluster_window == most_common_cluster);
        date_indices_in_cluster = i + cluster_indices_in_window - 1; % Convert to global indices
        
        % Check if these points meet the proximity condition
        distances_to_labeled_cluster_points = pdist2(nmds_coords(date_indices_in_cluster, :), nmds_coords(labeled_indices, :));
        min_distances = min(distances_to_labeled_cluster_points, [], 2);
        
        % If 4 or more points in the cluster meet the proximity threshold
        if sum(min_distances < proximity_threshold) >= 4
            % Assign the color based on the month-year of the criteria meeting points
            month_year_label = month_year_labels{date_indices_in_cluster(1)}; % Get the month-year of the first point in the cluster
            month_year_idx = find(strcmp(valid_month_years, month_year_label)); % Find the corresponding index for color
            
            % Only create a scatter handle for the first occurrence of a new month-year
            if ~ismember(month_year_label, used_month_years)
                % Color the first criteria meeting points by their month-year and add to the legend
                ptsz=140;
                scatter_handle = scatter(nmds_coords(date_indices_in_cluster(1), 1), nmds_coords(date_indices_in_cluster(1), 2), ...
                                         ptsz, colors_month_year(month_year_idx, :), 'filled','diamond');
                scatter_handles_criteria = [scatter_handles_criteria, scatter_handle];
                
                % Add the legend entry for the first occurrence
                legend_labels_criteria{end+1} = month_year_label; % Add month-year to the legend
                used_month_years{end+1} = month_year_label; % Track the used month-year
            end
            
            % Plot the rest of the points for this month-year without adding to the legend
            for j = 2:length(date_indices_in_cluster)
                date_idx = date_indices_in_cluster(j);
                
                % Color the subsequent criteria meeting points by their month-year
                scatter(nmds_coords(date_idx, 1), nmds_coords(date_idx, 2), ptsz, colors_month_year(month_year_idx, :), 'filled','diamond');
            end
        end
    end
end

% Step 8: Label points of interest in black
black_color = [0 0 0]; % Color for points of interest
for i = labeled_indices
    scatter(nmds_coords(i, 1), nmds_coords(i, 2), 60, black_color, 'filled');
    
    % Label the points of interest
    text(nmds_coords(i, 1), nmds_coords(i, 2), datestr(filtered_datetime_values(i), 'dd-mmm-yyyy'), ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'FontSize', 16, 'Rotation', 45, 'Color', black_color);
end

% Step 9: Draw ellipses around clusters using proper cluster colors
for i = 1:optimal_clusters
    % Get the coordinates for points in the current cluster
    cluster_points = nmds_coords(cluster_labels == i, :);
    
    % Plot an ellipse around the cluster points using the cluster color
    plot_ellipse(cluster_points(:, 1), cluster_points(:, 2), cluster_centers(i, :), colors_cluster(i, :)); % Use cluster color for ellipse
end

% Step 10: Manually add legend for clusters and month-year points
full_legend_labels = [legend_labels_clusters, legend_labels_criteria]; % Combine cluster and month-year labels
scatter_handles = [scatter_handles_clusters', scatter_handles_criteria]; % Combine handles for clusters and month-year points
legend(scatter_handles, full_legend_labels, 'Location', 'bestoutside','AutoUpdate','Off'); % Add the legend after plotting everything



%% Step 1: Perform k-means clustering with different numbers of clusters on NMDS coordinates
max_clusters = 10; % Maximum number of clusters to test
silhouette_avg = zeros(max_clusters-1, 1); % Store average silhouette scores

for k = 2:max_clusters
    % Apply k-means clustering
    cluster_labels = kmeans(nmds_coords, k, 'Replicates', 10, 'Distance', 'sqeuclidean');
    
    % Calculate silhouette scores for current number of clusters
    silhouette_avg(k-1) = mean(silhouette(nmds_coords, cluster_labels)); % Silhouette score
end

%% Step 2: Plot the silhouette scores to find the optimal number of clusters
figure;
plot(2:max_clusters, silhouette_avg, '-o');
title('Silhouette Scores for Different Numbers of Clusters');
xlabel('Number of Clusters');
ylabel('Average Silhouette Score');

% Step 3: Use the elbow method to find the optimal number of clusters
sum_of_squared_distances = zeros(max_clusters-1, 1);
for k = 2:max_clusters
    % Apply k-means clustering
    [~, ~, sumd] = kmeans(nmds_coords, k, 'Replicates', 10, 'Distance', 'sqeuclidean');
    sum_of_squared_distances(k-1) = sum(sumd); % Total within-cluster sum of squares
end

% Plot SSE (Elbow Method)
figure;
plot(2:max_clusters, sum_of_squared_distances, '-o');
title('Elbow Method for Optimal Number of Clusters');
xlabel('Number of Clusters');
ylabel('Sum of Squared Distances');

% Step 4: Based on the silhouette and elbow method plots, select the optimal number of clusters
optimal_clusters = 5; % Example: select optimal number of clusters (adjust based on results)

% Step 5: Apply k-means with the optimal number of clusters
[cluster_labels, cluster_centers] = kmeans(nmds_coords, optimal_clusters, 'Replicates', 10, 'Distance', 'sqeuclidean');

%% Step 6: Plot the NMDS results and color by cluster
close all
% Define updated date range for labeling points of interest
label_start_date = datetime(2020, 8, 21);
label_end_date = datetime(2020, 9, 2);

% Define 30-day buffer period
buffer_start_date = label_start_date - days(30);
buffer_end_date = label_end_date + days(30);

% Proximity threshold for NMDS space
proximity_threshold = 0.1;

% Step 1: Identify points in the given date range and 30-day buffer
labeled_indices = find(filtered_datetime_values >= label_start_date & filtered_datetime_values <= label_end_date);
buffer_indices = find(filtered_datetime_values >= buffer_start_date & filtered_datetime_values <= buffer_end_date);

% Step 2: Perform k-means clustering
optimal_clusters = 4; % Number of clusters from elbow method
[cluster_labels, cluster_centers] = kmeans(nmds_coords, optimal_clusters, 'Replicates', 10, 'Distance', 'sqeuclidean');

% Step 3: Plot NMDS points - color all points by cluster using consistent colors
figure;
all_indices = 1:length(filtered_datetime_values);
colors_cluster = lines(optimal_clusters); % Use consistent colormap for clusters

% Plot all points colored by cluster using consistent cluster colors
scatter_handles_clusters = gscatter(nmds_coords(:, 1), nmds_coords(:, 2), cluster_labels, colors_cluster, '.', 20);
hold on;

% Step 4: Identify proximal dates outside the buffer period
outside_buffer_indices = setdiff(all_indices, buffer_indices); % Exclude buffer points

% Find proximal points (within proximity to labeled dates)
proximal_indices = [];
for i = 1:length(outside_buffer_indices)
    distances_to_labeled = pdist2(nmds_coords(outside_buffer_indices(i), :), nmds_coords(labeled_indices, :));
    if min(distances_to_labeled) < proximity_threshold
        proximal_indices = [proximal_indices, outside_buffer_indices(i)];
    end
end

% Step 5: Extract the month-year from the dates for the criteria meeting points
month_year_labels = cellstr(datestr(filtered_datetime_values, 'mmm-yyyy')); % Extract month-year in 'MMM-YYYY' format

% Initialize separate handles and labels for the clusters and month-year
scatter_handles_criteria = [];
legend_labels_clusters = cell(1, optimal_clusters); % For clusters
legend_labels_criteria = {}; % For month-year points
used_colors = zeros(0, 3); % Track used colors for legend creation
used_month_years = {}; % Track used month-years to avoid duplicates in legend
valid_month_years = {}; % Track month-years that meet criteria

% Add cluster labels to the legend
for k = 1:optimal_clusters
    legend_labels_clusters{k} = ['Cluster ', num2str(k)];
end

% Step 6: Identify valid month-years that meet the criteria
for i = 1:length(filtered_datetime_values) - 7 + 1
    % Extract 7 days of data
    date_window = filtered_datetime_values(i:i + 7 - 1);
    cluster_window = cluster_labels(i:i + 7 - 1);
    
    % Skip the window if there is a gap of more than 7 days
    if max(diff(date_window)) > days(7)
        continue; % Skip this window due to large gap
    end
    
    % Check if 4 or more days are in the same cluster
    most_common_cluster = mode(cluster_window);
    count_in_cluster = sum(cluster_window == most_common_cluster);
    
    % If 4 or more days are in the same cluster
    if count_in_cluster >= 4
        % Find the points in this window that belong to the most common cluster
        cluster_indices_in_window = find(cluster_window == most_common_cluster);
        date_indices_in_cluster = i + cluster_indices_in_window - 1; % Convert to global indices
        
        % Check if these points meet the proximity condition
        distances_to_labeled_cluster_points = pdist2(nmds_coords(date_indices_in_cluster, :), nmds_coords(labeled_indices, :));
        min_distances = min(distances_to_labeled_cluster_points, [], 2);
        
        % If 4 or more points in the cluster meet the proximity threshold
        if sum(min_distances < proximity_threshold) >= 4
            % Add the month-year of these points to the list of valid month-years
            month_year_label = month_year_labels{date_indices_in_cluster(1)}; % Get the month-year of the first point
            if ~ismember(month_year_label, valid_month_years)
                valid_month_years{end+1} = month_year_label; % Add the month-year to the valid list
            end
        end
    end
end

% Now generate colors only for valid month-years
num_valid_month_years = length(valid_month_years);

% Step 7: Plot and color only criteria meeting points by month-year
for i = 1:length(filtered_datetime_values) - 7 + 1
    % Extract 7 days of data
    date_window = filtered_datetime_values(i:i + 7 - 1);
    cluster_window = cluster_labels(i:i + 7 - 1);
    
    % Skip the window if there is a gap of more than 7 days
    if max(diff(date_window)) > days(7)
        continue; % Skip this window due to large gap
    end
    
    % Check if 4 or more days are in the same cluster
    most_common_cluster = mode(cluster_window);
    count_in_cluster = sum(cluster_window == most_common_cluster);
    
    % If 4 or more days are in the same cluster
    if count_in_cluster >= 4
        % Find the points in this window that belong to the most common cluster
        cluster_indices_in_window = find(cluster_window == most_common_cluster);
        date_indices_in_cluster = i + cluster_indices_in_window - 1; % Convert to global indices
        
        % Check if these points meet the proximity condition
        distances_to_labeled_cluster_points = pdist2(nmds_coords(date_indices_in_cluster, :), nmds_coords(labeled_indices, :));
        min_distances = min(distances_to_labeled_cluster_points, [], 2);
        
        % If 4 or more points in the cluster meet the proximity threshold
        if sum(min_distances < proximity_threshold) >= 4
            % Assign the color based on the month-year of the criteria meeting points
            month_year_label = month_year_labels{date_indices_in_cluster(1)}; % Get the month-year of the first point in the cluster
            month_year_idx = find(strcmp(valid_month_years, month_year_label)); % Find the corresponding index for color
            
            % Only create a scatter handle for the first occurrence of a new month-year
            if ~ismember(month_year_label, used_month_years)
                % Color the first criteria meeting points by their month-year and add to the legend
                ptsz=140;
                scatter_handle = scatter(nmds_coords(date_indices_in_cluster(1), 1), nmds_coords(date_indices_in_cluster(1), 2), ...
                                         ptsz, colors_month_year(month_year_idx, :), 'filled','diamond');
                scatter_handles_criteria = [scatter_handles_criteria, scatter_handle];
                
                % Add the legend entry for the first occurrence
                legend_labels_criteria{end+1} = month_year_label; % Add month-year to the legend
                used_month_years{end+1} = month_year_label; % Track the used month-year
            end
            
            % Plot the rest of the points for this month-year without adding to the legend
            for j = 2:length(date_indices_in_cluster)
                date_idx = date_indices_in_cluster(j);
                
                % Color the subsequent criteria meeting points by their month-year
                scatter(nmds_coords(date_idx, 1), nmds_coords(date_idx, 2), ptsz, colors_month_year(month_year_idx, :), 'filled','diamond');
            end
        end
    end
end

% Step 8: Label points of interest in black
black_color = [0 0 0]; % Color for points of interest
for i = labeled_indices
    scatter(nmds_coords(i, 1), nmds_coords(i, 2), 60, black_color, 'filled');
    
    % Label the points of interest
    text(nmds_coords(i, 1), nmds_coords(i, 2), datestr(filtered_datetime_values(i), 'dd-mmm-yyyy'), ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'FontSize', 16, 'Rotation', 45, 'Color', black_color);
end

% Step 9: Draw ellipses around clusters using proper cluster colors
for i = 1:optimal_clusters
    % Get the coordinates for points in the current cluster
    cluster_points = nmds_coords(cluster_labels == i, :);
    
    % Plot an ellipse around the cluster points using the cluster color
    plot_ellipse(cluster_points(:, 1), cluster_points(:, 2), cluster_centers(i, :), colors_cluster(i, :)); % Use cluster color for ellipse
end

% Step 10: Manually add legend for clusters and month-year points
full_legend_labels = [legend_labels_clusters, legend_labels_criteria]; % Combine cluster and month-year labels
scatter_handles = [scatter_handles_clusters', scatter_handles_criteria]; % Combine handles for clusters and month-year points
legend(scatter_handles, full_legend_labels, 'Location', 'bestoutside','AutoUpdate','Off'); % Add the legend after plotting everything



%% Now not using NMDS but instead assembalge occurrence

% Define significant taxa (as provided)
significant_taxa = {'Asterionellopsis', 'Centric', 'Hemiaulus', 'Leptocylindrus', 'Skeletonema', 'Thalassionema', 'Thalassiosira'};

% Step 1: Extract the datetime and taxa columns 2 to 52
datetime_values = ifcbbray.datetime;
taxa_data = ifcbbray{:, 2:52}; % Extract columns 2 to 52 (the whole community)

% Step 2: Normalize the entire community abundance for each date
total_abundance = sum(taxa_data, 2);  % Sum of all taxa for each date
relative_abundance_data = taxa_data ./ total_abundance;  % Relative abundance for each date

% Step 3: Identify the significant taxa columns within the taxa_data (columns 2:52)
significant_indices = find(ismember(ifcbbray.Properties.VariableNames(2:52), significant_taxa));
relative_abundance_significant_taxa = relative_abundance_data(:, significant_indices);  % Extract relative abundance for significant taxa

% Step 4: Define the target event date range (8/21/2020 - 8/31/2020)
target_start_date = datetime(2020, 8, 20);
target_end_date = datetime(2020, 9, 2);
event_indices = find(datetime_values >= target_start_date & datetime_values <= target_end_date);

% Step 5: Compute the average relative abundance for the target event window (across all target dates)
target_event = mean(relative_abundance_significant_taxa(event_indices, :), 1);

% Step 6: Calculate Bray-Curtis dissimilarity between each date's significant taxa assemblage and the target event
n_samples = size(relative_abundance_significant_taxa, 1);
bray_curtis_dissimilarity = zeros(n_samples, 1);

for i = 1:n_samples
    numerator = sum(abs(relative_abundance_significant_taxa(i, :) - target_event));
    denominator = sum(relative_abundance_significant_taxa(i, :) + target_event);
    bray_curtis_dissimilarity(i) = numerator / denominator;
end

% Step 7: Set a threshold for high similarity (low Bray-Curtis dissimilarity)
threshold = prctile(bray_curtis_dissimilarity, 5); % 10th percentile as threshold for high similarity

% Step 8: Find dates where the dissimilarity is below the threshold
similarity_indices = find(bray_curtis_dissimilarity <= threshold);
similar_dates = datetime_values(similarity_indices);

% Step 9: Identify periods where at least 4/7 days meet the threshold
selected_dates = []; % Store final dates that meet the 4/7 day criteria

min_days = 4; % Minimum number of days in a 7-day window that must meet the similarity threshold
max_gap = 7;  % Maximum allowable gap between consecutive dates

% Step 10: Loop through similarity_indices to check for windows of 7 days with at least 4 days that meet criteria
for i = 1:length(similarity_indices)
    % Create a 7-day moving window
    for j = i:min(i+6, length(similarity_indices))
        window = similarity_indices(i:j); % Current window of indices
        
        % Calculate the gaps between consecutive dates in this window
        date_diffs = days(diff(datetime_values(window)));
        
        % Check if the maximum gap between consecutive days is <= max_gap
        if all(date_diffs <= max_gap)
            % Check if the window has at least 4 days that meet the similarity threshold
            if length(window) >= min_days
                selected_dates = [selected_dates; datetime_values(window)];
            end
        end
    end
end

% Remove duplicates from selected_dates
selected_dates = unique(selected_dates);

% Step 11: Display the dates that meet the continuous criteria
disp('Dates with continuous reappearance of the taxa assemblage:');
disp(selected_dates);


%% Initialize variables to store unique periods
periods = [];
mean_percentiles = [];

% Step 1: Calculate the percentile of each Bray-Curtis dissimilarity score
bray_percentiles = tiedrank(bray_curtis_dissimilarity) / length(bray_curtis_dissimilarity) * 100;

% Step 2: Loop through the selected_dates to identify continuous periods
start_idx = 1;
for i = 2:length(selected_dates)
    % Check if the gap between consecutive dates is more than 7 days
    if days(selected_dates(i) - selected_dates(i-1)) > 7
        % End the current period and compute the mean percentile for the current period
        end_idx = i - 1;
        current_period = selected_dates(start_idx:end_idx);
        
        % Find indices corresponding to the current period in datetime_values
        period_indices = ismember(datetime_values, current_period);
        
        % Calculate the mean percentile for the current period
        mean_percentile = mean(bray_percentiles(period_indices));
        
        % Store the start date, end date, and mean percentile
        periods = [periods; selected_dates(start_idx), selected_dates(end_idx)];
        mean_percentiles = [mean_percentiles; mean_percentile];
        
        % Update the start index for the next period
        start_idx = i;
    end
end

% Add the final period if there's a remaining continuous segment
if start_idx <= length(selected_dates)
    current_period = selected_dates(start_idx:end);
    period_indices = ismember(datetime_values, current_period);
    mean_percentile = mean(bray_percentiles(period_indices));
    periods = [periods; selected_dates(start_idx), selected_dates(end)];
    mean_percentiles = [mean_percentiles; mean_percentile];
end

% Step 3: Create a 3-column table with start date, end date, and mean Bray-Curtis percentile
result_table = table(periods(:,1), periods(:,2), mean_percentiles, ...
                     'VariableNames', {'StartDate', 'EndDate', 'MeanPercentile'});

% Display the result table
disp(result_table);


%% Functions



