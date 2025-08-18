% Load your data (adjust according to your data structure and format)
load('Q:\Dante\data\MB_Wildfire_Obs\gliders\binnedCUGN66_1ce6_95ff_4953.mat', 'binnedCUGN66')



%% Convert to daily TT
binnedCUGN66.datetime=datetime(binnedCUGN66.time, 'ConvertFrom', 'posixtime');
binnedCUGN66_t=struct2table(binnedCUGN66);
binnedCUGN66_tt=table2timetable(binnedCUGN66_t);
binnedCUGN66_tt.profile=double(binnedCUGN66_tt.profile);
binnedCUGN66_tt.mission=double(binnedCUGN66_tt.mission);

binnedCUGN66_tt_daily=retime(binnedCUGN66_tt,'daily','mean');
%%
cugn_sel=binnedCUGN66_t(binnedCUGN66_t.datetime.Year==2020 & binnedCUGN66_t.datetime.Month==8,:);

%%
u=cugn_sel.u;  
depth=cugn_sel.depth;

%%

% Convert longitude to offshore distance if necessary
% Assuming that the distance is a simple linear transformation of longitude:
minLongitude = min(cugn_sel.longitude);
distance = (cugn_sel.longitude - minLongitude) * -111; % Approx conversion factor for km/degree at mid-latitudes

% Depth and geostrophic velocity (mean of u and v for simplicity or use one of them)
depth = cugn_sel.depth;
u_velocity = mean([cugn_sel.u, cugn_sel.v], 2);  % Mean velocity or choose one

% Create the meshgrid (assuming distance and depth are already suitable for this)
[DistanceGrid, DepthGrid] = meshgrid(unique(distance), unique(depth));

% Grid the velocity data using griddata for interpolation
VelocityGrid = griddata(distance, depth, u_velocity, DistanceGrid, DepthGrid);

% Plotting
figure;
contourf(DistanceGrid, DepthGrid, VelocityGrid, 20, 'LineColor', 'none');  % More levels for smoother gradient
colormap(jet);
colorbar;
caxis([-0.1 0.1]);  % adjust based on actual range of your data
title('Geostrophic Velocity Distribution');
xlabel('Offshore Distance (km)');
ylabel('Depth (m)');
set(gca, 'YDir', 'reverse');

% Adding density or temperature contours if available
hold on;
[C, h] = contour(DistanceGrid, DepthGrid, griddata(distance, depth, cugn_sel.temperature, DistanceGrid, DepthGrid), [26 27], 'LineColor', 'k');
clabel(C, h);


%%

%% Filter for August 2020
cugn_sel = binnedCUGN66_t(binnedCUGN66_t.datetime.Year == 2020 & binnedCUGN66_t.datetime.Month == 8, :);

%% Extract unique days from the datetime
unique_days = unique(day(cugn_sel.datetime));  % Get unique days in August

%% Loop through each day and call the modified plotOceanData function
for i = 1:length(unique_days)
    current_day = unique_days(i);
    
    % Filter data for the current day
    day_sel = cugn_sel(day(cugn_sel.datetime) == current_day, :);
    
    % Plot using the modified function
    figure(1);
    plotOceanData(day_sel, 'salinity', 'depth', 'longitude', 'latitude');  % Now using both latitude and longitude
    pause(0.5)
    % Set title to include the specific day
    title(['Geostrophic Velocity Distribution for August ', num2str(current_day), ' (Within 10km of Shore)']);
end

%% Modified plotOceanData function
function plotOceanData(data, variableName, depthName, longitudeName, latitudeName)
    % This function plots a contour plot of the specified variable against depth and distance from the shore.
    % It calculates the offshore distance from the fixed point (36.8874425, -122.0059075).
    %
    % Inputs:
    %   data - The table containing oceanographic data
    %   variableName - The name of the variable to plot (e.g., 'temperature', 'salinity', 'u')
    %   depthName - The name of the column in data that contains depth information
    %   longitudeName - The name of the column in data for longitude
    %   latitudeName - The name of the column in data for latitude

    % Fixed point coordinates (latitude and longitude) from which the distance will be calculated
    ref_lat = 36.8874425;
    ref_lon = -122.0059075;

    % Calculate great-circle distance from the fixed point (haversine formula)
    R = 6371; % Radius of the Earth in km
    lat1 = deg2rad(data.(latitudeName));
    lon1 = deg2rad(data.(longitudeName));
    lat2 = deg2rad(ref_lat);
    lon2 = deg2rad(ref_lon);
    
    dlat = lat1 - lat2;
    dlon = lon1 - lon2;
    
    a = sin(dlat / 2).^2 + cos(lat1) .* cos(lat2) .* sin(dlon / 2).^2;
    c = 2 * atan2(sqrt(a), sqrt(1 - a));
    
    distance = R * c;  % Distance in km from the fixed reference point
    
    % Filter data to only include points within 10km offshore
    within_10km_idx = distance <= 10;
    data_within_10km = data(within_10km_idx, :);
    
    % Ensure that there is data left after filtering
    if isempty(data_within_10km)
        disp('No data within 10km of shore.');
        return;
    end

    % Extract the depth and the variable
    depth = data_within_10km.(depthName);
    variable = data_within_10km.(variableName);
    distance = distance(within_10km_idx);  % Filtered distances
    
    % Create the meshgrid
    [DistanceGrid, DepthGrid] = meshgrid(unique(distance), unique(depth));

    % Grid the variable data using griddata for interpolation
    VariableGrid = griddata(distance, depth, variable, DistanceGrid, DepthGrid);

    % Plotting
    contourf(DistanceGrid, DepthGrid, VariableGrid, 20, 'LineColor', 'none');
    colormap(jet);
    colorbar;
    title(['Distribution of ', variableName, ' within 10km of shore']);
    xlabel('Offshore Distance (km)');
    ylabel('Depth (m)');
    set(gca, 'YDir', 'reverse');

    % Check if temperature contours should be added
    if contains(variableName, {'temperature', 'salinity'})
        hold on;
        [C, h] = contour(DistanceGrid, DepthGrid, VariableGrid, 'LineColor', 'k');
        clabel(C, h);
    end
end