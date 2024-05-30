% Load your data (adjust according to your data structure and format)
load('binnedCUGN66_1ce6_95ff_4953.mat', 'binnedCUGN66')



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

clf 
plotOceanData(cugn_sel, 'salinity', 'depth', 'longitude')

%%

function plotOceanData(data, variableName, depthName, distanceName)
    % This function plots a contour plot of the specified variable against depth and distance.
    % 
    % Inputs:
    %   data - The table containing oceanographic data
    %   variableName - The name of the variable to plot (e.g., 'temperature', 'salinity', 'u')
    %   depthName - The name of the column in data that contains depth information
    %   distanceName - The name of the column in data to calculate offshore distance

    % Convert longitude to offshore distance (simple linear transformation, could be improved)
    minLongitude = min(data.(distanceName));
    distance = (data.(distanceName) - minLongitude) * -111;  % km/degree conversion

    % Extract the depth and the variable
    depth = data.(depthName);
    variable = data.(variableName);

    % Create the meshgrid
    [DistanceGrid, DepthGrid] = meshgrid(unique(distance), unique(depth));

    % Grid the variable data using griddata for interpolation
    VariableGrid = griddata(distance, depth, variable, DistanceGrid, DepthGrid);

    % Plotting
    figure
    contourf(DistanceGrid, DepthGrid, VariableGrid, 20, 'LineColor', 'none');
    colormap(jet);
    colorbar;
    title(['Distribution of ', variableName]);
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
