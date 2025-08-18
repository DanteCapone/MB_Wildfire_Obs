%MATLAB Script to analyze patterns in IFCB Phytoplankton abundances in
%association wildfire smoke and with other drivers

%Load in the data

%Shore Station
load('mlm_sst.mat')

%IFCB
load('ifcb_ucsc_all_mean.mat')
%Aeronet
load('monterey_lvl_1.5.mat')
load('monterey_lvl_2.0.mat')
ind=Monterey_lvl_1_5.AOD_500nm == -999.0;
Monterey_lvl_1_5.AOD_500nm(ind)=nan;



% 
% %% Compile pm2.5 into one times series
% 
% % Directory containing the data files
% dataDir = 'Q:\Dante\data\MB_Wildfire_Obs\pm2_5';
% 
% % Years to process
% years = 2008:2023;
% 
% % Initialize an empty table
% sc_pm2_5_daily = table();
% 
% % Loop through each year and read in the data files
% for year = years
%     % Generate the filename
%     filename = fullfile(dataDir, sprintf('sc_pm2.5_daily_%d.csv', year));
% 
%     % Read the CSV file
%     yearlyData = readtable(filename);
% 
%     % Concatenate the data
%     sc_pm2_5_daily = [sc_pm2_5_daily; yearlyData];
% end
% 
% sc_pm2_5_daily.datetime=sc_pm2_5_daily.Date;
% sc_pm2_5_daily.date=sc_pm2_5_daily.Date;
% sc_pm2_5_daily.Date=[];
% 
% sc_pm2_5_daily=sc_pm2_5_daily(sc_pm2_5_daily.LocalSiteName=="San Lorenzo Valley Middle School",:);
% sc_pm2_5_daily.pm2_5=sc_pm2_5_daily.pm2_5;
% sc_pm2_5_daily.pm2_5=[];
% % Save the concatenated table as a .mat file
% outputFilename = fullfile(dataDir, 'sc_pm2.5_daily.mat');
% % outputFilename = fullfile(dataDir, 'sc_pm2.5_daily_all.mat');
% 
% save(outputFilename, 'sc_pm2_5_daily');




%%
%% Load data
load('ifcb_bray.mat')



% Corss correlation between phytoplankton and AOD
plotting=1;
% ifcbbray=ifcbbray(ifcbbray.dateTime.Year==2020,:);
%Add AOD anomaly
ifcbbray.datetime=ifcbbray.datetime;
ifcbbray.dayOfYear=day(ifcbbray.datetime, 'dayofyear');
aod_climatology=climatology(ifcbbray,'AOD_500nm');




% Goal: Plot with taxa on y, time lag on the x and 

% Extract the taxa columns and AOD_500nm column
taxa_data = table2array(ifcbbray(:, 3:53));
aod_500nm =ifcbbray.AOD_500nm;

% Initialize variables to store coefficients and lags
cc = zeros(size(taxa_data, 1)*2-1, size(taxa_data, 2));
lags = zeros(size(taxa_data, 1)*2-1, size(taxa_data, 2));  % Corrected dimensions

taxa = ifcbbray.Properties.VariableNames(3:53);
clf;
for i = 1:size(taxa_data, 2)
    % Compute climatology and anomaly
    eval([taxa{i}, '_climatology = climatology(ifcbbray, taxa{i});']);
    eval(['ifcbbray = join(ifcbbray, ', taxa{i}, '_climatology, ''Keys'', ''dayOfYear'');']);
    eval(['ifcbbray.', taxa{i}, '_anomaly = ifcbbray.', taxa{i}, ' - ifcbbray.Fun_', taxa{i}, ';']);

    % Cross-correlate
    eval(['[cc(:, i), temp_lags] = xcorr(fillmissing(aod_500nm, ''linear''), ifcbbray.', taxa{i}, '_anomaly, ''coeff'');']);
    lags(:, i) = temp_lags;  % Correct the assignment to maintain consistency in dimensions

    if plotting == 1
        figure(1);
        plot(lags(:, i), cc(:, i), 'LineStyle', '-', 'LineWidth', 2);
        hold on;
        title("Cross-correlation between AOD 500nm &" + newline + "Phytoplankton Groups");
        xlabel('Lags (Days)');
        ylabel('Normalized Cross-correlation');
        xlim([-50, 50]);
    end
end

if plotting == 1
    legend(taxa_labels, 'Location', 'best'); % Add legend with dynamic taxa names
    set(gcf, 'Position', [-50, 1000, 1600, 800]);  % Adjust the figure position and size
end

%%
load('sc_pm2.5_daily_all.mat')
sc_pm2_5_daily=sc_pm2_5_daily(sc_pm2_5_daily.LocalSiteName=="San Lorenzo Valley Middle School",:);
% Initialize parameters
plotting = 0;

sc_pm2_5_daily_join=sc_pm2_5_daily(:,["date","pm2_5"]);
% Process data for the year 2020
% sc_pm2_5_daily = sc_pm2_5_daily(sc_pm2_5_daily.datetime.Year == 2020, :);

% Merge the sc_pm2_5_daily data with ifcbbray based on datetime
ifcb_pm25 = outerjoin(ifcbbray, sc_pm2_5_daily_join, 'Keys', 'date', 'MergeKeys', true, 'Type', 'left');


% Assuming ifcb_pm25 is already loaded into the workspace

% Extract the taxa columns and AOD_500nm column
taxa_data = table2array(ifcb_pm25(:, 3:53));
pm25 = ifcb_pm25.pm2_5;

% Initialize variables to store coefficients and lags
cc = zeros(size(taxa_data, 1)*2-1, size(taxa_data, 2));
lags = zeros(size(taxa_data, 1)*2-1, size(taxa_data, 2));

taxa = ifcb_pm25.Properties.VariableNames(3:53);
plotting = 1; % Assuming you want to plot

% Loop through each taxon
for i = 1:size(taxa_data, 2)
    % Cross-correlation
    eval(['[cc(:, i), temp_lags] = xcorr(fillmissing(pm25, ''linear''), ifcb_pm25.', taxa{i}, '_anomaly, ''coeff'');']);
    lags(:, i) = temp_lags;

    if plotting == 1
        figure(1);
        plot(lags(:, i), cc(:, i), 'LineStyle', '-', 'LineWidth', 2);
        hold on;
        title("Cross-correlation between AOD 500nm & Phytoplankton Groups");
        xlabel('Lags (Days)');
        ylabel('Normalized Cross-correlation');
        xlim([-50, 50]);
    end
end

if plotting == 1
    legend(taxa, 'Location', 'best');
    set(gcf, 'Position', [-50, 1000, 1600, 800]);
end

%% Find taxa with significant CC values at positive lags
threshold = 0.35;

significant_taxa = {};
plot_data = {};

for i = 1:size(taxa_data, 2)
    % Compute the cross-correlation and lags if not precomputed
    [cc(:, i), temp_lags] = xcorr(ifcb_pm25.([taxa{i}, '_anomaly']), fillmissing(pm25, 'linear'), 'coeff');
    lags(:, i) = temp_lags;

    % Compute first derivative of cc values
    cc_derivative = diff(cc(:, i));
    lags_midpoints = (lags(1:end-1, i) + lags(2:end, i)) / 2;
    
    % Find indices where derivative changes from positive to negative (relative maxima)
    maxima_indices = find(cc_derivative(1:end-1) > 0 & cc_derivative(2:end) < 0) + 1;
    maxima_lags = lags(maxima_indices, i);
    
    % Filter maxima within -50 to 0 lag range
    relevant_maxima_indices = maxima_indices(maxima_lags <= 50 & maxima_lags >= -1);

    % Check if the conditions are met: significant local maximum
    for j = 1:length(relevant_maxima_indices)
        if cc(relevant_maxima_indices(j), i) > threshold
            significant_taxa = [significant_taxa, taxa{i}];
            plot_data{end+1} = struct('cc', cc(:, i), 'lags', lags(:, i), 'label', taxa{i});
            break;
        end
    end
end

% Display the list of significant taxa
disp('Significant taxa with relative maximum cc between 50 and 0:');
disp(significant_taxa);

%% Plotting the results for significant taxa

figure;

subplot(2,1,1); % Lower plot for time series data
hold on;

start_date = datetime(2020, 7, 1);
end_date = datetime(2020, 11, 1);

yyaxis left;
colors_95 = lines(length(significant_taxa));

for k = 1:length(significant_taxa)
    taxa_data = ifcb_pm25.(significant_taxa{k});
    plot(ifcb_pm25.date, taxa_data, 'DisplayName', significant_taxa{k}, 'LineWidth', 2, 'LineStyle', '-', 'Color', colors_95(k, :));
    ylabel('Abundance');
end

yyaxis right;
pm25_plot = plot(ifcb_pm25.date, fillmissing(pm25, 'linear'), 'DisplayName', 'AOD 500 nm', 'LineStyle', '-', 'LineWidth', 2, 'Color', 'k');
ylabel('PM 2.5');

title('Time Series Data for Significant Taxa');
legend('show', 'Location', 'best');
datetick('x', 'yyyy-mm-dd', 'keeplimits', 'keepticks');
xlim([start_date end_date]);

shade_start = datetime(2020,8,16);
shade_end = datetime(2020,9,22);
hold on;
yyaxis right;
fill([shade_start, shade_start, shade_end, shade_end], [0, max(ylim) * 1.1, max(ylim) * 1.1, 0], [0, 0.8, 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

legend(flipud(findobj(gca, 'Type', 'Line')), 'Location','eastoutside');

ax = gca;
ax.XTick = start_date:calmonths(1):end_date;
datetick('x', 'mmm yyyy', 'keeplimits', 'keepticks');

subplot(2,1,2); % Upper plot for cross-correlations
hold on;
for k = 1:length(plot_data)
    plot(plot_data{k}.lags, plot_data{k}.cc, 'DisplayName', plot_data{k}.label, 'LineWidth', 2);
end
xlim([-50 50]);
ylim([0 0.8]);
title('Cross-Correlation between AOD 500nm and Significant Phytoplankton Groups');
xlabel('Lags (Days after AOD anomaly)');
ylabel('Normalized Cross-correlation');
legend('show', 'Location','eastoutside');

set(gcf, 'Position', [0, 0, 1600, 900]);

saving = 0;
if saving == 1
    saveas(gcf, "xcorr_time_series_plot.png");
end



%% Create a heatmap

xmax=30;
xmin=xmax*-1;
taxa=ifcbbray.Properties.VariableNames(3:53);
taxa = strrep(taxa, '_', ' ');
lag=lags(:,1);

%Plotting 
XLabels = lag;
% Convert each number in the array into a string
CustomXLabels = string(XLabels);
% Replace all but the fifth elements by spaces
CustomXLabels(mod(XLabels,xmax/10) ~= 0) = " ";
% Set the 'XDisplayLabels' property of the heatmap 
% object 'h' to the custom x-axis tick labels

figure()
hBar=heatmap(lag,taxa,cc');
% heatmap(ifcbbray.Properties.VariableNames(1:51), lag,cc')
xlim([-xmin,xmax])
% squeeze(correlationMatrix(:,:,numLags+1))
colorbar; % Add colorbar to show strength of coefficients
title('Cross-correlation between PM 2.5 Anomaly and IFCB Phytoplankton Abundance Anomaly');
xlabel('Lag (days after AOD anomaly)');
ylabel('Taxa')

colormap('jet')
caxis([0 0.8])
xlim([0 28]);

hBar.XDisplayLabels = CustomXLabels;
hBar.GridVisible = 'off';
s = struct(hBar);
s.XAxis.TickLabelRotation = 0;  % vertical
s.XAxis.FontSize=14;
s.YAxis.FontSize=14;

set(gcf,'Position',[0 0 900 1000])

saving=1;
if saving==1
    % Save figure as a PNG file in a specific directory
    directory = 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\xcorr'; % Specify your directory
    filename = ['xcorr_heatmap_',num2str(xmax),'_lags_anomaly_pm2.5.png']; % Specify the file name
    fullpath = fullfile(directory, filename); % Create full path
    saveas(gcf, fullpath); % Save the figure
end


%% Heatmap of Xcorr of variables 

% Create a heatmap

xmax=30;
xmin=xmax*-1;
taxa=ifcbbray.Properties.VariableNames(3:53);
taxa = strrep(taxa, '_', ' ');
lag=lags(:,1);

%Plotting 
XLabels = lag;
% Convert each number in the array into a string
CustomXLabels = string(XLabels);
% Replace all but the fifth elements by spaces
CustomXLabels(mod(XLabels,xmax/10) ~= 0) = " ";
% Set the 'XDisplayLabels' property of the heatmap 
% object 'h' to the custom x-axis tick labels

figure()
hBar=heatmap(lag,taxa,cc');
% heatmap(ifcbbray.Properties.VariableNames(1:51), lag,cc')
xlim([-xmin,xmax])
% squeeze(correlationMatrix(:,:,numLags+1))
colorbar; % Add colorbar to show strength of coefficients
title('Cross-correlation between PM 2.5 Anomaly and IFCB Phytoplankton Abundance Anomaly');
xlabel('Lag (days after PM 2.5 anomaly)');
ylabel('Taxa')

colormap('jet')
caxis([0 0.8])
xlim([0 28]);

hBar.XDisplayLabels = CustomXLabels;
hBar.GridVisible = 'off';
s = struct(hBar);
s.XAxis.TickLabelRotation = 0;  % vertical
s.XAxis.FontSize=32;
s.YAxis.FontSize=28;

set(gcf,'Position',[0 0 1400 1000])

saving=1;
if saving==1
    % Save figure as a PNG file in a specific directory
    directory = 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\xcorr'; % Specify your directory
    filename = ['xcorr_heatmap_',num2str(xmax),'_lags_anomaly_pm2.5.png']; % Specify the file name
    fullpath = fullfile(directory, filename); % Create full path
    saveas(gcf, fullpath); % Save the figure
end

%% Just plot PM 2.5 time series
clf
gcf
figure(1)
plot(sc_pm2_5_daily.date(sc_pm2_5_daily.date.Year==2020),fillmissing(sc_pm2_5_daily.pm2_5(sc_pm2_5_daily.date.Year==2020),'linear'))


%% Anomaly

start_date=datenum(datetime(2020,7,1));
end_date=datenum(datetime(2020,11,1));

shade_start = datenum(datetime(2020,8,1));
shade_end = datenum(datetime(2020,10,1));

figure(1)
shade_anomaly(datenum(sc_pm2_5_daily.date(sc_pm2_5_daily.date.Year==2020)),sc_pm2_5_daily.pm2_5(sc_pm2_5_daily.date.Year==2020))


fill([shade_start, shade_start, shade_end, shade_end], ...
     [min(ylim) * 1.1, max(ylim) * 1.1, max(ylim) * 1.1, min(ylim) * 1.1], ...
     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4); % Adjust 'FaceAlpha' for transparency

% Set x-axis ticks and labels
datetick('x')
% xticks(date_range);
% xticklabels(datestr(date_range, 'mmm, yyyy'));
xtickangle(30);  % Angle the dates for better visibility
xlim([datenum(start_date),datenum(end_date)])


%% Regresss pm 2.5 on taxa
colorz = [
    0.2, 0.6, 0.8;  % Light blue
    0.9, 0.5, 0.1;  % Orange
    0.4, 0.7, 0.2;  % Light green
    0.8, 0.2, 0.4;  % Pink
    0.2, 0.4, 0.8;  % Dark blue
    0.7, 0.2, 0.9   % Purple
];

% Calculate cross-correlation
% Define the range of negative lags
negative_lag_range = -50:0;
% Initialize arrays to store the optimal lags
optimal_lags = zeros(1, length(significant_taxa));

for i = 1:length(significant_taxa)
    % Get the anomaly data for the current taxon
    taxon_anomaly = ifcb_pm25.([significant_taxa{i}, '_anomaly']);
    % Compute cross-correlation for negative lags
    [cc, lags] = xcorr(fillmissing(pm25, 'linear'), taxon_anomaly, 'coeff');
    negative_indices = lags <= 0;
    cc_negative = cc(negative_indices);
    lags_negative = lags(negative_indices);
    % Find the maximum cross-correlation and corresponding lag
    [~, max_index] = max(cc_negative);
    optimal_lags(i) = lags_negative(max_index);
end


% Regression and plotting
clf;
figure(1);

for k = 1:length(significant_taxa)
    subplot(length(significant_taxa)/2, 2, k);
    taxa_data = ifcb_pm25.(significant_taxa{k});
    pm25_data = ifcb_pm25.pm2_5;
    dates = ifcb_pm25.datetime;

    % Shift AOD data by optimal lag
    optimal_lag = optimal_lags(k);
    pm25_data_shifted = circshift(pm25_data, abs(optimal_lag));

    % Filter out NaNs and zeros
    validIndices = ~isnan(pm25_data_shifted) & ~isnan(taxa_data) & pm25_data_shifted > 0 & taxa_data > 0;
    pm25_data_shifted = pm25_data_shifted(validIndices);
    taxa_data = taxa_data(validIndices);
    valid_dates = dates(validIndices);

    % Convert dates to datenum for easier comparison
    date_start = datetime('17-Aug-2020');
    date_end = datetime('22-Sep-2020');

    if length(pm25_data_shifted) > 1 && length(taxa_data) > 1
        % Log transform the data
        logAOD = (pm25_data_shifted);
        logAbundance = (taxa_data);

        % Scatter plot
        scatter(logAOD, logAbundance, 30, 'k', 'filled', 'DisplayName', significant_taxa{k}, 'MarkerFaceAlpha', 0.5);
        hold on;

        % Linear regression for all data
        p_all = polyfit(logAOD, logAbundance, 1);
        yfit_all = polyval(p_all, logAOD);
        % plot(logAOD, yfit_all, 'k-', 'LineWidth', 2);

        % Calculate R² for all data
        SS_tot_all = sum((logAbundance - mean(logAbundance)).^2);
        SS_res_all = sum((logAbundance - yfit_all).^2);
        R2_all = 1 - (SS_res_all / SS_tot_all);

        % Display R² for all data on the plot
        % text(max(logAOD) - max(logAOD) / 10, max(logAbundance) - max(logAbundance) / 10, sprintf('R² = %.2f', R2_all), 'FontSize', 12, 'Color', colorz(k, :));

        % Highlight points between 08/2020 and 09/2020
        highlight_indices = valid_dates >= date_start & valid_dates <= date_end;
        highlighted_dates = valid_dates(highlight_indices);
        scatter(logAOD(highlight_indices), logAbundance(highlight_indices), 40, colorz(k, :), 'filled', 'Marker', 'diamond', 'MarkerFaceAlpha', 0.5);

        % Add text labels for the highlighted dates
        
        % for j = 1:length(highlighted_dates)
        %     idx=find(valid_dates==highlighted_dates(j));
        %     highlight_date = highlighted_dates(j);
        %     highlight_logAOD = logAOD(idx);
        %     highlight_logAbundance = logAbundance(idx);
        %     text(highlight_logAOD, highlight_logAbundance, datestr(highlight_date, 'mm/dd/yyyy'), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', colorz(k, :));
        % end

    else
        % Not enough data points for regression
        scatter(logAOD, logAbundance, 10, 'DisplayName', significant_taxa{k});
        text(0.1, 0.9, 'Insufficient data', 'Units', 'normalized', 'FontSize', 12);
    end

    % ylabel(['log Abundance [Cells/L]']);
    % xlabel('log AOD 500nm (shifted)');
    ylabel(['Abundance [Cells/L]']);
    xlabel('PM 2.5 (shifted)');
    title(significant_taxa{k})
    hold off;
    sgtitle(['Regression between PM 2.5 (shifted to lag with max Xcorr) and Phytoplankton Abundance',newline,'Colored Points between ',datestr(date_start),'-',datestr(date_end)])
end

set(gcf, 'Position', [0,1000, 1600, 800]);  % Adjust the figure position and size

% set(gcf, 'Position', [-50, 1000, 1600, 800]);  % Adjust the figure position and size

saving=1;
if saving==1
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\ifcb_vs_aeronet\taxon_pm25_shifted_regression.png"]);
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\ifcb_vs_aeronet\taxon_pm25_shifted_regression.pdf"]);

end

% saving=1;
% if saving==1
%     saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\ifcb_vs_aeronet\taxon_aod_shifted_regression\taxon_aod_shifted_regression_log.png"]);
%     saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\ifcb_vs_aeronet\taxon_aod_shifted_regression\taxon_aod_shifted_regression_log.pdf"]);
% 
% end