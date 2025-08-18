%% Chlorophyll-a Analysis

% Add paths
addpath(genpath('Q:\Dante\Wildfire_Obs'));
addpath(genpath('Q:\Dante\data\MB_Wildfire_Obs\processed_data'));
addpath('Q:\Dante\data\satellite_chl_CCE\shadow\chl\')

%Plot settings
ftsz=24;
linewdt=3;
ftname='Helvetica';

%% Load the data

%satellite
load('Q:\Dante\data\MB_Wildfire_Obs\processed_data\satellite\new_climatologies_anomalies\MB_satellite_climatologies_struct.mat')
load('satellite_chl_MB_5day_climatology_data.mat')


%HABS
load('Q:\Dante\data\MB_Wildfire_Obs\shore_stations\habs_scwharf.mat')

%Shore station
load('Q:\Dante\data\MB_Wildfire_Obs\shore_stations\mlml_station_daily.mat')





%% HABS
HABs_SantaCruzWharf.weekofyear=week(HABs_SantaCruzWharf.datetime,'weekofyear');

clf
yearss=unique(HABs_SantaCruzWharf.datetime.Year);
for i=1:length(unique(HABs_SantaCruzWharf.datetime.Year))
    plot(HABs_SantaCruzWharf.weekofyear(HABs_SantaCruzWharf.datetime.Year==yearss(i)),HABs_SantaCruzWharf.Chl2(HABs_SantaCruzWharf.datetime.Year==yearss(i)))
hold on
end
datetick('x','mm/dd', 'keepticks')

%% Shading parameters for plots to show wildfire period and important dates

lag_days = 6;
start_date = datetime(2020, 8, 16);
end_date = datetime(2020, 9, 22);
lag_date = datetime(2020, 8, 21 + lag_days);
august_date = datetime(2020, 9, 5);
czu_date=datetime(2020, 8, 21);
%By day
start_date_day = day(start_date, 'dayofyear');
end_date_day = day(end_date, 'dayofyear');
lag_date_day = day(lag_date, 'dayofyear');
august_date_day = day(august_date, 'dayofyear');
czu_date_day=day(czu_date,'dayofyear');
%By week
start_date_week = week(start_date, 'weekofyear');
end_date_week = week(end_date, 'weekofyear');
lag_date_week = week(lag_date, 'weekofyear');
august_date_week = week(august_date, 'weekofyear');
czu_date_week = week(czu_date, 'weekofyear');


% By 5-day window
% Each 5-day window starts at day 1, 6, 11, 16, ..., so we divide by 5 and round up
start_date_5day = ceil(start_date_day / 5);
end_date_5day = ceil(end_date_day / 5);
lag_date_5day = ceil(lag_date_day / 5);
august_date_5day = ceil(august_date_day / 5);
czu_date_5day = ceil(czu_date_day / 5);

%Years
start_year=2016;
end_year=2023;


%For shading
shade_start_day = start_date_day;
shade_end_day = end_date_day;

shade_start_5day = start_date_5day;
shade_end_5day = end_date_5day;

shade_start_week = start_date_week;
shade_end_week = end_date_week;

HABs_SantaCruzWharf_clean = HABs_SantaCruzWharf(HABs_SantaCruzWharf.Avg_Chloro < 50, :);

%% 
figure(1)
clf
weekly_climatology(HABs_SantaCruzWharf_clean, 'Avg_Chloro')
hold on
% Add shaded area to the plot
xline(start_date_week, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
fill([shade_start_week, shade_start_week, shade_end_week, shade_end_week], ...
    [min(ylim) * 4.5, max(ylim) * 1.5, max(ylim) * 1.5 min(ylim) * 4.5], ...
    [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.3); 
xline(end_date_week, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');


hold on

% Add lines for 8 day lag and August Fire Smoke
% xline(lag_date_week, 'Color', '#cf572b', 'LineWidth', linewdt, 'LineStyle', ':'); 
% xline(august_date_week, 'Color', 'k', 'LineWidth', linewdt, 'LineStyle', ':');

%
ylim([-5 18])
% Add text labels at 30-degree angles
text(lag_date_week, 15, '8-day lag', 'FontSize', ftsz-4, 'FontName', 'Helvetica', ...
    'Rotation', 0, 'HorizontalAlignment', 'left','Color', '#cf572b');
text(august_date_week, 12, 'August Fire Smoke', 'FontSize', ftsz-4, 'FontName', 'Helvetica', ...
    'Rotation', 0, 'HorizontalAlignment', 'left','Color', 'k');

hold on
ylabel('Average Chlorophyll mg m^{-3}', 'FontSize', ftsz, 'FontName', ftname)

% The rest of your plotting code remains the same
yyaxis right
plot(HABs_SantaCruzWharf_clean.weekofyear(HABs_SantaCruzWharf_clean.datetime.Year == 2020), ...
    HABs_SantaCruzWharf_clean.Avg_Chloro(HABs_SantaCruzWharf_clean.datetime.Year == 2020), ...
    'LineStyle', '-', 'LineWidth', linewdt, 'Marker', 'none','Color','#285c36');
hold on


xlim([1 51])
ylabel('Chlorophyll mg m^{-3}', 'FontSize', ftsz, 'FontName', ftname)
xlabel('')
title('HABS')
legend({'Climatology', '2020'})

ax = gca;  % Get current axes
ax.FontSize = 16;
ax.FontName = 'Helvetica';

set(gcf, 'Position', [0 100 1600 600])



%% MLML Shore station fluorometer

%Compare different averaging periods

close all;
%Daily
subplot(4,1,1)
plotClimatology(mlml_tt(mlml_tt.datetime.Year>2017,:),'fluorescence')
xlim([1 365])

%5day
subplot(4,1,2)
plotClimatology5day(mlml_tt(mlml_tt.datetime.Year>2017,:),'fluorescence')
xlim([1 365])
%15 Day
subplot(4,1,3)
plotClimatology15day(mlml_tt(mlml_tt.datetime.Year>2017,:),'fluorescence'L)
xlim([1 365])

subplot(4,1,4)
plotClimatologyMonthly(mlml_tt(mlml_tt.datetime.Year>2017,:),'fluorescence')



%%
close all;
figure()
plot(mlml_tt.datetime(mlml_tt.fluorescence_qc_agg==1),mlml_tt.fluorescence(mlml_tt.fluorescence_qc_agg==1))

%%
mlml_tt=mlml_tt(~(isnan(mlml_tt.dayofyear)),:);
mlml_climatology=climatology_se(mlml_tt,'fluorescence');
mlml_tt_join = join(mlml_tt, mlml_climatology, 'Keys', 'dayofyear');
mlml_tt_join.anomaly_fluorescence=mlml_tt_join.fluorescence-mlml_tt_join.Mean_fluorescence;


%%
clf;
close all;
plotClimatology(mlml_tt,'fluorescence')
hold on
yyaxis right
plot(mlml_tt.dayofyear(mlml_tt.datetime.Year==2020),mlml_tt.fluorescence(mlml_tt.datetime.Year==2020),'LineStyle','-',...
    'LineWidth', linewdt)
set(gcf,'Position',[0 100 1600 600])

%% Satellite

% For axis labels
month_labs = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};

%Compare Bay to Upwelling Shadow
clf; figure(1)

%Full Bay
subplot(2,1,1)
xvals = 1:length(satellite_chl_MB_5day_climatology_daily_data.climatology);
fill([xvals, fliplr(xvals)], [satellite_chl_MB_5day_climatology_daily_data.climatology' + satellite_chl_MB_5day_climatology_daily_data.error', fliplr(satellite_chl_MB_5day_climatology_daily_data.climatology' - satellite_chl_MB_5day_climatology_daily_data.error')], 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
hold on;
plot(xvals, satellite_chl_MB_5day_climatology_daily_data.climatology, '-','LineWidth', linewdt,'Color','#979da6');
xlabel('Month');
ylabel('Average Chlorophyll-a [mg/m^3]');
xlim([0 73])
set(gca, 'XTick', 1:73/12:73, 'XTickLabel', month_labs, 'XTickLabelRotation', 0);  % Use the same x-axis labels
grid on;

%2020 data
year_sel=2020;
year_ind=find(satellite_chl_MB_5day_climatology_daily_data.years==year_sel);
yearLabel = sprintf('%d', year_sel);  % Create a label for the current year
plot(1:size(satellite_chl_MB_5day_climatology_daily_data.dailyData, 1), satellite_chl_MB_5day_climatology_daily_data.dailyData(:, year_ind), '-', 'Color',...
    'r', 'DisplayName', yearLabel, 'LineWidth', linewdt);

ylabel('Full Bay Chlorophyll-a [mg/m^3]');
legend(["Climatology";"2020"],'AutoUpdate','Off')

% Add lines for 8 day lag and August Fire Smoke
% xline(lag_date_5day, 'Color', '#cf572b', 'LineWidth', linewdt, 'LineStyle', ':'); 
% xline(august_date_5day, 'Color', 'k', 'LineWidth', linewdt, 'LineStyle', ':');




% Add shaded area to the plot
xline(start_date_5day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
fill([shade_start_5day, shade_start_5day, shade_end_5day, shade_end_5day], ...
    [min(ylim) * 4.5, max(ylim) * 1.5, max(ylim) * 1.5 min(ylim) * 4.5], ...
    [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.3); 
xline(end_date_5day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');

title('Full Bay');


grid on;
ax = gca;  % Get current axes
ax.FontSize = 16;
ax.FontName = 'Helvetica';


%Shadow
subplot(2,1,2)
xvals = 1:length(satellite_chl_MB_5day_climatology_08_23_shadow.climatology);
fill([xvals, fliplr(xvals)], [satellite_chl_MB_5day_climatology_08_23_shadow.climatology' + satellite_chl_MB_5day_climatology_08_23_shadow.error', fliplr(satellite_chl_MB_5day_climatology_08_23_shadow.climatology' - satellite_chl_MB_5day_climatology_08_23_shadow.error')], 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
hold on;
plot(xvals, satellite_chl_MB_5day_climatology_08_23_shadow.climatology, '-','LineWidth', linewdt,'Color','#979da6');
xlim([0 73])
ylim([0 2])
set(gca, 'XTick', 1:73/12:73, 'XTickLabel', month_labs, 'XTickLabelRotation', 0);  % Use the same x-axis labels
grid on;
hold on
year_ind=find(satellite_chl_MB_5day_climatology_08_23_shadow.years==year_sel);
yearLabel = sprintf('%d', year_sel);  % Create a label for the current year
plot(1:size(satellite_chl_MB_5day_climatology_08_23_shadow.dailyData, 1), satellite_chl_MB_5day_climatology_08_23_shadow.dailyData(:, year_ind), '-', 'Color',...
    '#66b398', 'DisplayName', yearLabel, 'LineWidth', linewdt);
ylabel('Upwelling Shadow Chlorophyll-a [mg/m^3]');
legend(["Climatology";"2020"],'AutoUpdate','Off')


% Add lines for 8 day lag and August Fire Smoke
% xline(lag_date_5day, 'Color', '#cf572b', 'LineWidth', linewdt, 'LineStyle', ':'); 
% xline(august_date_5day, 'Color', 'k', 'LineWidth', linewdt, 'LineStyle', ':');


% Add shaded area to the plot
xline(start_date_5day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
fill([shade_start_5day, shade_start_5day, shade_end_5day, shade_end_5day], ...
    [min(ylim) * 4.5, max(ylim) * 1.5, max(ylim) * 1.5 min(ylim) * 4.5], ...
    [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.3); 
xline(end_date_5day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');


% Update legend to only include selected years
% plot(xlim, [0 0], 'k-','LineWidth', linewdt,'Color', [0,0,0, 0.5],'LineStyle','-');  % Adding a zero line
% legend(legendInfo(1:index-1), 'Location', 'northwest');  % Adjust legend to show only included years
title('Upwelling Shadow');
% xlabel('Month');
ylabel('Chlorophyll-a [mg/m^3]');
xlim([0 73])
ylim([0 1.5])
set(gca, 'XTick', 1:73/12:73, 'XTickLabel', month_labs, 'XTickLabelRotation', 0);  % Use the same x-axis labels
% set(gca, 'XTick', 1:30:366, 'XTickLabel', month_labs, 'XTickLabelRotation', 45);  % Use the same x-axis labels

grid on;
ax = gca;  % Get current axes
ax.FontSize = 16;
ax.FontName = 'Helvetica';

saving=1;
if saving==1
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\satellite_chl\fullbay_vs_shadow.png"]);
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\satellite_chl\fullbay_vs_shadow.pdf"]);
end

%% 2020 satellite


% Read the normal 2020 data from 'chl5d_EOF2shadow.csv'
% Make sure the file has columns like SYear, SDay, Mean for each time point
filename = 'Q:\Dante\data\satellite_chl_CCE\shadow\chl\chl5d_EOF2shadow.csv';  % Path to the CSV file
data_2020 = readtable(filename);

% Filter data for the year 2020
data_2020 = data_2020(data_2020.SYear == 2020, :);

% Map the day of the year (SDay) to 5-day periods (1-73)
dayOfYear_2020 = data_2020.SDay;
periodIndex_2020 = ceil(dayOfYear_2020 / 5);  % Convert day of the year to 5-day periods

close all; clf
% Plot 2020 normal data on the right y-axis
yyaxis right;

% Plot the 2020 data
plot(periodIndex_2020, data_2020.Mean, '-', 'Color', '#66b398', 'DisplayName', '2020', 'LineWidth', linewdt);
hold on
% Customize the right y-axis
ylabel('Chlorophyll-a (2020) [mg/m^3]');
legend(["Climatology", "2020"], 'AutoUpdate', 'Off');


% Add lines for 8 day lag and August Fire Smoke
xline(lag_date_5day, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
xline(august_date_5day, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':');




% Add shaded area to the plot
xline(start_date_5day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
fill([shade_start_5day, shade_start_5day, shade_end_5day, shade_end_5day], ...
    [min(ylim) * 4.5, max(ylim) * 1.5, max(ylim) * 1.5 min(ylim) * 4.5], ...
    [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.3); 
xline(end_date_5day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');

% Final plot customization
xlim([0 73]);
ylim([0 25])

set(gca, 'XTick', 1:73/12:73, 'XTickLabel', month_labs, 'XTickLabelRotation', 0);
grid on;

% Adjust font and other aesthetics
ax = gca;
ax.FontSize = 16;
ax.FontName = 'Times';

title('Satellite Chlorophyll-a Data for EOF2 (2020 vs Climatology)');

dual_monitor=0;
if dual_monitor == 1
    set(gcf, 'Position', [0 1400 1400 600])
else
    set(gcf, 'Position', [0 100 1400 1000])
end


%% Plot Satellite, CalHABMAP, and Fluorometer with Helvetica Font Size 10 and Adjusted Figure Dimensions

close all;

% Set publication settings
ftsz = 10;  % Font size for publication
ftname = 'Helvetica';  % Font name
linewdt = 1.5;  % Line width for plotting

% Fluorometer
subplot(3,1,1);
weekly_climatology(mlml_tt(mlml_tt.datetime.Year > 2017, :), 'fluorescence', 2);

% ylabel('Average Chlorophyll mg m^{-3}', 'FontSize', ftsz, 'FontName', ftname);
ylim([0 20]);
hold on;
mlml_tt_plt = retime(mlml_tt, 'weekly', 'mean');
mlml_tt_plt.weekofyear = week(mlml_tt_plt.datetime, 'weekofyear');
plot(mlml_tt_plt.weekofyear(mlml_tt_plt.datetime.Year == 2020), mlml_tt_plt.fluorescence(mlml_tt_plt.datetime.Year == 2020), ...
    'LineStyle', '-', 'LineWidth', linewdt, 'Color', '#3e8a52');
% ylabel('Chlorophyll mg m^{-3}', 'FontSize', ftsz, 'FontName', ftname);
title('MLML Fluorometer', 'FontSize', ftsz, 'FontName', ftname);
xline(start_date_week, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
fill([shade_start_week, shade_start_week, shade_end_week, shade_end_week], ...
    [min(ylim) * 4.5, max(ylim) * 1.5, max(ylim) * 1.5, min(ylim) * 4.5], ...
    [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
xline(end_date_week, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
xline(lag_date_week, 'Color', '#cf572b', 'LineWidth', linewdt, 'LineStyle', ':');
xline(czu_date_week, 'Color', 'k', 'LineWidth', linewdt, 'LineStyle', ':');
ylim([0 20]);
ax = gca;
ax.FontSize = ftsz;
ax.FontName = ftname;

% Shore Station Water Sample
subplot(3,1,2);
weekly_climatology(HABs_SantaCruzWharf_clean, 'Avg_Chloro', 2);
hold on
% ylabel('Average Chlorophyll mg m^{-3}', 'FontSize', ftsz, 'FontName', ftname);
ylim([-5 30]);
xline(start_date_week, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
fill([shade_start_week, shade_start_week, shade_end_week, shade_end_week], ...
    [min(ylim) * 4.5, max(ylim) * 1.5, max(ylim) * 1.5, min(ylim) * 4.5], ...
    [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
xline(end_date_week, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
xline(lag_date_week, 'Color', '#cf572b', 'LineWidth', linewdt, 'LineStyle', ':');
xline(czu_date_week, 'Color', 'k', 'LineWidth', linewdt, 'LineStyle', ':');
plot(HABs_SantaCruzWharf_clean.weekofyear(HABs_SantaCruzWharf_clean.datetime.Year == 2020), ...
    HABs_SantaCruzWharf_clean.Avg_Chloro(HABs_SantaCruzWharf_clean.datetime.Year == 2020), ...
    'LineStyle', '-', 'LineWidth', linewdt, 'Color', '#50e678');
xlim([1 52]);
% ylabel('Chlorophyll mg m^{-3}', 'FontSize', ftsz, 'FontName', ftname);
title('Shore Station Water Sample', 'FontSize', ftsz, 'FontName', ftname);
ax = gca;
ax.FontSize = ftsz;
ax.FontName = ftname;

% Satellite
subplot(3,1,3);
month_labs = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
xvals = 1:73;
fill([xvals, fliplr(xvals)], ...
     [MB_satellite_climatologies_struct.EOF2shadow.d5.Mean' + MB_satellite_climatologies_struct.EOF2shadow.d5.CI95', ...
      fliplr(MB_satellite_climatologies_struct.EOF2shadow.d5.Mean' - MB_satellite_climatologies_struct.EOF2shadow.d5.CI95')], ...
      'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
hold on;
plot(xvals, MB_satellite_climatologies_struct.EOF2shadow.d5.Mean, '-', 'LineWidth', linewdt, 'Color', '#979da6');
hold on
% ylabel('Average Chlorophyll mg m^{-3}', 'FontSize', ftsz, 'FontName', ftname);
xlim([0 73]);
ylim([0 25]);
set(gca, 'XTick', 1:73/12:73, 'XTickLabel', month_labs, 'XTickLabelRotation', 45);
grid on;
plot(periodIndex_2020, data_2020.Mean, '-', 'Color', '#66b398', 'LineWidth', linewdt);
xline(lag_date_5day, 'Color', '#cf572b', 'LineWidth', linewdt, 'LineStyle', ':');
xline(czu_date_5day, 'Color', 'k', 'LineWidth', linewdt, 'LineStyle', ':');
xline(start_date_5day, 'Color', 'k', 'LineWidth', linewdt, 'LineStyle', '-');
fill([shade_start_5day, shade_start_5day, shade_end_5day, shade_end_5day], ...
    [min(ylim) * 4.5, max(ylim) * 1.5, max(ylim) * 1.5, min(ylim) * 4.5], ...
    [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
xline(end_date_5day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
title('Satellite', 'FontSize', ftsz, 'FontName', ftname);

ax = gca;
ax.FontSize = ftsz;
ax.FontName = ftname;


% Create an invisible axis for the shared label (for yyaxis left)
han_left = axes(gcf, 'visible', 'off');
han_left.YLabel.Visible = 'on';

% Move the left y-axis label further to the left by adjusting 'Position' property
ylabel(han_left, 'Chlorophyll-a mg m^{-3}', 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold');

% Create an invisible axis for the shared label (for yyaxis left)
han_bottom = axes(gcf, 'visible', 'off');
han_bottom.XLabel.Visible = 'on';

% Move the bottom y-axis label further to the bottom by adjusting 'Position' property
xlabel(han_bottom, 'Month', 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold');
set(get(han_bottom, 'XLabel'), 'Position', [0.5, -0.09, 0], 'VerticalAlignment', 'bottom');


% Set figure size for double-column width and aspect ratio (for L&O)
set(gcf, 'Units', 'inches', 'Position', [0, 0, 4, 6]);  % 7.16 inches width and height to maintain aspect ratio
set(gcf, 'PaperPositionMode', 'auto');  % Ensure the figure fits the paper size

% Save figure with the required font and size for publication
saving = 0;
if saving == 1
    print(gcf, '-dpng', '-r300', 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\final\figure_2_v0.png');
    print(gcf, '-dpdf', '-r300', 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\final\figure_2_v0.pdf');
end
saving = 0;
