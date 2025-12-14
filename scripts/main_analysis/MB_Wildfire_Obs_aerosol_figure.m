%% Main Analysis for Capone et al "Extreme wildfire conditions shift coastal phytoplankton community structure in California"



%% Add paths



%Project

addpath(genpath('data\'));
addpath(genpath('functions\'));





%% Load data



%I. Wildfire Aerosols

%Aeronet

load('data\aeronet\monterey_lvl_1.5.mat')



% Load PM2.5 dataset

load('data\pm2_5\sc_pm2.5_daily_all.mat')

% Santa Cruz

sc_pm2_5_daily_sc=sc_pm2_5_daily(sc_pm2_5_daily.LocalSiteName=="Santa Cruz",:);

% % SLV Middle

sc_pm2_5_daily_slv=sc_pm2_5_daily(sc_pm2_5_daily.LocalSiteName=="San Lorenzo Valley Middle School",:);



% MERRA-2 averaged for the monterey Bay

load('data\merra-2\merra_MB_bc_avg_daily_plt.mat')



% II. Chlorophyll



%satellite

load('data\satellite\new_climatologies_anomalies\MB_satellite_climatologies_struct.mat')

load('data\satellite\satellite_chl_MB_5day_climatology_data.mat')





%HABS

load('data\shore_stations\habs_scwharf.mat')

HABs_SantaCruzWharf.weekofyear=week(HABs_SantaCruzWharf.datetime,'weekofyear');

HABs_SantaCruzWharf_clean = HABs_SantaCruzWharf(HABs_SantaCruzWharf.Avg_Chloro < 50, :);



%Shore station

load('data\shore_stations\mlml_station_daily.mat')





% III. IFCB



% Taxonomy

load('data\ifcb_processed\ifcb_ucsc_all_mean.mat')

load('data\ifcb_processed\ifcb_bray.mat')





% IV. Physical Data

load('data\joined_physical_drivers\physical_forcings_biology_table_alltime.mat')











%% Set plotting specs

% Set publication settings

ftsz = 10;  % Font size for publication

ftname = 'Helvetica';  % Font name

linewdt = 2;  % Line width for plots





% Dates



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



shade_start= start_date;

shade_end= end_date;





shade_start_day = start_date_day;

shade_end_day = end_date_day;



shade_start_5day = start_date_5day;

shade_end_5day = end_date_5day;


%% I. Figure 1. Wildfire Aerosols Analysis

%Select fire  years to exlcude

fire_years=[2008, 2016,2020,2021];



legend_on = 0;

num_subplots = 3;

saving=0;

colorz = [

    0.75, 0.1, 0;        % Darker shade of Red

    0.3333, 0.6588, 0.4078;  % Green

    0.1882, 0.2784, 0.7;        % Blue

    0.7882, 0.4588, 0.9216];   % Pastel Purple





ftsz=8;

linewdt=1;



% Plot

close all;

figure()



% PM 2.5 Plot

figure()

t = tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile



yyaxis left;

plot(sc_pm2_5_daily_slv.dayOfYear(sc_pm2_5_daily_slv.datetime.Year == 2020), ...
    sc_pm2_5_daily_slv.pm2_5(sc_pm2_5_daily_slv.datetime.Year == 2020), ...
    'Color', colorz(3,:), 'LineStyle', '-', 'LineWidth', linewdt);

hold on;

%Add shading

plot(sc_pm2_5_daily_sc.dayOfYear(sc_pm2_5_daily_sc.datetime.Year == 2020), ...
    sc_pm2_5_daily_sc.pm2_5(sc_pm2_5_daily_sc.datetime.Year == 2020), ...
    'Color', colorz(2,:), 'LineStyle', '-', 'LineWidth', linewdt);

ylim([0 400])

ylabel(["PM 2.5 (\mu g m ^{-3})"], 'FontSize', ftsz, 'FontName', ftname);



yyaxis right;

plotClimatology5day(sc_pm2_5_daily_slv(sc_pm2_5_daily_slv.datetime.Year ~= 2020, :), 'pm2_5')

alpha(0.5)

ylim([0 30])

xlim([1 365])

hold on;

ylabel(["Average PM 2.5" + newline + "(\mu g m^{-3})"], 'FontSize', ftsz, 'FontName', ftname);



%Add shading

% xline(start_date_day, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-'); 

% fill([shade_start_day, shade_start_day, shade_end_day, shade_end_day], ...

%     [0, max(ylim) * 1.5, max(ylim) * 1.5 0], ...

%     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.2); 

% xline(end_date_day, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');



ax = gca;

ax.FontSize = ftsz-1;

ax.FontName = ftname;

ax.XTickLabelRotation=45;

set(gca, 'XTickLabel', []);





if legend_on == 1

    legend({'Climatology', '2020'}, 'FontSize', ftsz, 'FontName', ftname,'Location','northwest')

end



% Add "a)" label

text(-0.21, 1.1, 'c', 'Units', 'normalized', 'FontSize', ftsz+2, 'FontName', ftname,'FontWeight','bold');



% AOD Plot

nexttile

yyaxis left;

plot(Monterey_lvl_1_5.Day_of_Year(Monterey_lvl_1_5.datetime.Year == 2020), ...
    Monterey_lvl_1_5.AOD_500nm(Monterey_lvl_1_5.datetime.Year == 2020), ...
    'Color', colorz(1,:), 'LineStyle', '-', 'LineWidth', 1);

hold on



ylabel('AOD 500nm', 'FontSize', ftsz, 'FontName', ftname);

xlim([0 365])



yyaxis right;

plotClimatology5day(Monterey_lvl_1_5(~ismember(Monterey_lvl_1_5.datetime.Year, fire_years), :), 'AOD_500nm')

ylabel(["Average AOD" + newline + "500nm"], 'FontSize', ftsz, 'FontName', ftname);



alpha(0.1)

ylim([0 0.2])

xlim([1 365])

hold on;

%Add shading

% xline(start_date_day, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-'); 

% fill([shade_start_day, shade_start_day, shade_end_day, shade_end_day], ...

%     [0, max(ylim) * 1.5, max(ylim) * 1.5 0], ...

%     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.2); 

% xline(end_date_day, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');





ax = gca;

ax.FontSize = ftsz-1;

ax.FontName = ftname;

ax.XTickLabelRotation=45;

set(gca, 'XTickLabel', []);



if legend_on == 1

    legend({'Climatology (2003-present)', '2020'}, 'FontSize', ftsz, 'FontName', ftname)

end

ax = gca;

ax.FontSize = ftsz-1;

ax.FontName = ftname;



% Add "b)" label

text(-0.21, 1, 'd', 'Units', 'normalized', 'FontSize', ftsz+2, 'FontName', ftname,'FontWeight','bold');



% Black Carbon Mass Plot

nexttile



yyaxis left;

plot(merra_bc_plt.dayOfYear(merra_bc_plt.datetime.Year == 2020), ...
    merra_bc_plt.Median_BCCMASS(merra_bc_plt.datetime.Year == 2020), ...
    'Color', colorz(4,:), 'LineStyle', '-', 'LineWidth', linewdt);

ylabel(["Black Carbon" + newline + "Concentration (µg m^{-2})"], 'FontSize', ftsz, 'FontName', ftname);

hold on







yyaxis right;

plotClimatology5day(merra_bc_plt(~ismember(merra_bc_plt.datetime.Year, fire_years), :), 'Median_BCCMASS')

ylabel(["Average Black Carbon" + newline + "Concentration (µg m^{-2})"], 'FontSize', ftsz, 'FontName', ftname);

alpha(0.2)

ylim([3e-7 12e-7])

xlim([1 365])

hold on;

%Add shading

% xline(start_date_day, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-'); 

% fill([shade_start_day, shade_start_day, shade_end_day, shade_end_day], ...

%     [0, max(ylim) * 1.5, max(ylim) * 1.5 0], ...

%     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.2); 

% xline(end_date_day, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');



ax = gca;

ax.FontSize = ftsz-1;

ax.FontName = ftname;

ax.XTickLabelRotation=45;







if legend_on == 1

    legend({'Climatology (2003-present)', '2020'}, 'FontSize', ftsz, 'FontName', ftname)

end



% Add "c)" label

text(-0.21, 1, 'e', 'Units', 'normalized', 'FontSize', ftsz+2, 'FontName', ftname,'FontWeight','bold');



% Adjust figure size for publication (double-column width)

set(gcf, 'Units', 'inches', 'Position', [0, 0, 2.5, 3.6]);  % 5 inches width and 5 inches height for subplots

set(gcf, 'PaperPositionMode', 'auto');  % Ensure the figure fits the paper size



% Save the figure if required

saving = 1;

if saving == 1

    print(gcf, '-dpng', '-r1200', 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\final\figure_1_v0.png');

    print(gcf, '-dpdf', '-r1200', 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\final\figure_1_v0.pdf');

end

saving = 0;








