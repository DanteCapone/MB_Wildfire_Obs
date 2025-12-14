%% Main Analysis for Capone et al "Extreme wildfire conditions shift coastal phytoplankton community structure in California"



%% Add paths



%Project

addpath(genpath('data\'));





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


%% III. Figure 3. IFCB Morphometric Analysis 

% Calculate necessary Climatologies



%Load in the compiled morphology table

load('data\ifcb_processed\ifcb_nass_morpho_all_2016_2023.mat')



ifcb_nass_morpho_nass_climatology=climatology_se(ifcb_nass_morpho_2016_2023,'MeanSlopeAbundance','week');

ifcb_biov_morpho_nass_climatology=climatology_se(ifcb_nass_morpho_2016_2023,'MeanBiomassSum','week');



%Add week for joining

ifcb_nass_morpho_2016_2023.week=week(ifcb_nass_morpho_2016_2023.datetime,'weekofyear');



%Join

ifcb_nass_morpho_2016_2023_join=join(ifcb_nass_morpho_2016_2023,ifcb_biov_morpho_nass_climatology,'Keys',{'week'});

ifcb_nass_morpho_2016_2023_join.anomaly_biovolumesum=ifcb_nass_morpho_2016_2023_join.MeanBiomassSum-ifcb_nass_morpho_2016_2023_join.MeanBiomassSum_climatology;



ifcb_nass_morpho_2016_2023_join=join(ifcb_nass_morpho_2016_2023_join,ifcb_nass_morpho_nass_climatology,'Keys',{'week'});





%% Compute Size Class Proportions

% Add new columns to the DataFrame

ifcb_nass_morpho_2016_2023.SumBinBiovolume1_2 = ifcb_nass_morpho_2016_2023.BinBiovolume1 + ifcb_nass_morpho_2016_2023.BinBiovolume2;

ifcb_nass_morpho_2016_2023.SumBinBiovolume3_4 = ifcb_nass_morpho_2016_2023.BinBiovolume3 + ifcb_nass_morpho_2016_2023.BinBiovolume4;

ifcb_nass_morpho_2016_2023.SumBinBiovolume5_6 = ifcb_nass_morpho_2016_2023.BinBiovolume5 + ifcb_nass_morpho_2016_2023.BinBiovolume6;



% Calculate total sum for each row

ifcb_nass_morpho_2016_2023.TotalSum = ifcb_nass_morpho_2016_2023.SumBinBiovolume1_2 + ...
                                      ifcb_nass_morpho_2016_2023.SumBinBiovolume3_4 + ...
                                      ifcb_nass_morpho_2016_2023.SumBinBiovolume5_6;



% Calculate proportions

ifcb_nass_morpho_2016_2023.PropBinBiovolume1_2 = ifcb_nass_morpho_2016_2023.SumBinBiovolume1_2 ./ ifcb_nass_morpho_2016_2023.TotalSum;

ifcb_nass_morpho_2016_2023.PropBinBiovolume3_4 = ifcb_nass_morpho_2016_2023.SumBinBiovolume3_4 ./ ifcb_nass_morpho_2016_2023.TotalSum;

ifcb_nass_morpho_2016_2023.PropBinBiovolume5_6 = ifcb_nass_morpho_2016_2023.SumBinBiovolume5_6 ./ ifcb_nass_morpho_2016_2023.TotalSum;



% Extract data for the year 2020

year_filter = (ifcb_nass_morpho_2016_2023.datetime.Year == 2020);

dayofyear_2020 = ifcb_nass_morpho_2016_2023.dayofyear(year_filter);



% Extract the proportions for the year 2020

prop1_2_2020 = ifcb_nass_morpho_2016_2023.PropBinBiovolume1_2(year_filter);

prop3_4_2020 = ifcb_nass_morpho_2016_2023.PropBinBiovolume3_4(year_filter);

prop5_6_2020 = ifcb_nass_morpho_2016_2023.PropBinBiovolume5_6(year_filter);



% Generate a full range of days for the year 2020

full_dayofyear = (min(dayofyear_2020):max(dayofyear_2020))';



% Replace Inf values with NaN in the proportions before interpolation

prop1_2_2020(prop1_2_2020 == Inf) = NaN;

prop3_4_2020(prop3_4_2020 == Inf) = NaN;

prop5_6_2020(prop5_6_2020 == Inf) = NaN;



% Interpolate missing data (NaN values are automatically handled by interp1)

interp_prop1_2 = interp1(dayofyear_2020, prop1_2_2020, full_dayofyear, 'linear', 'extrap');

interp_prop3_4 = interp1(dayofyear_2020, prop3_4_2020, full_dayofyear, 'linear', 'extrap');

interp_prop5_6 = interp1(dayofyear_2020, prop5_6_2020, full_dayofyear, 'linear', 'extrap');



% Prepare interpolated data for the stacked bar plot

interpolated_proportions_2020 = [interp_prop1_2, interp_prop3_4, interp_prop5_6];







%% New groupings

% Add new columns to the DataFrame

ifcb_nass_morpho_2016_2023.SumBinBiovolume1 = ifcb_nass_morpho_2016_2023.BinBiovolume1;

ifcb_nass_morpho_2016_2023.SumBinBiovolume2 = ifcb_nass_morpho_2016_2023.BinBiovolume2;

ifcb_nass_morpho_2016_2023.SumBinBiovolume3_4 = ifcb_nass_morpho_2016_2023.BinBiovolume3 + ifcb_nass_morpho_2016_2023.BinBiovolume4;

ifcb_nass_morpho_2016_2023.SumBinBiovolume5_6 = ifcb_nass_morpho_2016_2023.BinBiovolume5 + ifcb_nass_morpho_2016_2023.BinBiovolume6;



% Calculate total sum for each row

ifcb_nass_morpho_2016_2023.TotalSum = ifcb_nass_morpho_2016_2023.SumBinBiovolume1 + ...
                                      ifcb_nass_morpho_2016_2023.SumBinBiovolume2 + ...
                                      ifcb_nass_morpho_2016_2023.SumBinBiovolume3_4 + ...
                                      ifcb_nass_morpho_2016_2023.SumBinBiovolume5_6;



% Calculate proportions

ifcb_nass_morpho_2016_2023.PropBinBiovolume1 = ifcb_nass_morpho_2016_2023.SumBinBiovolume1 ./ ifcb_nass_morpho_2016_2023.TotalSum;

ifcb_nass_morpho_2016_2023.PropBinBiovolume2 = ifcb_nass_morpho_2016_2023.SumBinBiovolume2 ./ ifcb_nass_morpho_2016_2023.TotalSum;

ifcb_nass_morpho_2016_2023.PropBinBiovolume3_4 = ifcb_nass_morpho_2016_2023.SumBinBiovolume3_4 ./ ifcb_nass_morpho_2016_2023.TotalSum;

ifcb_nass_morpho_2016_2023.PropBinBiovolume5_6 = ifcb_nass_morpho_2016_2023.SumBinBiovolume5_6 ./ ifcb_nass_morpho_2016_2023.TotalSum;



% Extract data for the year 2020

year_filter = (ifcb_nass_morpho_2016_2023.datetime.Year == 2020);

dayofyear_2020 = ifcb_nass_morpho_2016_2023.dayofyear(year_filter);



% Extract the proportions for the year 2020

prop1_2020 = ifcb_nass_morpho_2016_2023.PropBinBiovolume1(year_filter);

prop2_2020 = ifcb_nass_morpho_2016_2023.PropBinBiovolume2(year_filter);

prop3_4_2020 = ifcb_nass_morpho_2016_2023.PropBinBiovolume3_4(year_filter);

prop5_6_2020 = ifcb_nass_morpho_2016_2023.PropBinBiovolume5_6(year_filter);



% Generate a full range of days for the year 2020

full_dayofyear = (min(dayofyear_2020):max(dayofyear_2020))';



% Interpolate missing data

interp_prop1 = interp1(dayofyear_2020, prop1_2020, full_dayofyear, 'linear', 'extrap');

interp_prop2 = interp1(dayofyear_2020, prop2_2020, full_dayofyear, 'linear', 'extrap');

interp_prop3_4 = interp1(dayofyear_2020, prop3_4_2020, full_dayofyear, 'linear', 'extrap');

interp_prop5_6 = interp1(dayofyear_2020, prop5_6_2020, full_dayofyear, 'linear', 'extrap');



% Prepare interpolated data for the stacked bar plot

interpolated_proportions_2020 = [interp_prop1, interp_prop2, interp_prop3_4, interp_prop5_6];



%% New groupings for Abundance

% Add new columns for summed abundance bins (similar to how it's done for Biovolume)

ifcb_nass_morpho_2016_2023.SumBinAbundance1 = ifcb_nass_morpho_2016_2023.BinAbundance1;

ifcb_nass_morpho_2016_2023.SumBinAbundance2 = ifcb_nass_morpho_2016_2023.BinAbundance2;

ifcb_nass_morpho_2016_2023.SumBinAbundance3_4 = ifcb_nass_morpho_2016_2023.BinAbundance3 + ifcb_nass_morpho_2016_2023.BinAbundance4;

ifcb_nass_morpho_2016_2023.SumBinAbundance5_6 = ifcb_nass_morpho_2016_2023.BinAbundance5 + ifcb_nass_morpho_2016_2023.BinAbundance6;



% Calculate total sum for each row (for abundance)

ifcb_nass_morpho_2016_2023.TotalSumAbundance = ifcb_nass_morpho_2016_2023.SumBinAbundance1 + ...
                                               ifcb_nass_morpho_2016_2023.SumBinAbundance2 + ...
                                               ifcb_nass_morpho_2016_2023.SumBinAbundance3_4 + ...
                                               ifcb_nass_morpho_2016_2023.SumBinAbundance5_6;



% Calculate proportions (for abundance)

ifcb_nass_morpho_2016_2023.PropBinAbundance1 = ifcb_nass_morpho_2016_2023.SumBinAbundance1 ./ ifcb_nass_morpho_2016_2023.TotalSumAbundance;

ifcb_nass_morpho_2016_2023.PropBinAbundance2 = ifcb_nass_morpho_2016_2023.SumBinAbundance2 ./ ifcb_nass_morpho_2016_2023.TotalSumAbundance;

ifcb_nass_morpho_2016_2023.PropBinAbundance3_4 = ifcb_nass_morpho_2016_2023.SumBinAbundance3_4 ./ ifcb_nass_morpho_2016_2023.TotalSumAbundance;

ifcb_nass_morpho_2016_2023.PropBinAbundance5_6 = ifcb_nass_morpho_2016_2023.SumBinAbundance5_6 ./ ifcb_nass_morpho_2016_2023.TotalSumAbundance;



%Extract abundance data for the year 2020

year_filter = (ifcb_nass_morpho_2016_2023.datetime.Year == 2020);

dayofyear_2020_abundance = ifcb_nass_morpho_2016_2023.dayofyear(year_filter);



% Extract the abundance proportions for the year 2020

prop1_abundance_2020 = ifcb_nass_morpho_2016_2023.PropBinAbundance1(year_filter);

prop2_abundance_2020 = ifcb_nass_morpho_2016_2023.PropBinAbundance2(year_filter);

prop3_4_abundance_2020 = ifcb_nass_morpho_2016_2023.PropBinAbundance3_4(year_filter);

prop5_6_abundance_2020 = ifcb_nass_morpho_2016_2023.PropBinAbundance5_6(year_filter);



% Generate a full range of days for the year 2020 (same as biovolume)

full_dayofyear_abundance = (min(dayofyear_2020_abundance):max(dayofyear_2020_abundance))';



% Interpolate missing abundance data

interp_prop1_abundance = interp1(dayofyear_2020_abundance, prop1_abundance_2020, full_dayofyear_abundance, 'linear', 'extrap');

interp_prop2_abundance = interp1(dayofyear_2020_abundance, prop2_abundance_2020, full_dayofyear_abundance, 'linear', 'extrap');

interp_prop3_4_abundance = interp1(dayofyear_2020_abundance, prop3_4_abundance_2020, full_dayofyear_abundance, 'linear', 'extrap');

interp_prop5_6_abundance = interp1(dayofyear_2020_abundance, prop5_6_abundance_2020, full_dayofyear_abundance, 'linear', 'extrap');



% Prepare interpolated data for the stacked bar plot (abundance)

interpolated_abundance_proportions_2020 = [interp_prop1_abundance, interp_prop2_abundance, interp_prop3_4_abundance, interp_prop5_6_abundance];





%% Figure 3. IFCB Mophometics

load('data\ifcb_processed\temp_processed\ifcb_nass_morpho_part1_1_2016_2023.mat')



% Bin edges and legend labels

% Define Bin Widths for ESD

bin_min = 15;

bin_max = 150;



num_bins = 7;



% Generate logarithmically spaced bin edges

binEdges_esd = round(logspace(log10(bin_min), log10(bin_max), num_bins));



% Update bin labels with adjusted upper bounds and rounded values

binLabels = cell(1, length(binEdges_esd)-1);

for i = 1:length(binEdges_esd)-1

    if i < length(binEdges_esd)-1

        % All bins except the last: [e1, e2 - 0.1], rounded

        lower_bound = round(binEdges_esd(i));

        upper_bound = round(binEdges_esd(i+1) - 1);

        binLabels{i} = sprintf('%d - %d', lower_bound, upper_bound);

    else

        % Last bin: [en-1, en], rounded

        lower_bound = round(binEdges_esd(i));

        upper_bound = round(binEdges_esd(i+1));

        binLabels{i} = sprintf('%d - %d', lower_bound, upper_bound);

    end

end





start_date = datetime(2016, 5, 1);

end_date = datetime(2023, 12, 1);

start_date_doy=day(start_date,'dayofyear');

end_date_doy=day(end_date,'dayofyear');



peak_fire = datetime(2020,8,21);

lag_days = 6;

start_date = datetime(2020, 8, 16);

end_date = datetime(2020, 9, 22);

lag_date = datetime(2020, 8, 21 + lag_days);

august_date = datetime(2020, 9, 11);



% Define days and datenum variables for shading

start_date_day = day(start_date, 'dayofyear');

end_date_day = day(end_date, 'dayofyear');

lag_date_day = day(lag_date, 'dayofyear');

peak_fire_day = day(peak_fire, 'dayofyear');

august_date_day = day(august_date, 'dayofyear');



% For shading

shade_start_day = start_date_day;

shade_end_day = end_date_day;



% Set up publication settings

ftsz = 8;  % Font size for publication

ftname = 'Helvetica';  % Font name

linewdt = 1.5;  % Line width for plotting







% Setup tiled layout

num_subplots = 5;

t = tiledlayout(num_subplots, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact');



% Define custom colors for the plot

colors_for_plot = {'#8bade8','#e69e5c','#73bd86'};



% Define a contrasting color palette

colors = [

    0.0000, 0.4470, 0.7410;  % Bright Blue

    0.8500, 0.3250, 0.0980;  % Bright Red

    0.9290, 0.6940, 0.1250;  % Bright Yellow

    0.4940, 0.1840, 0.5560;  % Purple

    0.4660, 0.6740, 0.1880;  % Bright Green

    0.3010, 0.7450, 0.9330;  % Cyan

    0.6350, 0.0780, 0.1840;  % Dark Red

    1.0000, 0.0000, 0.0000;  % Pure Red

    0.0000, 1.0000, 0.0000;  % Pure Green

    1.0000, 0.0000, 1.0000;  % Magenta

    0.0000, 1.0000, 1.0000;  % Bright Cyan

    1.0000, 1.0000, 0.0000;  % Bright Yellow

];





close all



% List of subplot letters

subplot_labels = {'a', 'b', 'c', 'd', 'e'};



% Setup tiled layout

tiledlayout(num_subplots, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact')





% Plot 1: Size-class proportions

nexttile

hBar = bar(full_dayofyear, interpolated_proportions_2020, 'stacked', 'BarWidth', 1.1);

hold on

colors_short = [0.0000, 0.4470, 0.7410;  % Bright Blue

            0.8500, 0.3250, 0.0980;  % Bright Red

          0.8509, 0.5922, 0.9098;

          0.7373, 0.8824, 0.4941];

for k = 1:length(hBar)

    hBar(k).FaceColor = colors_short(k, :);

end

hold on

% Add shaded area to the plot and add lines for 8 day lag and August Fire Smoke

xline(peak_fire_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':'); 

xline(lag_date_day, 'Color', '#cf572b', 'LineWidth', 1.5, 'LineStyle', ':'); 

% xline(august_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':');

xline(start_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 

fill([shade_start_day, shade_start_day, shade_end_day, shade_end_day], ...
    [min(ylim) * 1.5, max(ylim) * 1.5, max(ylim) * 1.5 min(ylim) * 1.5], ...
    [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off'); 

xline(end_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');



ylabel(["Proportion" + newline + "Total Biovolume"], 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold');

legend({'15-21 \mum','22-32 \mum','33-69 \mum','70-150 \mum'}, 'Location', 'bestoutside', 'FontSize', ftsz, 'AutoUpdate', 'off');

ylim([0 1])

xlim([start_date_doy end_date_doy])



% Set common x-axis properties

% xticks(xticks_days)

datetick('x','mmm','keeplimits')



text(-0.08, 1.25, subplot_labels{1}, 'Units', 'normalized', 'FontSize', 10, ...
     'FontName', ftname, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');





%Plot 2: Total Biovolume Anomaly

nexttile

%Multiply by 1000 to convert to um^3/L

conversion_factor=1e9;

conversion_factor=1;

shade_anomaly(movmean(ifcb_nass_morpho_2016_2023_join.dayofyear(year_filter),5), conversion_factor*movmean(ifcb_nass_morpho_2016_2023_join.anomaly_biovolumesum(year_filter),5))

% plot(movmean(ifcb_nass_morpho_2016_2023_join.dayofyear(year_filter),5), movmean(ifcb_nass_morpho_2016_2023_join.anomaly_biovolumesum(year_filter),5), ...

% '-','LineWidth', 1.5, 'Color', '#979da6');

hold on



% Add shaded area to the plot and add lines for 8 day lag and August Fire Smoke

xline(peak_fire_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':'); 

xline(lag_date_day, 'Color', '#cf572b', 'LineWidth', 1.5, 'LineStyle', ':'); 

% xline(august_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':');

xline(start_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 

xline(end_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');



yline(0, '--k', 'LineWidth', 1);

ylabel(["Biovolume" + newline + "(\mum^3 L^{-1})"], 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold');



% Set common x-axis properties

datetick('x','mmm','keeplimits')

xlim([start_date_doy end_date_doy])

ylim([-1.5e7*conversion_factor 3e7*conversion_factor])



fill([shade_start_day, shade_start_day, shade_end_day, shade_end_day], ...
    [min(ylim) * 1.5, max(ylim) * 1.5, max(ylim) * 1.5 min(ylim) * 1.5], ...
    [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off'); 

text(-0.08, 1.25, subplot_labels{2}, 'Units', 'normalized', 'FontSize', 10, ...
     'FontName', ftname, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');















% Plot 2: Total Sized Biovolume

nexttile

year_filter = (ifcb_nass_morpho_2016_2023.datetime.Year == 2020);

% yyaxis left

for i = 8:13

        %Multiply by 1000 to convert to mm^3/L

        conversion_factor=1e9;

        conversion_factor=1;

        plot(ifcb_nass_morpho_2016_2023.dayofyear(year_filter), conversion_factor*ifcb_nass_morpho_2016_2023{year_filter, i}, '-', ...
            'LineWidth', 1.5, 'Color', colors(i-7, :));

        ylim([0 4e7*conversion_factor])

        hold on

        ylabel(["Sized Biovolume" + newline + "(\mum^3 L^{-1})"], 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold');



end

legend(binLabels, 'Location', 'bestoutside', 'FontSize', ftsz-1, 'NumColumns', 1, ...
       'AutoUpdate', 'off');

title(legend, 'Size class [\mum]');



% Add shaded area to the plot and add lines for 6 day lag and August Fire Smoke

xline(peak_fire_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':'); 

xline(lag_date_day, 'Color', '#cf572b', 'LineWidth', 1.5, 'LineStyle', ':'); 

% xline(august_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':');



xline(start_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 

fill([shade_start_day, shade_start_day, shade_end_day, shade_end_day], ...
    [min(ylim) * 1.5, max(ylim) * 1.5, max(ylim) * 1.5 min(ylim) * 1.5], ...
    [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.2); 

xline(end_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');







xlim([start_date_doy end_date_doy])



% Set common x-axis properties

% xticks(xticks_days)

datetick('x','mmm','keeplimits')



text(-0.08, 1.25, subplot_labels{3}, 'Units', 'normalized', 'FontSize', 10, ...
     'FontName', ftname, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');







% Plot 3: Total Sized Abundance

nexttile

year_filter = (ifcb_nass_morpho_part1_1_2016_2023.datetime.Year == 2020);



for i = 14:19

    if i==14

        ax = gca; % Get current axes

        yyaxis left 

        %Multiply by 1000 to convert to /L

        plot(ifcb_nass_morpho_part1_1_2016_2023.dayofyear(year_filter), (ifcb_nass_morpho_part1_1_2016_2023{year_filter, i})*1000, '-', ...
            'LineWidth', 1.5, 'Color', colors(i-13, :));

        hold on

        ylabel(["cells L^{-1}" + newline + "(15-21 \mum)"], 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold','Color', 'k');

        ylim([0 2e6])

        fill([shade_start_day, shade_start_day, shade_end_day, shade_end_day], ...
        [min(ylim) * 2, max(ylim) * 2, max(ylim) * 2 min(ylim) * 2], ...
        [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off'); 

        ax.YColor = 'k'; % Set left y-axis tick labels and axis color to black



    else 

        yyaxis right

        ylabel(["cells L^{-1}" + newline + "(>21 \mum)"], 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold','Color', 'k');

        plot(ifcb_nass_morpho_part1_1_2016_2023.dayofyear(year_filter), (ifcb_nass_morpho_part1_1_2016_2023{year_filter, i})*1000, '-', ...
            'LineWidth', 1.5, 'Color', colors(i-13, :));

        ylim([0 5e5])

        ax.YColor = 'k'; % Set left y-axis tick labels and axis color to black



    end



    hold on

end

legend(binLabels, 'Location', 'bestoutside', 'FontSize', ftsz, 'AutoUpdate', 'off');

title(legend, 'Size class [\mum]');

legend('hide')



% Add shaded area to the plot and add lines for 8 day lag and August Fire Smoke

xline(peak_fire_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':'); 

xline(lag_date_day, 'Color', '#cf572b', 'LineWidth', 1.5, 'LineStyle', ':'); 

% xline(august_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':');

xline(start_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 

xline(end_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');





xlim([start_date_doy end_date_doy])

datetick('x','mmm','keeplimits')



text(-0.08, 1.25, subplot_labels{4}, 'Units', 'normalized', 'FontSize', 10, ...
     'FontName', ftname, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');





% Plot 4: NASS

nexttile

plot(ifcb_nass_morpho_2016_2023.dayofyear(year_filter), real(ifcb_nass_morpho_2016_2023.MeanSlopeAbundance(year_filter)), ...
    '-','LineWidth', 1.5, 'Color', colors_for_plot{1});

hold on



%Plot climatology from previous datatable

climatology_plot=(ifcb_nass_morpho_2016_2023_join(:,{'dayofyear';'MeanSlopeAbundance_climatology';'StandardError_MeanSlopeAbundance'}));



% Find unique rows based on 'MeanBiomassSum_climatology' column

[~, uniqueIdx] = unique(climatology_plot.MeanSlopeAbundance_climatology);



% Select only the unique rows

climatology_plot = climatology_plot(uniqueIdx, :);

climatology_plot = sortrows(climatology_plot,"dayofyear","ascend");

plot(climatology_plot.dayofyear, climatology_plot.MeanSlopeAbundance_climatology, 'LineWidth', 1.5, 'Color', '#979da6');

ci95 = 1.96 * climatology_plot.StandardError_MeanSlopeAbundance;

fill([climatology_plot.dayofyear; flipud(climatology_plot.dayofyear)], ...
    [climatology_plot.MeanSlopeAbundance_climatology + ci95; flipud(climatology_plot.MeanSlopeAbundance_climatology - ci95)], ...
    'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

ylabel(['NASS'], 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold');

legend({'2020','Climatology'}, 'Location', 'bestoutside','FontSize', ftsz, 'AutoUpdate', 'off');

ylim([-5 0])



% Add shaded area to the plot and add lines for 8 day lag and August Fire Smoke

xline(peak_fire_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':'); 

xline(lag_date_day, 'Color', '#cf572b', 'LineWidth', 1.5, 'LineStyle', ':'); 

% xline(august_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':');

xline(start_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 

fill([shade_start_day, shade_start_day, shade_end_day, shade_end_day], ...
    [min(ylim) * 1.5, max(ylim) * 1.5, max(ylim) * 1.5 min(ylim) * 1.5], ...
    [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off'); 

xline(end_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');



hold on

line_dates = [day(datetime(2020, 8, 16), 'dayofyear'), day(datetime(2020, 9, 22), 'dayofyear')];

line_heights = [max(ylim), max(ylim)];

bar(line_dates, line_heights, 'BarWidth', 0.005, 'FaceColor', 'k', 'EdgeColor', 'none');





xlim([start_date_doy end_date_doy])



% Set common x-axis properties

% xticks(xticks_days)

datetick('x','mmm','keeplimits')



text(-0.08, 1.25, subplot_labels{5}, 'Units', 'normalized', 'FontSize', 10, ...
     'FontName', ftname, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');









% Set figure size and save options

set(gcf, 'Units', 'inches', 'Position', [0, 0, 5, 6]);  % Adjusted figure size

set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 5, 6]);  % Match paper size exactly

set(gcf, 'PaperPositionMode', 'auto');  % Auto-scale to ensure full fit



% Saving the figure with high resolution for publication

saving = 0;

if saving == 1

    print(gcf, 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\final\figure_3_v0.png', '-dpng', '-r1200');

    print(gcf, 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\final\figure_3_v0.pdf', '-dpdf', '-r1200');

end



%% **Unused IFCB visual** Abundance porporitons

% Add a new plot using the abundance data (similar to the final tile but for abundance)

% Continue using the tiled layout from the previous code

close all

figure

% Setup tiled layout

t = tiledlayout(2, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact'); % Increase the number of subplots by 1



% (Insert the existing plots here, if necessary, before the new one)



% Plot 6: Size-class proportions using Abundance

nexttile

hBarAbundance = bar(full_dayofyear, fillmissing(interpolated_abundance_proportions_2020,"linear"), 'stacked', 'BarWidth', 1.1);

hold on

% Define the colors for the bars

colors_abundance = [0.0000, 0.4470, 0.7410;  % Bright Blue

                    0.8500, 0.3250, 0.0980;  % Bright Red

                    0.8509, 0.5922, 0.9098;  % Light Purple

                    0.7373, 0.8824, 0.4941]; % Light Green

for k = 1:length(hBarAbundance)

    hBarAbundance(k).FaceColor = colors_abundance(k, :);

end



% Add the same shading and fire smoke indicators

xline(peak_fire_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':'); 

% xline(lag_date_day, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 

xline(start_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 

fill([shade_start_day, shade_start_day, shade_end_day, shade_end_day], ...
    [min(ylim) * 1.5, max(ylim) * 1.5, max(ylim) * 1.5 min(ylim) * 1.5], ...
    [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off'); 

xline(end_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');



% Set ylabel and other formatting for abundance

ylabel(["Proportion" + newline + "Total Abundance"], 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold');

legend({'15-22 \mum','22-32 \mum','33-69 \mum','70-150 \mum'}, 'Location', 'bestoutside');

ylim([0 1])

xlim([start_date_doy end_date_doy])



% Set common x-axis properties

datetick('x','mmm','keeplimits')



nexttile

hBar = bar(full_dayofyear, interpolated_proportions_2020, 'stacked', 'BarWidth', 1.1);

hold on

colors = [0.0000, 0.4470, 0.7410;  % Bright Blue

            0.8500, 0.3250, 0.0980;  % Bright Red

          0.8509, 0.5922, 0.9098;

          0.7373, 0.8824, 0.4941];

for k = 1:length(hBar)

    hBar(k).FaceColor = colors(k, :);

end

% line_dates = [day(datetime(2020, 8, 18), 'dayofyear'), day(datetime(2020, 9, 22), 'dayofyear')];

% line_heights = [1, 1];

% bar(line_dates, line_heights, 'BarWidth', 0.005, 'FaceColor', 'k', 'EdgeColor', 'none');

hold on

% Add shaded area to the plot and add lines for 8 day lag and August Fire Smoke

xline(peak_fire_day, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':'); 

% xline(lag_date_day, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 

% xline(august_date_day, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');

xline(start_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 

fill([shade_start_day, shade_start_day, shade_end_day, shade_end_day], ...
    [min(ylim) * 1.5, max(ylim) * 1.5, max(ylim) * 1.5 min(ylim) * 1.5], ...
    [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off'); 

xline(end_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');



ylabel(["Proportion" + newline + "Total Biovolume"], 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold');

legend({'15-22 \mum','22-32 \mum','33-69 \mum','70-150 \mum'}, 'Location', 'bestoutside');

ylim([0 1])

xlim([start_date_doy end_date_doy])



% Set common x-axis properties

% xticks(xticks_days)

datetick('x','mmm','keeplimits')



% Format figure size for dual or single monitor use (same as earlier)

dual_monitor=0;

if dual_monitor == 1

    set(gcf, 'Position', [0 1100 1400 600])

else

    set(gcf, 'Position', [0 100 1400 600])

end



% Saving the figure

saving=0;

if saving == 1

    saveas(gcf, 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\ifcb\morphometrics_proportions_plots.png');

    saveas(gcf, 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\ifcb\morphometrics_proportions_plots.pdf');

end



%% IV. Figure 4. Cross-correlation & Taxanomic Changes



% Parameters for significance thresholding

numtests = 5000;  % Increased number of permutations for better p-value estimation

alpha = 0.05;  % Significance level



% Initialize containers for results for all taxa

all_taxa = {};

all_max_correlations = [];

all_max_lags = [];

all_p_values = [];

all_adjusted_p_values = [];



% Loop through each taxon

for i = 1:size(taxa_data, 2)

    % Get the current taxon data and PM2.5 data (after filling missing values)

    series1 = ifcb_pm25{:, taxa{i}};

    series2 = fillmissing(pm25, 'linear');

    

    % Compute the actual cross-correlation over the full lag range

    [cc, lags] = xcorr(series1, series2, 'coeff');

    

    % Filter to keep only lags between -10 and 30 (window of interest)

    lag_filter = (lags >= -10) & (lags <= 30);

    filtered_cc = cc(lag_filter);

    filtered_lags = lags(lag_filter);

    

    % Find the maximum cross-correlation and corresponding lag within the window of interest

    [max_cc, max_idx] = max(abs(filtered_cc));

    max_lag = filtered_lags(max_idx);

    

    % Perform random permutations using Fourier surrogate method (across all lags)

    max_rand_ccs = zeros(numtests, 1);

    for test = 1:numtests

        rand_series1 = generate_fourier_surrogate(series1);  % Generate Fourier surrogate for series1

        

        % Compute cross-correlation across all lags

        [xc_rand, ~] = xcorr(rand_series1, series2, 'coeff');

        

        % Now filter to keep only lags between -10 and 30 (same window of interest)

        rand_filtered_cc = xc_rand(lag_filter);

        

        % Store the maximum cross-correlation from the window of interest

        max_rand_ccs(test) = max(abs(rand_filtered_cc));

    end



    % Dynamic threshold based on the 95th percentile of null distribution

    dynamic_threshold = prctile(max_rand_ccs, 95);



    % Calculate p-value as the proportion of random max cross-correlations that exceeded the observed max

    tolerance = 1e-5;

    p_value = mean(max_rand_ccs >= (max_cc - tolerance));

    

    % Apply a small minimum threshold to avoid p-value = 0

    p_value = max(p_value, 1 / numtests);  % Ensure p-value is not exactly 0



    % Save the results for all taxa, whether significant or not

    all_taxa = [all_taxa, taxa{i}];  % Save taxon name

    all_max_correlations = [all_max_correlations, max_cc];  % Save max correlation

    all_max_lags = [all_max_lags, max_lag];  % Save corresponding lag

    all_p_values = [all_p_values, p_value];  % Save non-adjusted p-value

end







%% Now apply Benjamini-Hochberg FDR correction to all p-values

[sorted_p_values, sort_idx] = sort(all_p_values);  % Sort p-values in ascending order

num_taxa = length(all_p_values);  % Total number of taxa



% Adjust the p-values according to the Benjamini-Hochberg procedure

adjusted_p_values = sorted_p_values .* num_taxa ./ (1:num_taxa);

adjusted_p_values(adjusted_p_values > 1) = 1;  % Ensure adjusted p-values do not exceed 1



% Re-map the adjusted p-values back to their original order

unsorted_adjusted_p_values = zeros(size(all_p_values));

unsorted_adjusted_p_values(sort_idx) = adjusted_p_values;  % Remap sorted p-values to original order

all_adjusted_p_values = unsorted_adjusted_p_values;  % Store the adjusted p-values for all taxa



% Prepare the output with significance markings

taxa_with_marks = all_taxa;  % Copy of taxa names

for i = 1:num_taxa

    % Mark taxa with significant non-adjusted p-value

    if all_p_values(i) < alpha

        taxa_with_marks{i} = [taxa_with_marks{i}, '*'];  % Add * to the taxon name

    end

    % Mark taxa with significant adjusted p-value

    if all_adjusted_p_values(i) < alpha

        taxa_with_marks{i} = [taxa_with_marks{i}, '*'];  % Add ** to the taxon name

    end

end



% Round the p-values and max CC to 3 decimal places

rounded_p_values = round(all_p_values, 3);

rounded_adjusted_p_values = round(all_adjusted_p_values, 3);

rounded_max_correlations = round(all_max_correlations, 3);



% Create a table with the columns: taxa, max CC, lag, non-adjusted p-value, adjusted p-value

results_table = table(taxa_with_marks', rounded_max_correlations', all_max_lags', rounded_p_values', rounded_adjusted_p_values', ...
    'VariableNames', {'Taxa', 'Max_CC', 'Lag', 'P_value', 'Adjusted_p_value'});



% Sort the table by Max_CC in descending order

sorted_results_table = sortrows(results_table, 'Max_CC', 'descend');



% Save the sorted table to a CSV file

csv_filename = 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\tables\all_taxa_with_p_values.csv';

writetable(sorted_results_table, csv_filename);



% Display the sorted table to the user

disp('Results saved to all_taxa_with_p_values.csv');

disp(sorted_results_table);



% Plot the null distribution of maximum cross-correlations for visualization

figure;

histogram(all_max_rand_ccs, 30);

hold on;

xline(dynamic_threshold, 'r', 'Threshold (95th percentile)');

title('Null Distribution of Max Cross-Correlations');

xlabel('Max Cross-Correlation (Null)');

ylabel('Frequency');

legend('Random Surrogate CCs', 'Threshold');

hold off;





%% Now all taxa for plotting

numtests = 5000;  % Increased number of permutations for better p-value estimation

alpha = 0.05;  % Significance level



% Initialize containers for significant results

significant_taxa = {};

significant_p_values = [];

significant_max_correlations = [];

significant_max_lags = [];

all_max_rand_ccs = [];  % Store all random max cross-correlations for thresholding



% Store cross-correlation values for the heatmap

all_cc = [];  % Will store cross-correlation values for all taxa

all_lags = lags;  % Store all possible lags (shared across taxa)

significant_mask = zeros(length(all_lags), size(taxa_data, 2));  % Initialize mask for significant values



% Loop through each taxon

for i = 1:size(taxa_data, 2)

    % Get the current taxon data and PM2.5 data (after filling missing values)

    series1 = ifcb_pm25{:, taxa{i}};

    series2 = fillmissing(pm25, 'linear');

    

    % Compute the actual cross-correlation over the full lag range

    [cc, lags] = xcorr(series1, series2, 'coeff');

    

    % Store the cross-correlation values for this taxon

    all_cc(:, i) = cc;  % Add to all cross-correlation values for the heatmap

    

    % Filter to keep only lags between -10 and 30 (window of interest)

    lag_filter = (lags >= -10) & (lags <= 30);

    filtered_cc = cc(lag_filter);

    filtered_lags = lags(lag_filter);

    

    % Find the maximum cross-correlation and corresponding lag within the window of interest

    [max_cc, max_idx] = max(abs(filtered_cc));

    max_lag = filtered_lags(max_idx);

    

    % Perform random permutations using Fourier surrogate method (across all lags)

    max_rand_ccs = zeros(numtests, 1);

    for test = 1:numtests

        rand_series1 = generate_fourier_surrogate(series1);  % Generate Fourier surrogate for series1

        

        % Compute cross-correlation across all lags

        [xc_rand, ~] = xcorr(rand_series1, series2, 'coeff');

        

        % Now filter to keep only lags between -10 and 30 (same window of interest)

        rand_filtered_cc = xc_rand(lag_filter);

        

        % Store the maximum cross-correlation from the window of interest

        max_rand_ccs(test) = max(abs(rand_filtered_cc));

    end

    all_max_rand_ccs = [all_max_rand_ccs; max_rand_ccs];  % Collect max values from all tests



    % Dynamic threshold based on the 95th percentile of null distribution

    dynamic_threshold = prctile(max_rand_ccs, 95);

    

    % Only consider the taxon if its max CC exceeds the dynamic threshold

    if max_cc > dynamic_threshold

        % Calculate p-value as the proportion of random max cross-correlations that exceeded the observed max

        tolerance = 1e-5;

        p_value = mean(max_rand_ccs >= (max_cc - tolerance));

        

        % Apply a small minimum threshold to avoid p-value = 0

        p_value = max(p_value, 1 / numtests);  % Ensure p-value is not exactly 0

        

        % Check if the p-value is significant

        if p_value < alpha

            % Mark the cross-correlation as significant in the mask

            significant_mask(lag_filter, i) = (abs(filtered_cc) > dynamic_threshold);  % Mark * for significant points

        end

    end

end



%% Figure 4 Final. Combined Heatmap and PDF Plots Side-by-Side

close all;



subplot_labels = {'\bf{b}','\bf{c}','\bf{d}','\bf{e}','\bf{f}','\bf{g}','\bf{h}'};



italic_taxa = {

    'Akashiwo', 'Asterionellopsis','Boreadinium', 'Ceratium','Chaetoceros', ...

    'Cochlodinium', 'Corethron', 'Cyl Nitz', ...

    'Det Cer Lau', 'Dictyocha', 'Dinophysis', ...

    'Ditylum', 'Entomoneis', 'Eucampia', ...

    'Gymnodinium','Gyrodinium', 'Hemiaulus', 'Leptocylindrus', ...

    'Peridinium','Phaeocystis', 'Polykrikos','Protocentrum','Protoperidinium', 'Pseudo nitzschia', ...

    'Pyramimonas','Skeletonema', 'Thalassionema', 'Thalassiosira', ...

    'Tiarina', 'Tontonia', 'Torodinium', ...

    'Tropidoneis', 'Vicicitus'

};







ftsz = 8;

figure('Units', 'inches', 'Position', [0, 0, 5, 6]);



% --- Left plot: Heatmap using imagesc --- %

subplot('Position', [0.22, 0.12, 0.42, 0.8]);  % Wider and aligned with PDFs



imagesc(lags, 1:numel(taxa), data);

colormap(jet);

caxis([0, 0.8]);

xlim([0 15]);



% Taxa label formatting

taxa_clean = strrep(taxa, '_', ' ');

set(gca, 'YTick', 1:numel(taxa_clean));

set(gca, 'YTickLabel', taxa_clean);



% Italicize only taxa in italic_taxa

set(gca, 'TickLabelInterpreter', 'tex');  % Switch from LaTeX back to default interpreter

yticklabels = taxa_clean;



for i = 1:numel(taxa_clean)

    base_name = strrep(taxa_clean{i}, '*', '');

    

    if ismember(base_name, italic_taxa)

        yticklabels{i} = ['\it ', taxa_clean{i}];  % TeX-style italics

    end

end



set(gca, 'YTick', 1:numel(taxa_clean));

set(gca, 'YTickLabel', yticklabels);

set(gca, 'TickLabelInterpreter', 'tex');  % Use TeX, which respects \it

set(gca, 'FontName', ftname);  % Restore your figure's font

set(gca, 'FontSize', ftsz);  % Restore your figure's font





% Axes styling

ax = gca;

ax.YAxis.FontSize = ftsz;

ax.YAxis.FontName = ftname;

ax.XAxis.FontSize = ftsz;

ax.XAxis.FontName = ftname;

ax.XTick = lags;

ax.XTickLabelRotation = 0;

ax.XTickLabels=CustomXLabels;

xlabel('Lag (days after PM 2.5 anomaly)', 'FontSize', ftsz, 'FontName', ftname);



% Top colorbar aligned with left plot

cb = colorbar('Location', 'northoutside');

cb.Position = [0.21, 0.93, 0.42, 0.015];

cb.FontSize = ftsz;

cb.FontName = ftname;



cb.Label.String = '\it Cross correlation';  % Bold italic C

cb.Label.FontSize = ftsz;

cb.Label.FontName = ftname;

cb.Label.Interpreter = 'tex';  % Use TeX for bold/italic formatting

cb.Label.Rotation = 0;         % Keep horizontal

cb.Label.HorizontalAlignment = 'center';

cb.Label.VerticalAlignment = 'middle';



% Shift label to the right of the colorbar

% cb.Label.Position(1) = cb.Position(1) + cb.Position(3) + 0.1;  % x + width + small offset

cb.Label.Position(2) = cb.Position(2)*4.5;    % slightly above vertical center



% Subplot label 'a'

annotation('textbox', [0.08, 0.92, 0.03, 0.05], ...
           'String', '\bf{a}', ...
           'FontSize', ftsz + 2, ...
           'FontName', ftname, ...
           'LineStyle', 'none', ...
           'VerticalAlignment', 'bottom');





% --- Right plot: PDFs --- %

num_taxa = numel(unique_taxa);

subplot_height = 0.8 / num_taxa;



for ii = 1:num_taxa

    subplot('Position', [0.68, 0.12 + (num_taxa - ii)*subplot_height, 0.2, subplot_height - 0.01]);



    % Density for full dataset

    data_full = log10(ifcbbray_dates_filt.(unique_taxa{ii}));

    data_full = data_full(~isnan(data_full) & ~isinf(data_full));

    [f_full, x_full] = ksdensity(data_full);

    plot(x_full, f_full, 'LineWidth', 2, 'Color', 'black'); hold on;

    area(x_full, f_full, 'FaceColor', 'black', 'FaceAlpha', 0.5);



    % Density for windowed subset

    data_window = log10(ifcb_2020_sel.(unique_taxa{ii}));

    data_window = data_window(~isnan(data_window) & ~isinf(data_window));

    [f_window, x_window] = ksdensity(data_window);

    plot(x_window, f_window, 'LineWidth', 2, 'Color', colorz(ii, :));

    area(x_window, f_window, 'FaceColor', colorz(ii, :), 'FaceAlpha', 0.7);



    xlim([-2 8]);

    ylim([0 0.9]);

    yticks([0.3, 0.7]);



    % Only add x-label to bottom plot

    if ii == num_taxa

        xlabel('Log_{10} cells L^{-1}', ...
               'FontSize', ftsz, ...
               'FontName', ftname);

    else

        set(gca, 'XTickLabel', []);  % Remove x-tick labels

    end



    % Taxa name (top-left inside plot)

    x_text = min(xlim) + 0.02 * range(xlim);

    y_text = min(ylim) + 0.95 * range(ylim);

    if ismember(unique_taxa{ii}, italic_taxa)

        text(x_text, y_text, unique_taxa{ii}, ...
             'FontSize', ftsz, ...
             'FontName', ftname, ...
             'FontAngle', 'italic', ...
             'HorizontalAlignment', 'left', ...
             'VerticalAlignment', 'top');

    else

        text(x_text, y_text, unique_taxa{ii}, ...
             'FontSize', ftsz, ...
             'FontName', ftname, ...
             'HorizontalAlignment', 'left', ...
             'VerticalAlignment', 'top');

    end



    % Subplot label (outside top-left)

    pos = get(gca, 'Position');

    annotation('textbox', [pos(1) - 0.04, pos(2) + pos(4) - 0.005, 0.03, 0.03], ...
               'String', subplot_labels{ii}, ...
               'FontWeight', 'bold', ...
               'FontSize', ftsz + 2, ...
               'FontName', ftname, ...
               'LineStyle', 'none', ...
               'VerticalAlignment', 'top');



    % Axes formatting

    ax = gca;

    ax.FontSize = ftsz;

    ax.FontName = ftname;

    ax.XAxisLocation = 'bottom';

    ax.XTickLabelRotation = 0;

    ax.XDir = 'normal';

    ax.YAxisLocation = 'right';

end



% --- Shared y-axis label on right --- %

han = axes(gcf, 'Visible', 'off', 'Position', [0.98, 0.14, 0.01, 0.76]);

han.YLabel.Visible = 'on';

ylabel(han, shared_y_label, ...
       'FontWeight','bold', ...
       'FontSize', ftsz, ...
       'FontName', ftname, ...
       'Rotation', -90, ...
       'HorizontalAlignment', 'center', ...
       'VerticalAlignment', 'bottom');



% Export options

set(gcf, 'PaperPositionMode', 'auto');



saving=0;

if saving==1

    print(gcf, '-dpng', '-r1200', "C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\final\figure5_v1.png");

    print(gcf, '-dpdf', '-r1200', "C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\final\figure5_v1.pdf");

end





