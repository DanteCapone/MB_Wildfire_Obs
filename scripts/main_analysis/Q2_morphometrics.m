%IFCB Analysis: 
%% Script to conduct morphometric analysis on processed data

%Add paths
addpath('Q:\Dante\Wildfire_Obs\functions\')
addpath('Q:\Dante\Wildfire_Obs\')
addpath(genpath('Q:\Dante\data\MB_Wildfire_Obs\processed_data'))

load('ifcb_nass_morpho_all_2016_2023.mat')




%%
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



% %% Plot all data all time
% 
% ftsz=16;
% xtick_fontsize=14;
% clf;
% figure(2);
% hold on;
% 
% % Define the start and end dates for x-axis limits
% start_date = datetime(2016, 1, 1);
% end_date = datetime(2023, 12, 31);
% start_date_doy=day(start_date,'dayofyear');
% end_date_doy=day(end_date,'dayofyear');
% xticks_days=start_date:caldays(10):end_date;
% xticks_days=start_date:caldays(10):end_date;
% 
% 
% % Define the shading window
% shade_start = datetime(2020,8,16);
% shade_end = datetime(2020,9,22);
% shade_start_doy = day(shade_start, 'dayofyear');
% shade_end_doy = day(shade_end, 'dayofyear');
% 
% % Bin edges and legend labels
% % Define Bin Widths for ESD
% bin_min = 15;
% bin_max = 150;
% 
% num_bins = 7;
% 
% % Generate logarithmically spaced bin edges
% binEdges_esd = round(logspace(log10(bin_min), log10(bin_max), num_bins), 1);
% 
% % Update bin labels with adjusted upper bounds
% binLabels = cell(1, length(binEdges_esd)-1);
% for i = 1:length(binEdges_esd)-1
%     if i < length(binEdges_esd)-1
%         % All bins except the last: [e1, e2 - 0.1]
%         binLabels{i} = sprintf('%.1f - %.1f', binEdges_esd(i), binEdges_esd(i+1) - 0.1);
%     else
%         % Last bin: [en-1, en]
%         binLabels{i} = sprintf('%.1f - %.1f', binEdges_esd(i), binEdges_esd(i+1));
%     end
% end
% 
% % Plot BinBiovolume1 to BinBiovolume9
% % Define a contrasting color palette
% colors = [
%     0.0000, 0.4470, 0.7410;  % Bright Blue
%     0.8500, 0.3250, 0.0980;  % Bright Red
%     0.9290, 0.6940, 0.1250;  % Bright Yellow
%     0.4940, 0.1840, 0.5560;  % Purple
%     0.4660, 0.6740, 0.1880;  % Bright Green
%     0.3010, 0.7450, 0.9330;  % Cyan
%     0.6350, 0.0780, 0.1840;  % Dark Red
%     1.0000, 0.0000, 0.0000;  % Pure Red
%     0.0000, 1.0000, 0.0000;  % Pure Green
%     1.0000, 0.0000, 1.0000;  % Magenta
%     0.0000, 1.0000, 1.0000;  % Bright Cyan
%     1.0000, 1.0000, 0.0000;  % Bright Yellow
% ];
% 
% subplot(3,1,1)
% for i = 9:14
%     colorIndex = mod(i-9, size(colors, 1)) + 1;
%     plot(ifcb_nass_morpho_2016_2023.datetime, ifcb_nass_morpho_2016_2023{:, i}, '-', 'LineWidth', 2, 'Color', colors(colorIndex, :));
%     hold on
% end
% 
% % xlabel('Date','FontSize',ftsz,'FontName',ftname,'FontWeight','bold')
% ylabel('Total Sized Biovolume (um^3)','FontSize',ftsz,'FontName',ftname,'FontWeight','bold')
% 
% 
% % Add legend
% legend(binLabels, 'Location', 'best', 'FontSize', ftsz-4,'AutoUpdate', 'off');
% title(legend, 'Size class [um]');
% 
% % Add shaded area to the plot
% fill([shade_start, shade_start, shade_end, shade_end], ...
%      [min(ylim) 5e7 5e7 min(ylim)], ...
%      [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
% ylim([0 4e7])
% set(gca, 'FontSize', xtick_fontsize)  % Set x-tick label size
% 
% xlim([start_date end_date])
% % xticks(xticks_days)
% % xticklabels(datestr(xticks_days))
% % 
% datetick('x','mm/yyyy','keeplimits','keepticks')
% 
% 
% subplot(3,1,2)
% yyaxis left
% plot(ifcb_nass_morpho_2016_2023.datetime, ifcb_nass_morpho_2016_2023.MeanESD, '-','LineWidth',2,'Color','#8bade8');
% hold on
% 
% %xlabel('Date','FontSize',ftsz,'FontName',ftname,'FontWeight','bold')
% ylabel("Mean Equivalent Spherical"+newline+"Diameter (um)",'FontSize',ftsz,'FontName',ftname,'FontWeight','bold')
% 
% % Add mean line
% meanESD = mean(ifcb_nass_morpho_2016_2023.MeanESD);
% line([start_date end_date], [meanESD meanESD], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5, 'Color', '#8bade8');
% 
% % Add shaded area to the plot
% fill([shade_start, shade_start, shade_end, shade_end], ...
%      [0 40 40 0], ...
%      [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
% set(gca, 'FontSize', xtick_fontsize)  % Set x-tick label size
% hold on
% yyaxis right
% plot(ifcb_nass_morpho_2016_2023.datetime, ifcb_nass_morpho_2016_2023.MeanBiomassSum, '-','LineWidth',2,'Color','#e69e5c');
% ylabel(["Total Community"+newline+"Biovolume (um^3)"],'FontSize',ftsz,'FontName',ftname,'FontWeight','bold')
% 
% % Add mean line
% meanbiovol= mean(ifcb_nass_morpho_2016_2023.MeanBiomassSum);
% line([start_date end_date], [meanbiovol meanbiovol], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5, 'Color', '#e69e5c');
% 
% xlim([start_date end_date])
% % xticks(xticks_days)
% % xticklabels(datestr(xticks_days))
% % 
% datetick('x','mm/yyyy','keeplimits','keepticks')
% 
% subplot(3,1,3)
% plot((ifcb_nass_morpho_2016_2023.datetime), real(ifcb_nass_morpho_2016_2023.MeanSlopeAbundance), '-','LineWidth',2,'Color','#73bd86');
% hold on
% xlim([start_date end_date])
% % xticks(xticks_days)
% % xticklabels(datestr(xticks_days))
% % 
% datetick('x','mm/yyyy','keeplimits','keepticks')
% xlabel('Date','FontSize',ftsz,'FontName',ftname,'FontWeight','bold')
% ylabel(['Normalized Abundance',newline,'Size Spectrum Slope'],'FontSize',ftsz,'FontName',ftname,'FontWeight','bold')
% 
% % Add mean line
% meanSlopeAbundance = mean(ifcb_nass_morpho_2016_2023.MeanSlopeAbundance);
% line([start_date end_date], [meanSlopeAbundance meanSlopeAbundance], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5, 'Color', [0 0 0 0.5]);
% % 
% % % Add shaded area to the plot
% fill([shade_start, shade_start, shade_end, shade_end], ...
%      [-6 2 2 -6], ...
%      [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
% set(gca, 'FontSize', xtick_fontsize)  % Set x-tick label size
% 
% 
% dual_monitor=1;
% if dual_monitor==1
%     set(gcf,'Position',[-400 1400 1600 1200])
% else
%     set(gcf,'Position',[0 100 1600 1200])
% end
% 
% 
% 
% 
% saving = 0;
% if saving == 1
%     saveas(gcf, 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\ifcb\ifcb_nass_2020.png');
%     saveas(gcf, 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\ifcb\ifcb_nass_2020.pdf');
% end
% 
% 
% 
% 
% %% Now by day of year
% 
% ifcb_nass_morpho_2016_2023=ifcb_nass_morpho_2016_2023;
% ftsz = 16;
% xtick_fontsize = 14;
% clf;
% figure(2);
% hold on;
% 
% % Calculate DayOfYear and add to ifcb_nass_morpho_2016_2023 DataFrame
% ifcb_nass_morpho_2016_2023.DayOfYear = day(ifcb_nass_morpho_2016_2023.datetime, 'dayofyear');
% 
% % Get unique years in the data
% unique_years = year(ifcb_nass_morpho_2016_2023.datetime);
% 
% % Bin edges and legend labels
% binEdges_esd = round(logspace(log10(bin_min), log10(bin_max), num_bins), 1);
% binLabels = cell(1, length(binEdges_esd) - 1);
% for i = 1:length(binEdges_esd) - 1
%     binLabels{i} = sprintf('%.1f - %.1f', binEdges_esd(i), binEdges_esd(i + 1));
% end
% 
% % Define a contrasting color palette
% colors = lines(length(unique(unique_years))); % Generate enough colors for each year
% 
% subplot(3,1,1)
% for y = 1:length(unique(unique_years))
%     year_filter = unique_years == unique(unique_years(y));
%     for i = 9:14
%         plot(ifcb_nass_morpho_2016_2023.DayOfYear(year_filter), ifcb_nass_morpho_2016_2023{year_filter, i}, '-', ...
%             'LineWidth', 2, 'Color', colors(y, :));
%         hold on
%     end
% end
% 
% ylabel('Total Sized Biovolume (um^3)', 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold');
% 
% % Add legend
% legend(binLabels, 'Location', 'best', 'FontSize', ftsz, 'AutoUpdate', 'off');
% title(legend, 'Size class [um]');
% 
% % Add shaded area to the plot (adjust if necessary based on DayOfYear)
% 
% fill([shade_start_doy, shade_start_doy, shade_end_doy, shade_end_doy], ...
%      [min(ylim) 5e7 5e7 min(ylim)], ...
%      [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
% ylim([0 4e7])
% set(gca, 'FontSize', xtick_fontsize)
% 
% xlim([1 366])
% xticks(0:30:365)
% xticklabels(datestr(datenum(0:30:365), 'dd-mmm'))
% 
% subplot(3,1,2)
% yyaxis left
% for y = 1:length(unique(unique_years))
%     year_filter = unique_years == unique(unique_years(y));
%     plot(ifcb_nass_morpho_2016_2023.DayOfYear(year_filter), ifcb_nass_morpho_2016_2023.MeanESD(year_filter), ...
%         '-','LineWidth', 2, 'Color', colors(y, :));
%     hold on
% end
% 
% ylabel("Mean Equivalent Spherical" + newline + "Diameter (um)", 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold');
% 
% % Add shaded area to the plot (adjust if necessary based on DayOfYear)
% fill([shade_start_doy, shade_start_doy, shade_end_doy, shade_end_doy], ...
%      [0 40 40 0], ...
%      [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
% set(gca, 'FontSize', xtick_fontsize)
% hold on
% yyaxis right
% for y = 1:length(unique(unique_years))
%     year_filter = unique_years == unique(unique_years(y));
%     plot(ifcb_nass_morpho_2016_2023.DayOfYear(year_filter), ifcb_nass_morpho_2016_2023.MeanBiomassSum(year_filter), ...
%         '-','LineWidth', 2, 'Color', colors(y, :));
%     hold on
% end
% ylabel(["Total Community" + newline + "Biovolume (um^3)"], 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold');
% 
% xlim([1 366])
% xticks(0:30:365)
% xticklabels(datestr(datenum(0:30:365), 'dd-mmm'))
% 
% subplot(3,1,3)
% for y = 1:length(unique(unique_years))
%     year_filter = unique_years == unique(unique_years(y));
%     plot(ifcb_nass_morpho_2016_2023.DayOfYear(year_filter), real(ifcb_nass_morpho_2016_2023.MeanSlopeAbundance(year_filter)), ...
%         '-','LineWidth', 2, 'Color', colors(y, :));
%     hold on
% end
% 
% xlim([1 366])
% xticks(0:30:365)
% xticklabels(datestr(datenum(0:30:365), 'dd-mmm'))
% xlabel('Day of Year', 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold')
% ylabel(['Normalized Abundance', newline, 'Size Spectrum Slope'], 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold')
% 
% % Add shaded area to the plot (adjust if necessary based on DayOfYear)
% fill([shade_start_doy, shade_start_doy, shade_end_doy, shade_end_doy], ...
%      [-6 2 2 -6], ...
%      [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
% set(gca, 'FontSize', xtick_fontsize)
% 
% % Positioning
% dual_monitor = 0;
% if dual_monitor == 1
%     set(gcf,'Position',[-400 1200 1200 1000])
% else
%     set(gcf,'Position',[0 100 1600 1200])
% end
% 
% % Saving
% saving = 0;
% if saving == 1
%     saveas(gcf, 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\ifcb\ifcb_nass_2020.png');
%     saveas(gcf, 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\ifcb\ifcb_nass_2020.pdf');
% end

%% Calculate necessary Climatologies
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


%% Final Plot for the 2020 window with 5 day moving mean window
load('ifcb_nass_morpho_part1_1_2016_2023.mat')

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
ftsz = 10;  % Font size for publication
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
xline(peak_fire_day, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':'); 
xline(lag_date_day, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% xline(august_date_day, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');
xline(start_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
fill([shade_start_day, shade_start_day, shade_end_day, shade_end_day], ...
    [min(ylim) * 1.5, max(ylim) * 1.5, max(ylim) * 1.5 min(ylim) * 1.5], ...
    [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off'); 
xline(end_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');

ylabel(["Proportion" + newline + "Total Biovolume"], 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold');
legend({'15-21 \mum','22-32 \mum','33-69 \mum','70-150 \mum'}, 'Location', 'bestoutside', 'FontSize', ftsz, 'AutoUpdate', 'off');
ylim([0 1])
xlim([start_date_doy end_date_doy])

% Set common x-axis properties
% xticks(xticks_days)
datetick('x','mmm','keeplimits')

text(0.01, 0.98, subplot_labels{1}, 'Units', 'normalized', 'FontSize', 16, ...
     'FontName', ftname, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');


%Plot 2: Total Biovolume Anomaly
nexttile
%Multiply by 1000 to convert to mm^3/L
conversion_factor=1000/1e9;
conversion_factor=1;
shade_anomaly(movmean(ifcb_nass_morpho_2016_2023_join.dayofyear(year_filter),5), conversion_factor*movmean(ifcb_nass_morpho_2016_2023_join.anomaly_biovolumesum(year_filter),5))
% plot(movmean(ifcb_nass_morpho_2016_2023_join.dayofyear(year_filter),5), movmean(ifcb_nass_morpho_2016_2023_join.anomaly_biovolumesum(year_filter),5), ...
% '-','LineWidth', 2, 'Color', '#979da6');
hold on

% Add shaded area to the plot and add lines for 8 day lag and August Fire Smoke
xline(peak_fire_day, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':'); 
xline(lag_date_day, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% xline(august_date_day, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');
xline(start_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
xline(end_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');

yline(0, '--k', 'LineWidth', 1);
ylabel(["Community Biovolume" + newline + "Anomaly (\mum^3 mL^{-1})"], 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold');

% Set common x-axis properties
datetick('x','mmm','keeplimits')
xlim([start_date_doy end_date_doy])
ylim([-1.5e7*conversion_factor 3e7*conversion_factor])

fill([shade_start_day, shade_start_day, shade_end_day, shade_end_day], ...
    [min(ylim) * 1.5, max(ylim) * 1.5, max(ylim) * 1.5 min(ylim) * 1.5], ...
    [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off'); 
text(0.01, 0.98, subplot_labels{2}, 'Units', 'normalized', 'FontSize', 16, ...
     'FontName', ftname, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');







% Plot 2: Total Sized Biovolume
nexttile
year_filter = (ifcb_nass_morpho_2016_2023.datetime.Year == 2020);
% yyaxis left
for i = 8:13
        %Multiply by 1000 to convert to mm^3/L
        conversion_factor=1000/1e9;
        conversion_factor=1;
        plot(ifcb_nass_morpho_2016_2023.dayofyear(year_filter), conversion_factor*ifcb_nass_morpho_2016_2023{year_filter, i}, '-', ...
            'LineWidth', 2, 'Color', colors(i-7, :));
        ylim([0 4e7*conversion_factor])
        hold on
        ylabel(["Total Sized" + newline + "Biovolume (\mum^3 mL^{-1})"], 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold');

end
legend(binLabels, 'Location', 'bestoutside', 'FontSize', ftsz, 'AutoUpdate', 'off');
title(legend, 'Size class [\mum]');

% Add shaded area to the plot and add lines for 6 day lag and August Fire Smoke
xline(peak_fire_day, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':'); 
xline(lag_date_day, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% xline(august_date_day, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');

xline(start_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
fill([shade_start_day, shade_start_day, shade_end_day, shade_end_day], ...
    [min(ylim) * 1.5, max(ylim) * 1.5, max(ylim) * 1.5 min(ylim) * 1.5], ...
    [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4); 
xline(end_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');



xlim([start_date_doy end_date_doy])

% Set common x-axis properties
% xticks(xticks_days)
datetick('x','mmm','keeplimits')

text(0.01, 0.98, subplot_labels{3}, 'Units', 'normalized', 'FontSize', 16, ...
     'FontName', ftname, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');



% Plot 3: Total Sized Abundance
nexttile
year_filter = (ifcb_nass_morpho_part1_1_2016_2023.datetime.Year == 2020);

for i = 14:19
    if i==14
        ax = gca; % Get current axes
        yyaxis left 
        %Multiply by 1000 to convert to /L
        plot(ifcb_nass_morpho_part1_1_2016_2023.dayofyear(year_filter), (ifcb_nass_morpho_part1_1_2016_2023{year_filter, i}), '-', ...
            'LineWidth', 2, 'Color', colors(i-13, :));
        hold on
        ylabel('Cells mL^{-1} (15-21 \mum)', 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold','Color', 'k');
        ylim([0 2e3])
        fill([shade_start_day, shade_start_day, shade_end_day, shade_end_day], ...
        [min(ylim) * 2, max(ylim) * 2, max(ylim) * 2 min(ylim) * 2], ...
        [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off'); 
        ax.YColor = 'k'; % Set left y-axis tick labels and axis color to black

    else 
        yyaxis right
        ylabel(["Cells mL^{-1}" + newline + "(>21 \mum)"], 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold','Color', 'k');
        plot(ifcb_nass_morpho_part1_1_2016_2023.dayofyear(year_filter), (ifcb_nass_morpho_part1_1_2016_2023{year_filter, i}), '-', ...
            'LineWidth', 2, 'Color', colors(i-13, :));
        ylim([0 5e2])
        ax.YColor = 'k'; % Set left y-axis tick labels and axis color to black

    end

    hold on
end
legend(binLabels, 'Location', 'bestoutside', 'FontSize', ftsz, 'AutoUpdate', 'off');
title(legend, 'Size class [\mum]');

% Add shaded area to the plot and add lines for 8 day lag and August Fire Smoke
xline(peak_fire_day, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':'); 
xline(lag_date_day, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% xline(august_date_day, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');
xline(start_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
xline(end_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');


xlim([start_date_doy end_date_doy])
datetick('x','mmm','keeplimits')

text(0.01, 0.98, subplot_labels{4}, 'Units', 'normalized', 'FontSize', 16, ...
     'FontName', ftname, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');


% Plot 4: NASS
nexttile
plot(ifcb_nass_morpho_2016_2023.dayofyear(year_filter), real(ifcb_nass_morpho_2016_2023.MeanSlopeAbundance(year_filter)), ...
    '-','LineWidth', 2, 'Color', colors_for_plot{1});
hold on

%Plot climatology from previous datatable
climatology_plot=(ifcb_nass_morpho_2016_2023_join(:,{'dayofyear';'MeanSlopeAbundance_climatology';'StandardError_MeanSlopeAbundance'}));

% Find unique rows based on 'MeanBiomassSum_climatology' column
[~, uniqueIdx] = unique(climatology_plot.MeanSlopeAbundance_climatology);

% Select only the unique rows
climatology_plot = climatology_plot(uniqueIdx, :);
climatology_plot = sortrows(climatology_plot,"dayofyear","ascend");
plot(climatology_plot.dayofyear, climatology_plot.MeanSlopeAbundance_climatology, 'LineWidth', 2, 'Color', '#979da6');
ci95 = 1.96 * climatology_plot.StandardError_MeanSlopeAbundance;
fill([climatology_plot.dayofyear; flipud(climatology_plot.dayofyear)], ...
    [climatology_plot.MeanSlopeAbundance_climatology + ci95; flipud(climatology_plot.MeanSlopeAbundance_climatology - ci95)], ...
    'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
ylabel(['Normalized Abundance', newline, 'Size Spectrum Slope'], 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold');
legend({'2020','Climatology'}, 'Location', 'bestoutside','FontSize', ftsz, 'AutoUpdate', 'off');
ylim([-5 0])

% Add shaded area to the plot and add lines for 8 day lag and August Fire Smoke
xline(peak_fire_day, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':'); 
xline(lag_date_day, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
% xline(august_date_day, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');
xline(start_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
fill([shade_start_day, shade_start_day, shade_end_day, shade_end_day], ...
    [min(ylim) * 1.5, max(ylim) * 1.5, max(ylim) * 1.5 min(ylim) * 1.5], ...
    [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off'); 
xline(end_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');

hold on
line_dates = [day(datetime(2020, 8, 16), 'dayofyear'), day(datetime(2020, 9, 22), 'dayofyear')];
line_heights = [max(ylim), max(ylim)];
bar(line_dates, line_heights, 'BarWidth', 0.005, 'FaceColor', 'k', 'EdgeColor', 'none');


xlim([start_date_doy end_date_doy])

% Set common x-axis properties
% xticks(xticks_days)
datetick('x','mmm','keeplimits')

text(0.01, 0.98, subplot_labels{5}, 'Units', 'normalized', 'FontSize', 16, ...
     'FontName', ftname, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');




% Set figure size and save options
set(gcf, 'Units', 'inches', 'Position', [0, 0, 7, 8.8]);  % Double-column width with scaled height
set(gcf, 'PaperPositionMode', 'auto');  % Ensure figure fits paper size

% Saving the figure with high resolution for publication
saving = 1;
if saving == 1
    print(gcf, 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\final\figure_3_v0.png', '-dpng', '-r300');
    print(gcf, 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\final\figure_3_v0.pdf', '-dpdf', '-r300');
end

%% Abundance porporitons
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
xline(peak_fire_day, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':'); 
% xline(lag_date_day, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':'); 
xline(start_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-'); 
fill([shade_start_day, shade_start_day, shade_end_day, shade_end_day], ...
    [min(ylim) * 1.5, max(ylim) * 1.5, max(ylim) * 1.5 min(ylim) * 1.5], ...
    [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off'); 
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
    [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off'); 
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
saving=1;
if saving == 1
    saveas(gcf, 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\ifcb\morphometrics_proportions_plots.png');
    saveas(gcf, 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\ifcb\morphometrics_proportions_plots.pdf');
end


%% Biovolume and abundance growth rate by size
% Filter the dataset for the year 2020
year_filter = (ifcb_nass_morpho_2016_2023.datetime.Year == 2020);
filtered_data_plot2 = ifcb_nass_morpho_2016_2023(year_filter, :);
filtered_data_plot3 = ifcb_nass_morpho_part1_1_2016_2023(year_filter, :);

% Calculate growth rates for Plot 2 (columns 8 to 13) and Plot 3 (columns 14 to 19)
filtered_data_plot2 = calculate_growth_rate(filtered_data_plot2, 'datetime');
filtered_data_plot3 = calculate_growth_rate(filtered_data_plot3, 'datetime');

% Create a new figure for growth rate plots
figure; close all;
t = tiledlayout(2, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Plot growth rates for Plot 2 (columns 8 to 13)
nexttile;
hold on;
num_colors = size(colors, 1); % Get the number of available colors

% Iterate through columns 8 to 13 to plot growth rates
for i = 8:13
    growth_rate_column = [filtered_data_plot2.Properties.VariableNames{i}, '_growthrate'];
    if ismember(growth_rate_column, filtered_data_plot2.Properties.VariableNames)
        plot(filtered_data_plot2.dayofyear, ...
             movmean(filtered_data_plot2{:, growth_rate_column}, 7, 'omitnan'), ...
             'LineWidth', 2, 'Color', colors(mod(i-8, num_colors) + 1, :), 'DisplayName', binLabels{i-7});
    end
end

% Formatting for Plot 2 Growth Rates
ylabel(["Growth Rate" + newline + "Biovolume (ln(A_{t+1}/A_t))"], 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold');
xlim([start_date_doy end_date_doy]);
legend('show', 'Location', 'bestoutside', 'FontSize', ftsz, 'AutoUpdate', 'off');
title(legend, 'Size class [um]');
yline(0, 'k', 'LineWidth', 1.5, 'LineStyle', '--'); % Adds a horizontal line at y = 0, in black ('k')

xline(peak_fire_day, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');
xline(lag_date_day, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':');
% xline(august_date_day, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');
xline(start_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
fill([shade_start_day, shade_start_day, shade_end_day, shade_end_day], ...
     [min(ylim) * 1.5, max(ylim) * 1.5, max(ylim) * 1.5 min(ylim) * 1.5], ...
     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
xline(end_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');

% Plot growth rates for Plot 3 (columns 14 to 19)
nexttile;
hold on;

% Iterate through columns 14 to 19 to plot growth rates
for i = 14:19
    growth_rate_column = [filtered_data_plot3.Properties.VariableNames{i}, '_growthrate'];
    if ismember(growth_rate_column, filtered_data_plot3.Properties.VariableNames)
        plot(filtered_data_plot3.dayofyear, ...
             movmean(filtered_data_plot3{:, growth_rate_column}, 7, 'omitnan'), ...
             'LineWidth', 2, 'Color', colors(mod(i-14, num_colors) + 1, :), 'DisplayName', binLabels{i-13});
    end
end

% Formatting for Plot 3 Growth Rates
ylabel(["Growth Rate" + newline + "Abundance (ln(A_{t+1}/A_t))"], 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold');
xlim([start_date_doy end_date_doy]);
legend('show', 'Location', 'bestoutside', 'FontSize', ftsz, 'AutoUpdate', 'off');
title(legend, 'Size class [um]');

yline(0, 'k', 'LineWidth', 1.5, 'LineStyle', '--'); % Adds a horizontal line at y = 0, in black ('k')

xline(peak_fire_day, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');
xline(lag_date_day, 'Color', '#cf572b', 'LineWidth', 2, 'LineStyle', ':');
% xline(august_date_day, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');
xline(start_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');
fill([shade_start_day, shade_start_day, shade_end_day, shade_end_day], ...
     [min(ylim) * 1.5, max(ylim) * 1.5, max(ylim) * 1.5 min(ylim) * 1.5], ...
     [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
xline(end_date_day, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-');

% Set common x-axis properties for both plots
datetick('x', 'mmm', 'keeplimits');
xlabel('Date', 'FontSize', ftsz, 'FontName', ftname, 'FontWeight', 'bold');



dual_monitor=0;
if dual_monitor == 1
    set(gcf, 'Position', [0 1100 1400 1000])
else
    set(gcf, 'Position', [0 100 1400 1000])
end








%% Scrap

%% Plot climatology

ifcb_nass_morpho_nass_climatology=climatology_se(ifcb_nass_morpho_2016_2023,'MeanSlopeAbundance');
% plotClimatology(ifcb_nass_morpho_2016_2023(ifcb_nass_morpho_2016_2023.datetime.Year ~= 2020,:),'MeanSlopeAbundance')
close all
figure
yyaxis right
plotClimatology(ifcb_nass_morpho_2016_2023(ifcb_nass_morpho_2016_2023.datetime.Year ~= 2020,:),'MeanSlopeAbundance')
hold on
yyaxis left
plot(movmean(ifcb_nass_morpho_nass_climatology.dayofyear,1),movmean(ifcb_nass_morpho_nass_climatology.Mean_MeanSlopeAbundance,10), ...
    'LineWidth',4,'Color','r')
datetick('x')
hold on

hold on

years_sel=2020:2023;
for i=1:length(years_sel)
plot(ifcb_nass_morpho_2016_2023.dayofyear(ifcb_nass_morpho_2016_2023.datetime.Year == years_sel(i),:), ...
    ifcb_nass_morpho_2016_2023.MeanSlopeAbundance(ifcb_nass_morpho_2016_2023.datetime.Year == years_sel(i),:),'LineWidth',2, ...
    'Color',colors(i,:),'Marker','none','LineStyle','-')
hold on
end
legend({"Climatology",string(years_sel)})


%% Functions

function taxa_table = calculate_growth_rate(taxa_table, datetime_column)
    % Function to calculate growth rate for each taxa in the table
    % taxa_table: input table with taxa abundances
    % datetime_column: the name of the datetime column in the table
    
    % Sort the table by the datetime column (to ensure chronological order)
    taxa_table = sortrows(taxa_table, datetime_column);
    
    % Get the list of taxa columns (exclude the datetime column)
    taxa_names = taxa_table.Properties.VariableNames;
    taxa_names(strcmp(taxa_names, datetime_column)) = [];
    
    % Extract datetime information as a datetime array
    datetimes = table2array(taxa_table(:, datetime_column));
    
    % Iterate over each taxa column
    for i = 1:length(taxa_names)
        taxa_name = taxa_names{i};
        
        % Extract the taxa values as an array
        taxa_values = table2array(taxa_table(:, taxa_name));
        
        % Initialize growth rate with NaN values
        growth_rates = NaN(size(taxa_values));
        
        % Calculate growth rate for each consecutive time step
        for t = 2:length(taxa_values)
            if days(datetimes(t) - datetimes(t-1)) == 1  % Ensure time difference is exactly 1 day
                if taxa_values(t-1) > 0  % To avoid division by zero or negative rates
                    growth_rates(t) = log(taxa_values(t) / taxa_values(t-1));
                end
            end
        end
        
        % Append growth rate column with the suffix '_growthrate'
        taxa_table.([taxa_name, '_growthrate']) = growth_rates;
    end
end
