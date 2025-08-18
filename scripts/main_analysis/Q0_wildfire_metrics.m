%Wildfire Metric Anlaysis

% Add paths
addpath(genpath('Q:\Dante\Wildfire_Obs'));
addpath(genpath('Q:\Dante\data\MB_Wildfire_Obs\processed_data'));
addpath('Q:\Dante\data\MB_Wildfire_Obs\aeronet')

%% Load the data
%Select fire  years to exlcude
fire_years=[2008, 2016,2020,2021];

%Aeronet
load('Q:\Dante\data\MB_Wildfire_Obs\aeronet\monterey_lvl_1.5.mat')

% %Compute and add anomaly
% Monterey_lvl_1_5.dayOfYear=day(Monterey_lvl_1_5.datetime, 'dayofyear');
% Monterey_lvl_1_5_climatology=climatology(Monterey_lvl_1_5,'AOD_500nm',1);
% Monterey_lvl_1_5 = join(Monterey_lvl_1_5, Monterey_lvl_1_5_climatology, 'Keys', 'dayOfYear');
% Monterey_lvl_1_5.anomaly=Monterey_lvl_1_5.AOD_500nm-Monterey_lvl_1_5.AOD_500nm_climatology;
% 
% save('Q:\Dante\data\MB_Wildfire_Obs\aeronet\monterey_lvl_1.5.mat','Monterey_lvl_1_5')

% Load PM2.5 dataset
load('Q:\Dante\data\MB_Wildfire_Obs\pm2_5\sc_pm2.5_daily_all.mat')

% Compute and add anomaly
% sc_pm2_5_daily.dayOfYear = day(sc_pm2_5_daily.datetime, 'dayofyear');
% sc_pm2_5_daily_climatology = climatology(sc_pm2_5_daily, 'pm2_5');
% sc_pm2_5_daily = join(sc_pm2_5_daily, sc_pm2_5_daily_climatology, 'Keys', 'dayOfYear');
% sc_pm2_5_daily.anomaly = sc_pm2_5_daily.pm2_5 - sc_pm2_5_daily.pm2_5_climatology;
% 
% save('Q:\Dante\data\MB_Wildfire_Obs\pm2_5\sc_pm2.5_daily_all.mat','sc_pm2_5_daily')
%Subset PM 2.5 by site
% Santa Cruz
sc_pm2_5_daily_sc=sc_pm2_5_daily(sc_pm2_5_daily.LocalSiteName=="Santa Cruz",:);
% % SLV Middle
sc_pm2_5_daily_slv=sc_pm2_5_daily(sc_pm2_5_daily.LocalSiteName=="San Lorenzo Valley Middle School",:);



% MERRA-2 averaged for the monterey Bay
load('merra_MB_bc_avg_daily_plt.mat')


%% AOD

clf
colorz = [
    0.75, 0.1, 0;        % Darker shade of Red
    0.3333, 0.6588, 0.4078;  % Green
    0.4, 0.6, 0.8;        % Blue
    0.7882, 0.4588, 0.9216   % Pastel Purple
];

close all;
yyaxis left;
plotClimatology5day(Monterey_lvl_1_5(~ismember(Monterey_lvl_1_5.datetime.Year, fire_years),:),'AOD_500nm')
alpha(0.1)
ylim([0 0.2])
hold on 
ax = gca;  % Get current axes
ax.FontSize = ftsz;


yyaxis right;
% plot(Monterey_lvl_1_5.Day_of_Year(Monterey_lvl_1_5.datetime.Year == 2008),Monterey_lvl_1_5.AOD_500nm(Monterey_lvl_1_5.datetime.Year == 2008),'Color',colorz(1,:),'LineStyle','-','LineWidth',2);
hold on
% plot(Monterey_lvl_1_5.Day_of_Year(Monterey_lvl_1_5.datetime.Year == 2016),Monterey_lvl_1_5.AOD_500nm(Monterey_lvl_1_5.datetime.Year == 2016),'Color',colorz(2,:),'LineStyle','-','LineWidth',2);
hold on
plot(Monterey_lvl_1_5.Day_of_Year(Monterey_lvl_1_5.datetime.Year == 2020),Monterey_lvl_1_5.AOD_500nm(Monterey_lvl_1_5.datetime.Year == 2020),'Color',colorz(1,:),'LineStyle','-','LineWidth',2.5);
hold on
% plot(Monterey_lvl_1_5.Day_of_Year(Monterey_lvl_1_5.datetime.Year == 2021),Monterey_lvl_1_5.AOD_500nm(Monterey_lvl_1_5.datetime.Year == 2021),'Color',colorz(4,:),'LineStyle','-','LineWidth',2)
ylabel(['AOD 500nm'],'FontSize',ftsz);
% title('Monterey Aeronet Aerosol Optical Depth');

% text
% text(148,1.5,["Indians Fire"+newline+"(Monterey, 2008)"],"Color",'r')
% text(175,3.8,["Soberanes Fire"+newline+"(Monterey, 2016)"],"Color",'g')
text(240, 300, ["CZU/SCU Lightning"+ newline+"Complex Fire" + newline + "(Santa Cruz, 2020)"], "Color", 'k')
text(260, 80, ["August Complex Fire" + newline + "(Santa Cruz, 2020)"], "Color", 'k')
% text(212,1.6,["Dixie Fire"+newline+"(2021)"],"Color",[0.7882, 0.4588, 0.9216])
legend({'Climatology (2003-present)','2008','2016','2020','2021'})
legend({'Climatology (2003-present)','2020'})
ax = gca;  % Get current axes
ax.FontSize = ftsz;

set(gcf,'Position',[0 100 1600 600])

saving=1;
if saving==1
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\wildfire_metrica\aeronet\monterey_aod_plot.png"]);
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\wildfire_metrica\aeronet\monterey_aod_plot.pdf"]);
end
saving=0;


%% PM 2.5 Analysis




% Compute and add anomaly
sc_pm2_5_daily.dayOfYear = day(sc_pm2_5_daily.datetime, 'dayofyear');
sc_pm2_5_daily_climatology = climatology(sc_pm2_5_daily, 'pm2_5');
sc_pm2_5_daily = join(sc_pm2_5_daily, sc_pm2_5_daily_climatology, 'Keys', 'dayOfYear');
sc_pm2_5_daily.anomaly = sc_pm2_5_daily.pm2_5 - sc_pm2_5_daily.Fun_pm2_5;

%%
clf
plot(sc_pm2_5_daily.datetime,sc_pm2_5_daily.pm2_5)
% %% PM2.5 Plot: Santa Cruz subplots
% 
num_subplots=2;
figure(1)
clf; close all


figure(1)
subplot(num_subplots,1,1)
sc_pm2_5_daily_sc=sc_pm2_5_daily(sc_pm2_5_daily.LocalSiteName=="Santa Cruz",:);
yyaxis left;
plotClimatology5day(sc_pm2_5_daily_sc(sc_pm2_5_daily_sc.datetime.Year ~= 2020,:),'pm2_5')
alpha(0.5)
% ylim([0 100])
hold on 

% ylabel(['PM 2.5 Climatology (\mu g/m^3)']);

% yyaxis right;
% plot(sc_pm2_5_daily_sc.dayOfYear(sc_pm2_5_daily_sc.datetime.Year == 2008), sc_pm2_5_daily_sc.pm2_5(sc_pm2_5_daily_sc.datetime.Year == 2008), 'Color', colorz(1,:), 'LineStyle', '-', 'LineWidth', 2);
hold on
% plot(sc_pm2_5_daily_sc.dayOfYear(sc_pm2_5_daily_sc.datetime.Year == 2016), sc_pm2_5_daily_sc.pm2_5(sc_pm2_5_daily_sc.datetime.Year == 2016), 'Color', colorz(2,:), 'LineStyle', '-', 'LineWidth', 2);
hold on
plot(sc_pm2_5_daily_sc.dayOfYear(sc_pm2_5_daily_sc.datetime.Year == 2020), sc_pm2_5_daily_sc.pm2_5(sc_pm2_5_daily_sc.datetime.Year == 2020), 'Color', colorz(3,:), 'LineStyle', '-', 'LineWidth', 2);
hold on
% plot(sc_pm2_5_daily_sc.dayOfYear(sc_pm2_5_daily_sc.datetime.Year == 2021), sc_pm2_5_daily_sc.pm2_5(sc_pm2_5_daily_sc.datetime.Year == 2021), 'Color', colorz(4,:), 'LineStyle', '-', 'LineWidth', 2);

% xlim([0 366])
ylim([0 100])

ylabel(['Daily PM 2.5 (\mu g/m^3)']);
title('Santa Cruz PM2.5 Concentration: Santa Cruz');

% Text annotations for significant events
% text(148, 35, ["Indians Fire" + newline + "(Monterey, 2008)"], "Color", colorz(1,:))
% text(175, 45, ["Soberanes Fire" + newline + "(Monterey, 2016)"], "Color", colorz(2,:))
text(195, 80, ["CZU/SCU Lightning"+ newline+"Complex Fire" + newline + "(Santa Cruz, 2020)"], "Color", colorz(3,:))
text(212,50,["Dixie Fire"+newline+"(2021)"],"Color",colorz(4,:))

legend({'Climatology (2008-present)','2020','2021'})



subplot(num_subplots,1,2)
% figure(2)
yyaxis left;
plotClimatology5day(sc_pm2_5_daily_slv(sc_pm2_5_daily_slv.datetime.Year ~= 2020,:),'pm2_5')
alpha(0.5)
ylim([0 30])
hold on 

ylabel(['PM 2.5 Climatology (\mu g/m^3)']);

yyaxis right;
% plot(sc_pm2_5_daily_slv.dayOfYear(sc_pm2_5_daily_slv.datetime.Year == 2008), sc_pm2_5_daily_slv.pm2_5(sc_pm2_5_daily_slv.datetime.Year == 2008), 'Color', colorz(1,:), 'LineStyle', '-', 'LineWidth', 2);
hold on
% plot(sc_pm2_5_daily_slv.dayOfYear(sc_pm2_5_daily_slv.datetime.Year == 2016), sc_pm2_5_daily_slv.pm2_5(sc_pm2_5_daily_slv.datetime.Year == 2016), 'Color', colorz(2,:), 'LineStyle', '-', 'LineWidth', 2);
% hold on
plot(sc_pm2_5_daily_slv.dayOfYear(sc_pm2_5_daily_slv.datetime.Year == 2020), sc_pm2_5_daily_slv.pm2_5(sc_pm2_5_daily_slv.datetime.Year == 2020), 'Color', colorz(1,:), 'LineStyle', '-', 'LineWidth', 2);
hold on
% plot(sc_pm2_5_daily_slv.dayOfYear(sc_pm2_5_daily_slv.datetime.Year == 2021), sc_pm2_5_daily_slv.pm2_5(sc_pm2_5_daily_slv.datetime.Year == 2021), 'Color', colorz(4,:), 'LineStyle', '-', 'LineWidth', 2);

xlim([0 366])

ylim([0 400])

ylabel(['Daily PM 2.5 (\mu g/m^3)']);
title('Santa Cruz PM2.5 Concentration: San Lorenzo Valley Middle');

% Text annotations for significant events
% text(148, 35, ["Indians Fire" + newline + "(Monterey, 2008)"], "Color", 'g')
% text(175, 45, ["Soberanes Fire" + newline + "(Monterey, 2016)"], "Color", colorz(2,:))
text(240, 300, ["CZU/SCU Lightning"+ newline+"Complex Fire" + newline + "(Santa Cruz, 2020)"], "Color", 'k')
text(260, 80, ["August Complex Fire" + newline + "(Santa Cruz, 2020)"], "Color", 'k')
% text(212,120,["Dixie Fire"+newline+"(2021)"],"Color",colorz(4,:))

legend({'Climatology (2016-present)','2020','2021'})

% 
% set(gcf,'Position',[0 100 1600 600])
% 
% saving=0;
% if saving==1
%     saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\wildfire_metrica\pm2.5\sc_pm25_plot.png"]);
%     saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\wildfire_metrica\pm2.5\sc_pm25_plot.pdf"]);
% end
% saving=0;

%% Santa Cruz and SLV Pm 2.5
figure(1)
clf; close all


sc_pm2_5_daily_sc=sc_pm2_5_daily(sc_pm2_5_daily.LocalSiteName=="Santa Cruz",:);
yyaxis left;
plotClimatology5day(sc_pm2_5_daily(sc_pm2_5_daily.datetime.Year ~= 2020,:),'pm2_5')

alpha(0.5)
ylim([0 16])
hold on 

ylabel(['PM 2.5 Climatology (\mu g/m^3)']);
ax = gca;  % Get current axes
ax.FontSize = ftsz;


% Plot fire years
yyaxis right;
% plot(sc_pm2_5_daily_sc.dayOfYear(sc_pm2_5_daily_sc.datetime.Year == 2008), sc_pm2_5_daily_sc.pm2_5(sc_pm2_5_daily_sc.datetime.Year == 2008), 'Color', colorz(1,:), 'LineStyle', '-', 'LineWidth', 2);
hold on
% plot(sc_pm2_5_daily_sc.dayOfYear(sc_pm2_5_daily_sc.datetime.Year == 2016), sc_pm2_5_daily_sc.pm2_5(sc_pm2_5_daily_sc.datetime.Year == 2016), 'Color', colorz(2,:), 'LineStyle', '-', 'LineWidth', 2);
hold on
plot(sc_pm2_5_daily_sc.dayOfYear(sc_pm2_5_daily_sc.datetime.Year == 2020), sc_pm2_5_daily_sc.pm2_5(sc_pm2_5_daily_sc.datetime.Year == 2020), 'Color', colorz(2,:), 'LineStyle', '-', 'LineWidth', 2);
hold on
plot(sc_pm2_5_daily_slv.da(sc_pm2_5_daily_slv.datetime.Year == 2020), sc_pm2_5_daily_slv.pm2_5(sc_pm2_5_daily_slv.datetime.Year == 2020), 'Color', colorz(3,:), 'LineStyle', '-', 'LineWidth', 2);

xlim([0 366])
ylim([0 400])

ylabel(['Daily PM 2.5 (\mu g/m^3)']);
% title('Santa Cruz PM2.5 Concentration: Santa Cruz');
ax = gca;  % Get current axes
ax.FontSize = ftsz;

% Text annotations for significant events
% text(148, 35, ["Indians Fire" + newline + "(Monterey, 2008)"], "Color", colorz(1,:))
% text(175, 45, ["Soberanes Fire" + newline + "(Monterey, 2016)"], "Color", colorz(2,:))
% text(240, 300, ["CZU/SCU Lightning"+ newline+"Complex Fire" + newline + "(Santa Cruz, 2020)"], "Color", 'k')
% text(260, 80, ["August Complex Fire" + newline + "(Santa Cruz, 2020)"], "Color", 'k')

% text(212,50,["Dixie Fire"+newline+"(2021)"],"Color",colorz(4,:))

legend({'Climatology (2008-present)','2020','2021'})


% SLV Middle
sc_pm2_5_daily_slv=sc_pm2_5_daily(sc_pm2_5_daily.LocalSiteName=="San Lorenzo Valley Middle School",:);


legend({'Climatology (2016-present)','Santa Cruz','San Lorenzo Valley Middle School'})


set(gcf,'Position',[0 100 1600 600])

saving=1;
if saving==1
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\wildfire_metrica\pm2.5\sc_pm25_plot.png"]);
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\wildfire_metrica\pm2.5\sc_pm25_plot.pdf"]);
end

%


%% BC Plot

clf; close all

figure()
yyaxis left;
plotClimatology5day(merra_bc_plt(~ismember(merra_bc_plt.datetime.Year, fire_years),:),'Median_BCCMASS')
ylabel(['Average Black Carbon Concentration [µg/m²]']);

alpha(0.2)
ylim([3e-7 12e-7])
hold on 
ax = gca;  % Get current axes
ax.FontSize = ftsz;


yyaxis right;
% plot(merra_bc_plt.dayOfYear(merra_bc_plt.datetime.Year == 2008), merra_bc_plt.Median_BCCMASS(merra_bc_plt.datetime.Year == 2008), 'Color', colorz(1,:), 'LineStyle', '-', 'LineWidth', 2);
hold on
% plot(merra_bc_plt.dayOfYear(merra_bc_plt.datetime.Year == 2016), merra_bc_plt.Median_BCCMASS(merra_bc_plt.datetime.Year == 2016), 'Color', colorz(2,:), 'LineStyle', '-', 'LineWidth', 2);
hold on
plot(merra_bc_plt.dayOfYear(merra_bc_plt.datetime.Year == 2020), merra_bc_plt.Median_BCCMASS(merra_bc_plt.datetime.Year == 2020), 'Color', colorz(4,:), 'LineStyle', '-', 'LineWidth', 2);
hold on
% plot(merra_bc_plt.dayOfYear(merra_bc_plt.datetime.Year == 2021), merra_bc_plt.Median_BCCMASS(merra_bc_plt.datetime.Year == 2021), 'Color', colorz(4,:), 'LineStyle', '-', 'LineWidth', 2);

ylabel(['Black Carbon Concentration [µg/m²]'],'FontSize',ftsz);
xlabel([''],'FontSize',ftsz);
ax = gca;  % Get current axes
ax.FontSize = ftsz;

% title('Monterey Bay MERRA-2 Black Carbon Surface Mass Concentration');

% Text annotations for significant events
% text(146, 0.5e-5, ["Indians Fire" + newline + "(Monterey, 2008)"], "Color",  colorz(1,:))
% text(181, 0.7e-5, ["Soberanes Fire" + newline + "(Monterey, 2016)"], "Color",  colorz(2,:))
% text(210, 2.2e-5, ["CZU/SCU Lightning"+ newline+"Complex Fire" + newline + "(Santa Cruz, 2020)"], "Color",  colorz(3,:))
% text(200, 1.7e-5,["Dixie Fire"+newline+"(2021)"],"Color",colorz(4,:))
% text(195, 0.5e-5, ["CZU/SCU Lightning"+ newline+"Complex Fire" + newline + "(Santa Cruz, 2020)"], "Color", 'k')
% text(260, 2.5e-5, ["August Complex Fire" + newline + "(Santa Cruz, 2020)"], "Color", 'k')

legend({'Climatology (2003-present)','2008','2016','2020','2021'})
legend({'Climatology (2003-present)','2020'})

set(gcf,'Position',[0 0 1600 600])

saving=1;
if saving==1
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\wildfire_metrica\merra-2\merra2_monterey_plot.png"]);
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\wildfire_metrica\merra-2\merra2_monterey_plot.pdf"]);
end
saving=0;


%% Plot with 3 Metrics as Subplots (Helvetica Font, Sized for Publication)

% Set publication settings
ftsz = 10;  % Font size for publication
ftname = 'Helvetica';  % Font name
linewdt = 2;  % Line width for plots
legend_on = 0;
num_subplots = 3;

close all;
figure()

% PM 2.5 Plot
subplot(num_subplots, 1, 1)
yyaxis left;
plot(sc_pm2_5_daily_slv.dayOfYear(sc_pm2_5_daily_slv.datetime.Year == 2020), ...
    sc_pm2_5_daily_slv.pm2_5(sc_pm2_5_daily_slv.datetime.Year == 2020), ...
    'Color', colorz(3,:), 'LineStyle', '-', 'LineWidth', linewdt);
hold on;
plot(sc_pm2_5_daily_sc.dayOfYear(sc_pm2_5_daily_sc.datetime.Year == 2020), ...
    sc_pm2_5_daily_sc.pm2_5(sc_pm2_5_daily_sc.datetime.Year == 2020), ...
    'Color', colorz(2,:), 'LineStyle', '-', 'LineWidth', linewdt);
ylim([0 400])
ylabel(["PM 2.5 (\mu g/m^3)"], 'FontSize', ftsz, 'FontName', ftname);

yyaxis right;
plotClimatology5day(sc_pm2_5_daily_slv(sc_pm2_5_daily_slv.datetime.Year ~= 2020, :), 'pm2_5')
alpha(0.5)
ylim([0 30])
xlim([1 365])
hold on;
ylabel(["Average PM 2.5 (\mu g/m^3)"], 'FontSize', ftsz, 'FontName', ftname);

ax = gca;
ax.FontSize = ftsz;
ax.FontName = ftname;
if legend_on == 1
    legend({'Climatology (2016-present)', '2020'}, 'FontSize', ftsz, 'FontName', ftname)
end

% Add "a)" label
text(0.02, 0.9, 'a)', 'Units', 'normalized', 'FontSize', ftsz+2, 'FontName', ftname,'FontWeight','bold');

% AOD Plot
subplot(num_subplots, 1, 2)
yyaxis left;
plot(Monterey_lvl_1_5.Day_of_Year(Monterey_lvl_1_5.datetime.Year == 2020), ...
    Monterey_lvl_1_5.AOD_500nm(Monterey_lvl_1_5.datetime.Year == 2020), ...
    'Color', colorz(1,:), 'LineStyle', '-', 'LineWidth', 2.5);
ylabel('AOD 500nm', 'FontSize', ftsz, 'FontName', ftname);
xlim([0 365])
if legend_on == 1
    legend({'Climatology (2003-present)', '2020'}, 'FontSize', ftsz, 'FontName', ftname)
end
ax = gca;
ax.FontSize = ftsz;
ax.FontName = ftname;

yyaxis right;
plotClimatology5day(Monterey_lvl_1_5(~ismember(Monterey_lvl_1_5.datetime.Year, fire_years), :), 'AOD_500nm')
alpha(0.1)
ylim([0 0.2])
xlim([1 365])
hold on;
ax = gca;
ax.FontSize = ftsz;
ax.FontName = ftname;

% Add "b)" label
text(0.02, 0.9, 'b)', 'Units', 'normalized', 'FontSize', ftsz+2, 'FontName', ftname,'FontWeight','bold');

% Black Carbon Mass Plot
subplot(num_subplots, 1, 3)
yyaxis left;
plot(merra_bc_plt.dayOfYear(merra_bc_plt.datetime.Year == 2020), ...
    merra_bc_plt.Median_BCCMASS(merra_bc_plt.datetime.Year == 2020), ...
    'Color', colorz(4,:), 'LineStyle', '-', 'LineWidth', linewdt);
ylabel(["Black Carbon" + newline + "Concentration [µg/m²]"], 'FontSize', ftsz, 'FontName', ftname);

ax = gca;
ax.FontSize = ftsz;
ax.FontName = ftname;

yyaxis right;
plotClimatology5day(merra_bc_plt(~ismember(merra_bc_plt.datetime.Year, fire_years), :), 'Median_BCCMASS')
ylabel(["Average Black Carbon" + newline + "Concentration [µg/m²]"], 'FontSize', ftsz, 'FontName', ftname);
alpha(0.2)
ylim([3e-7 12e-7])
xlim([1 365])
hold on;
ax = gca;
ax.FontSize = ftsz;
ax.FontName = ftname;

if legend_on == 1
    legend({'Climatology (2003-present)', '2020'}, 'FontSize', ftsz, 'FontName', ftname)
end

% Add "c)" label
text(0.02, 0.9, 'c)', 'Units', 'normalized', 'FontSize', ftsz+2, 'FontName', ftname,'FontWeight','bold');

% Adjust figure size for publication (double-column width)
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6.16, 8]);  % 7.16 inches width and 8 inches height for subplots
set(gcf, 'PaperPositionMode', 'auto');  % Ensure the figure fits the paper size

% Save the figure if required
saving = 1;
if saving == 1
    print(gcf, '-dpng', '-r300', 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\final\figure_1_v0.png');
    print(gcf, '-dpdf', '-r300', 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\final\figure_1_v0.pdf');
end
saving = 0;
