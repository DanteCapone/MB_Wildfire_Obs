%Wildfire Metric Anlaysis

%% Plotting parameters
ftsz=16;

%%


%Select fire  years to exlcude
fire_years=[2008, 2016,2020,2021];

%Aeronet
load('Q:\Dante\data\MB_Wildfire_Obs\aeronet\monterey_lvl_1.5.mat')
load('Q:\Dante\data\MB_Wildfire_Obs\aeronet\monterey_lvl_2.0.mat')
ind=Monterey_lvl_1_5.AOD_500nm == -999.0;
Monterey_lvl_1_5.AOD_500nm(ind)=nan;

%Compute and add anomaly
Monterey_lvl_1_5.dayOfYear=day(Monterey_lvl_1_5.datetime, 'dayofyear');
Monterey_lvl_1_5_climatology=climatology(Monterey_lvl_1_5,'AOD_500nm');
Monterey_lvl_1_5 = join(Monterey_lvl_1_5, Monterey_lvl_1_5_climatology, 'Keys', 'dayOfYear');
Monterey_lvl_1_5.anomaly=Monterey_lvl_1_5.AOD_500nm-Monterey_lvl_1_5.Fun_AOD_500nm;

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
plotClimatology(Monterey_lvl_1_5(~ismember(Monterey_lvl_1_5.datetime.Year, fire_years),:),'AOD_500nm')
alpha(0.1)
ylim([0 0.2])
hold on 
ax = gca;  % Get current axes
ax.FontSize = 16;


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
ax.FontSize = 16;

set(gcf,'Position',[0 100 1600 600])

saving=1;
if saving==1
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\wildfire_metrica\aeronet\monterey_aod_plot.png"]);
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\wildfire_metrica\aeronet\monterey_aod_plot.pdf"]);
end
saving=0;


%% PM 2.5 Analysis

% Load PM2.5 dataset
load('Q:\Dante\data\MB_Wildfire_Obs\pm2_5\sc_pm2.5_daily_all.mat')


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
% num_subplots=2;
% figure(1)
% clf; close all
% 
% 
% figure(1)
% subplot(num_subplots,1,1)
sc_pm2_5_daily_sc=sc_pm2_5_daily(sc_pm2_5_daily.LocalSiteName=="Santa Cruz",:);
yyaxis left;
plotClimatology(sc_pm2_5_daily_sc(sc_pm2_5_daily_sc.datetime.Year ~= 2020,:),'pm2_5')
alpha(0.5)
% ylim([0 100])
hold on 

% ylabel(['PM 2.5 Climatology (\mu g/m^3)']);

% yyaxis right;
plot(sc_pm2_5_daily_sc.dayOfYear(sc_pm2_5_daily_sc.datetime.Year == 2008), sc_pm2_5_daily_sc.pm2_5(sc_pm2_5_daily_sc.datetime.Year == 2008), 'Color', colorz(1,:), 'LineStyle', '-', 'LineWidth', 2);
hold on
plot(sc_pm2_5_daily_sc.dayOfYear(sc_pm2_5_daily_sc.datetime.Year == 2016), sc_pm2_5_daily_sc.pm2_5(sc_pm2_5_daily_sc.datetime.Year == 2016), 'Color', colorz(2,:), 'LineStyle', '-', 'LineWidth', 2);
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

% 
% % SLV Middle
% sc_pm2_5_daily_slv=sc_pm2_5_daily(sc_pm2_5_daily.LocalSiteName=="San Lorenzo Valley Middle School",:);
% 
% subplot(num_subplots,1,2)
% % figure(2)
% yyaxis left;
% plotClimatology(sc_pm2_5_daily_slv(sc_pm2_5_daily_slv.datetime.Year ~= 2020,:),'pm2_5')
% alpha(0.5)
% ylim([0 30])
% hold on 
% 
% ylabel(['PM 2.5 Climatology (\mu g/m^3)']);
% 
% yyaxis right;
% % plot(sc_pm2_5_daily_slv.dayOfYear(sc_pm2_5_daily_slv.datetime.Year == 2008), sc_pm2_5_daily_slv.pm2_5(sc_pm2_5_daily_slv.datetime.Year == 2008), 'Color', colorz(1,:), 'LineStyle', '-', 'LineWidth', 2);
% hold on
% % plot(sc_pm2_5_daily_slv.dayOfYear(sc_pm2_5_daily_slv.datetime.Year == 2016), sc_pm2_5_daily_slv.pm2_5(sc_pm2_5_daily_slv.datetime.Year == 2016), 'Color', colorz(2,:), 'LineStyle', '-', 'LineWidth', 2);
% % hold on
% plot(sc_pm2_5_daily_slv.dayOfYear(sc_pm2_5_daily_slv.datetime.Year == 2020), sc_pm2_5_daily_slv.pm2_5(sc_pm2_5_daily_slv.datetime.Year == 2020), 'Color', colorz(1,:), 'LineStyle', '-', 'LineWidth', 2);
% hold on
% % plot(sc_pm2_5_daily_slv.dayOfYear(sc_pm2_5_daily_slv.datetime.Year == 2021), sc_pm2_5_daily_slv.pm2_5(sc_pm2_5_daily_slv.datetime.Year == 2021), 'Color', colorz(4,:), 'LineStyle', '-', 'LineWidth', 2);
% 
% xlim([0 366])
% 
% ylim([0 400])
% 
% ylabel(['Daily PM 2.5 (\mu g/m^3)']);
% title('Santa Cruz PM2.5 Concentration: San Lorenzo Valley Middle');
% 
% % Text annotations for significant events
% % text(148, 35, ["Indians Fire" + newline + "(Monterey, 2008)"], "Color", 'g')
% % text(175, 45, ["Soberanes Fire" + newline + "(Monterey, 2016)"], "Color", colorz(2,:))
% text(240, 300, ["CZU/SCU Lightning"+ newline+"Complex Fire" + newline + "(Santa Cruz, 2020)"], "Color", 'k')
% text(260, 80, ["August Complex Fire" + newline + "(Santa Cruz, 2020)"], "Color", 'k')
% % text(212,120,["Dixie Fire"+newline+"(2021)"],"Color",colorz(4,:))
% 
% legend({'Climatology (2016-present)','2020','2021'})
% 
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
plotClimatology(sc_pm2_5_daily(sc_pm2_5_daily.datetime.Year ~= 2020,:),'pm2_5')

alpha(0.5)
ylim([0 16])
hold on 

ylabel(['PM 2.5 Climatology (\mu g/m^3)']);
ax = gca;  % Get current axes
ax.FontSize = 16;


% Plot fire years
yyaxis right;
% plot(sc_pm2_5_daily_sc.dayOfYear(sc_pm2_5_daily_sc.datetime.Year == 2008), sc_pm2_5_daily_sc.pm2_5(sc_pm2_5_daily_sc.datetime.Year == 2008), 'Color', colorz(1,:), 'LineStyle', '-', 'LineWidth', 2);
hold on
% plot(sc_pm2_5_daily_sc.dayOfYear(sc_pm2_5_daily_sc.datetime.Year == 2016), sc_pm2_5_daily_sc.pm2_5(sc_pm2_5_daily_sc.datetime.Year == 2016), 'Color', colorz(2,:), 'LineStyle', '-', 'LineWidth', 2);
hold on
plot(sc_pm2_5_daily_sc.dayOfYear(sc_pm2_5_daily_sc.datetime.Year == 2020), sc_pm2_5_daily_sc.pm2_5(sc_pm2_5_daily_sc.datetime.Year == 2020), 'Color', colorz(2,:), 'LineStyle', '-', 'LineWidth', 2);
hold on
plot(sc_pm2_5_daily_slv.dayOfYear(sc_pm2_5_daily_slv.datetime.Year == 2020), sc_pm2_5_daily_slv.pm2_5(sc_pm2_5_daily_slv.datetime.Year == 2020), 'Color', colorz(3,:), 'LineStyle', '-', 'LineWidth', 2);

xlim([0 366])
ylim([0 400])

ylabel(['Daily PM 2.5 (\mu g/m^3)']);
% title('Santa Cruz PM2.5 Concentration: Santa Cruz');
ax = gca;  % Get current axes
ax.FontSize = 16;

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

%% MERRA-2 averaged for the monterey Bay

load('merra_bc_plt.mat')


%% BC Plot

clf; close all

figure()
yyaxis left;
plotClimatology(merra_bc_plt(~ismember(merra_bc_plt.datetime.Year, fire_years),:),'Median_BCCMASS')
ylabel(['Average Black Carbon Concentration [µg/m²]']);

alpha(0.2)
ylim([3e-7 12e-7])
hold on 
ax = gca;  % Get current axes
ax.FontSize = 16;


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
ax.FontSize = 16;

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


