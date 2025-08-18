%Compile physical data

%% Load the data and add climatologies

addpath("Q:\Dante\data\MB_Wildfire_Obs\pm2_5\")
addpath("Q:\Dante\data\MB_Wildfire_Obs\other_forcings\")
addpath("Q:\Dante\data\MB_Wildfire_Obs\aeronet\")
addpath("Q:\Dante\data\MB_Wildfire_Obs\shore_stations\")
addpath("Q:\Dante\data\MB_Wildfire_Obs\processed_data\")

% Load datasets
load('habs_scwharf.mat')
load('mlm_sst.mat')
load('sanlorenzoutflow.mat')
load("BEUTIdaily_SC.mat")
load('mei.mat')
load('npgo.mat')
load('ndbc_46042.mat')
load("sc_pm2.5_daily_all.mat")
load('monterey_lvl_1.5.mat')
load('ifcb_bray.mat')
load('mlml_par_tt.mat')
load('sc_wharf_water.mat')
load('sc_wharf_weather.mat')
load('mlml_station_daily.mat', 'mlml_tt')
load('EOF2shadow_d5_daily_interpolated.mat')

% Compute anomaly
Monterey_lvl_1_5.dayOfYear = day(Monterey_lvl_1_5.datetime, 'dayofyear');


%% Select specific columns from PM2.5 dataset
sc_pm2_5_daily = sc_pm2_5_daily(:, {'pm2_5', 'LocalSiteName', 'datetime'});
sc_pm2_5_daily_slv = sc_pm2_5_daily(strcmp(sc_pm2_5_daily.LocalSiteName, 'San Lorenzo Valley Middle School'), :);
sc_pm2_5_daily_sc = sc_pm2_5_daily(strcmp(sc_pm2_5_daily.LocalSiteName, 'Santa Cruz'), :);

%% Convert tables to timetables
ifcb_bray_tt = table2timetable(ifcbbray,'RowTimes','datetime');
mei_tt = table2timetable(mei);
beuti_tt = table2timetable(BEUTI_SC);
slr_outflow_tt = table2timetable(sanlorenzodailyoutflow);
habs_tt = table2timetable(HABs_SantaCruzWharf);
sc_pm2_5_daily_sc_tt = table2timetable(sc_pm2_5_daily_sc);
sc_pm2_5_daily_slv_tt = table2timetable(sc_pm2_5_daily_slv);
ndbc46042_tt = table2timetable(ndbc46042);
mlm_sst_tt = table2timetable(mlm_sst);

%% Convert timezones to UTC
mlm_sst_tt.datetime.TimeZone = 'UTC';
ifcb_bray_tt.datetime.TimeZone = 'UTC';
mei_tt.datetime.TimeZone = 'UTC';
beuti_tt.datetime.TimeZone = 'UTC';
ndbc46042_tt.datetime.TimeZone = 'UTC';
slr_outflow_tt.datetime.TimeZone = 'UTC';
habs_tt.datetime.TimeZone = 'UTC';
sc_pm2_5_daily_sc_tt.datetime.TimeZone = 'UTC';
sc_pm2_5_daily_slv_tt.datetime.TimeZone = 'UTC';
mlmlpar_tt.datetime.TimeZone='UTC';
sc_wharf_weather_tt.datetime.TimeZone='UTC';
sc_wharf_water_tt.datetime.TimeZone='UTC';
mlml_tt.datetime.TimeZone='UTC';


%% Resample to daily values for only numeric columns
% Helper function to get numeric columns and apply retime
retime_numeric_only = @(tt) retime(tt(:, varfun(@isnumeric, tt, 'OutputFormat', 'uniform')), 'daily', 'mean');
retime_nondaily = @(tt) retime(tt(:, varfun(@isnumeric, tt, 'OutputFormat', 'uniform')), 'daily', 'nearest');

mlm_sst_tt = retime_numeric_only(mlm_sst_tt);
ifcb_bray_tt = retime_numeric_only(ifcb_bray_tt);
mei_tt = retime_nondaily(mei_tt);
beuti_tt = retime_numeric_only(beuti_tt);
ndbc46042_tt = retime_numeric_only(ndbc46042_tt);
slr_outflow_tt = retime_numeric_only(slr_outflow_tt);
habs_tt = retime_nondaily(habs_tt);
sc_pm2_5_daily_sc_tt = retime_numeric_only(sc_pm2_5_daily_sc_tt);
sc_pm2_5_daily_slv_tt = retime_numeric_only(sc_pm2_5_daily_slv_tt);
mlmlpar_tt=retime_numeric_only(mlmlpar_tt);
sc_wharf_weather_tt=retime_numeric_only(sc_wharf_weather_tt);
sc_wharf_water_tt=retime_numeric_only(sc_wharf_water_tt);
mlml_tt=retime_numeric_only(mlml_tt);

%% Subset to only desired columns
slr_outflow_tt=slr_outflow_tt(:,{'outflow'});
habs_tt=habs_tt(:,{'Chl_Volume_Filtered','Avg_Chloro','Phosphate',...
    'Nitrate','Nitrite','Nitrite_Nitrate','Silicate','Ammonium','pDA','Alexandrium_spp','Dinophysis_spp','Pseudo_nitzschia_seriata_group'});
ndbc46042_tt=ndbc46042_tt(:,{'wspd_across','wspd_along','wtmp','atmp','wvht'});
sc_wharf_water_tt=sc_wharf_water_tt(:,{'chlorophyll','oxygen','oxy_sat','temperature','turbidity'});
sc_wharf_weather_tt=sc_wharf_weather_tt(:,{'air_temperature_cm_time__mean_over_pt1h','relative_humidity_cm_time__mean_over_pt1h',...
    'solar_irradiance_cm_time__mean_over_pt1h','surface_downwelling_photosynthetic_radiative_flux_in_air','wind_speed_cm_time__mean_over_pt1h',...
    'wind_from_direction_cm_time__mean_over_pt1h'});

%% Merge all timetables into one
% Use outerjoin to preserve all dates across all datasets
merged_data = outerjoin(mlm_sst_tt, ifcb_bray_tt, 'Keys', 'datetime', 'MergeKeys', true);
merged_data = outerjoin(merged_data, mei_tt, 'Keys', 'datetime', 'MergeKeys', true);
merged_data = outerjoin(merged_data, beuti_tt, 'Keys', 'datetime', 'MergeKeys', true);
merged_data = outerjoin(merged_data, ndbc46042_tt, 'Keys', 'datetime', 'MergeKeys', true);
merged_data = outerjoin(merged_data, slr_outflow_tt, 'Keys', 'datetime', 'MergeKeys', true);
merged_data = outerjoin(merged_data, habs_tt, 'Keys', 'datetime', 'MergeKeys', true);
merged_data = outerjoin(merged_data, sc_pm2_5_daily_sc_tt, 'Keys', 'datetime', 'MergeKeys', true);
merged_data = outerjoin(merged_data, sc_pm2_5_daily_slv_tt, 'Keys', 'datetime', 'MergeKeys', true);
merged_data = outerjoin(merged_data, mlmlpar_tt, 'Keys', 'datetime', 'MergeKeys', true);
merged_data = outerjoin(merged_data, sc_wharf_weather_tt, 'Keys', 'datetime', 'MergeKeys', true);
merged_data = outerjoin(merged_data,sc_wharf_water_tt, 'Keys', 'datetime', 'MergeKeys', true);
merged_data = outerjoin(merged_data, EOF2shadow_d5_daily_interpolated , 'Keys', 'datetime', 'MergeKeys', true);


% merged_data=merged_data(merged_data.datetime > datetime(2016,1,1,'TimeZone','UTC'),:);
physical_forcings_biology_table=merged_data;

save('Q:\Dante\data\MB_Wildfire_Obs\processed_data\joined_physical_drivers\physical_forcings_biology_table_alltime.mat','physical_forcings_biology_table')
writetimetable(physical_forcings_biology_table,'Q:\Dante\data\MB_Wildfire_Obs\processed_data\joined_physical_drivers\physical_forcings_biology_table_alltime.csv')