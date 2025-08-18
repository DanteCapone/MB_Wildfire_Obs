%Process Data from Shore Station MLML

% load('Q:\Dante\data\MB_Wildfire_Obs\shore_stations\mlm_raw_5_2024.mat')
load('C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\Santa_Cruz\data\mlml_mlml_sea_d652_b4dc_f59f.mat')

%%
field_names = fieldnames(mlml_mlml_sea);
unique_field_names = matlab.lang.makeUniqueStrings(field_names);

% Create a new structure with unique field names
mlml_mlml_sea_unique = struct();
for k = 1:numel(field_names)
    mlml_mlml_sea_unique.(unique_field_names{k}) = mlml_mlml_sea.(field_names{k});
end

% Convert the unique structure to table and timetable
mlml_mlml_sea_t = struct2table(mlml_mlml_sea_unique);
mlml_mlml_sea_t.datetime=datetime(mlml_mlml_sea_t.time, 'ConvertFrom', 'posixtime', 'Format', 'MM/dd/yyyy');
mlml_tt = table2timetable(mlml_mlml_sea_t);


%% Add flur anomaly
mlml_chl=mlml_tt(mlml_tt.fluorescence_qc_agg==1,{'fluorescence'});
mlml_chl=retime(mlml_chl,'daily','median');
mlml_chl.dayOfYear=day(mlml_chl.datetime, 'dayofyear');
mlml_chl_climatology=climatology(mlml_chl,'fluorescence');
mlml_chl = join(mlml_chl, mlml_chl_climatology, 'Keys', 'dayOfYear');
mlml_chl.anomaly=mlml_chl.fluorescence-mlml_chl.nanmedian_fluorescence;


mlml_chl=mlml_tt(mlml_tt.fluorescence_qc_agg==1,{'fluorescence'});
mlml_chl=retime(mlml_chl,'daily','median');
mlml_chl.dayOfYear=day(mlml_chl.datetime, 'dayofyear');
mlml_chl_climatology=climatologyWeekly(mlml_chl,'fluorescence');
mlml_chl = join(mlml_chl, mlml_chl_climatology, 'Keys', 'dayOfYear');
mlml_chl.anomaly=mlml_chl.fluorescence-mlml_chl.nanmedian_fluorescence;


%% Plot

month_labs = {'Jan', 'Feb', 'Ma', 'Apr', 'May', 'Jun', ...
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};

% Define the date range
start_date = datetime(2010, 1, 1);
end_date = datetime(2023, 12, 31);
start_date_num = datenum(datetime(2020, 1, 1));
end_date_num = datenum(datetime(2020, 12, 31));

shade_start = day(datetime(2020,8,16),'dayofyear');
shade_end = day(datetime(2020,9,23),'dayofyear');

clf
figure(1)
subplot(2,1,1)
% plot(mlml_tt.datetime(mlml_tt.fluorescence_qc_agg==1),movmean(mlml_tt.fluorescence(mlml_tt.fluorescence_qc_agg==1),12*24,1,"omitmissing"))
% hold on
% plot(mlml_chl.datetime,mlml_chl.fluorescence)
% plot(mlml_chl.dayOfYear(1:365),mlml_chl.nanmedian_fluorescence(1:365))
% hold on
% plot(mlml_chl.dayOfYear(mlml_chl.datetime.Year==2020),mlml_chl.anomaly(mlml_chl.datetime.Year==2020))
yyaxis left
shade_anomaly(HABs_SantaCruzWharf.dayOfYear(HABs_SantaCruzWharf.datetime.Year==2020),HABs_SantaCruzWharf.anomaly_Avg_Chloro(HABs_SantaCruzWharf.datetime.Year==2020),'#cc78c9','#9a45f5')
ylim([-5 90])
hold on
ylabel(['Moss Landing Marine Lab',newline,' 2020 Chlorophyll-a Anomaly',newline,'[mg/m^3]'],'FontSize',14)

yyaxis right
shade_anomaly(mlml_chl.dayOfYear(mlml_chl.datetime.Year==2020),mlml_chl.anomaly(mlml_chl.datetime.Year==2020))
ylim([-5 90])

set(gca, 'XTick', 1:30:365, 'XTickLabel', month_labs, 'XTickLabelRotation', 45);  % Use the same x-axis labels
xlim([0 366])
datetick('x','keeplimits')
hold on
ylabel(['CalHABMap Weekly Chlorophyll-a Anomaly',newline,'[mg/m^3]'])

% fill([shade_start, shade_start, shade_end, shade_end], ...
%      [min(ylim) * 1, max(ylim) * 1.1, max(ylim) * 1.1, min(ylim) * 1], ...
%      [0.6509 0.8078 0.8902], 'EdgeColor', 'none', 'FaceAlpha', 0.2); % Adjust 'FaceAlpha' for transparency

% grid off
set(gcf,'Position',[0 0 1400 800])
saving=1;
if saving==1
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\shore_stations\mlml_and_calhabmap_chl.png"]);
    saveas(gcf, ["C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs\MB_Wildfire_Obs\figures\shore_stations\mlml_and_calhabmap_chl.pdf"]);
end

%% Plot climatology and anomaly

plotClimatology(mlml_chl(mlml_chl.datetime.Year ~= 2020,:),'fluorescence')




%%

joined_tt=synchronize(ifcb_bray_tt,mlml_tt,'first','nearest');


%%
start_date=datetime(2010,1,1,'TimeZone','America/Los_Angeles');
end_date=datetime(2024,12,31,'TimeZone','America/Los_Angeles');

joined_tt=joined_tt(joined_tt.datetime>=start_date & joined_tt.datetime <= end_date,:);
%QC
joined_tt.fluorescence(joined_tt.fluorescence_qc_agg~=1)=NaN;
joined_tt.mole_concentration_of_nitrate_i(joined_tt.mole_concentration_of_nitrate_i_1~=1)=NaN;

joined_tt.dayOfYear=day(joined_tt.datetime,'dayofyear');

% Extract the taxa columns and AOD_500nm column
taxa_data = table2array(joined_tt(:, 3:53));
sst=joined_tt.sea_water_temperature;
fluorescence=joined_tt.fluorescence;
nitrate=joined_tt.mole_concentration_of_nitrate_i;


% Initialize variables to store coefficients and lags
cc= zeros(size(taxa_data, 1)*2-1,size(taxa_data, 2));
lags = zeros(size(taxa_data, 2), size(taxa_data, 1)*2-1);

taxa=joined_tt.Properties.VariableNames(3:53);
taxa=significant_taxa;
for i=1:size(taxa_data,1)

    %Compute climatology and anomaly
    eval([taxa{i},'_climatology=climatology(joined_tt,taxa{i});']);
    eval(['joined_tt = join(joined_tt,',taxa{i},'_climatology, ''Keys'', ''dayOfYear'');']);
    eval(['joined_tt.',taxa{i},'_anomaly=joined_tt.',taxa{i},'-joined_tt.nanmedian_',taxa{i},'_joined_tt;'])

    %Cross-correlate
    eval(['[cc(:,i),lags(i,:)] = xcorr(fillmissing(nitrate, ''linear''),joined_tt.',taxa{i},'_anomaly,''coeff'');'])
    % Plot the cross-correlation
    plotting=1;
    if plotting==1
    figure(1)
    % stem(lags, cc);
    plot(lags(i,:),cc(:,i),'LineStyle','-','LineWidth',2)
    hold on
    title("Cross-correlation between AOD 500nm &"+newline+ "Phytoplankton Groups");
    xlabel('Lags (Days)');
    ylabel('Normalized Cross-correlation');
    end
    xlim([-110,110])
end
 
