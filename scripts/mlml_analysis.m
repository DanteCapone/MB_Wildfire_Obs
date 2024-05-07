%Process Data from Shore Station MLML

load('Q:\Dante\data\MB_Wildfire_Obs\shore_stations\mlm_raw_5_2024.mat')

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
mlml_tt = table2timetable(mlml_mlml_sea_t);
mlm_tt_daily=retime(mlml_tt,'daily');

%% Plot

% Define the date range
start_date = datetime(2020, 7, 1);
end_date = datetime(2020, 11, 1);

mlml_chl=mlml_tt(mlml_tt.fluorescence_qc_agg==1,{'fluorescence'});
mlml_chl=retime(mlml_chl,'daily');
plot(mlml_tt.datetime(mlml_tt.fluorescence_qc_agg==1),movmean(mlml_tt.fluorescence(mlml_tt.fluorescence_qc_agg==1),12*24,1,"omitmissing"))
hold on
xlim([start_date end_date])
% plot(mlml_chl.datetime,mlml_chl.fluorescence)
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
 
