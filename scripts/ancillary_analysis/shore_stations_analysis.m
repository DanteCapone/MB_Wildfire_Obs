%Shore Stations Analysis


%% Load Other shore stations

%Humboldt

%BML


%Tiburon

%MLML














data=BMLWTSc877a3e9397a;
data.time=char(data.time);
year=str2num(data.time(:,1:4));
month=str2num(data.time(:,6:7));
day=str2num(data.time(:,9:10));
hour=str2num(data.time(:,12:13));
minute=str2num(data.time(:,15:16));
second=str2num(data.time(:,18:19));
data.datetime=datetime(year,month,day,hour,minute,second);
data.datetime.TimeZone = 'America/Los_Angeles';

bml_obs_tt=table2timetable(data);


%%
plot(bmlseawaterfluorescence2018daily.datetime,bmlseawaterfluorescence2018daily.chla)

%%

plot(bml_obs_tt.datetime,bml_obs_tt.chl_conc)


%%
folderPath = 'C:\Users\Dante Capone\OneDrive\Desktop\Scripps_PhD\Wildfire_Obs';
filePattern = 'bml_seawater_fluorescence_*.csv'; % assuming CSV files
finalTable = concatenateFilesFromFolder(folderPath, filePattern, @readtable);



%%
% Generate date ticks every 4 months
dates=finalTable.datetime;
start_date = dates(1);
end_date = dates(end);
date_ticks = start_date:calmonths(2):end_date;

plot(finalTable.datetime,finalTable.chla)


% Set the XTick property and format using datetick
ax = gca;
ax.XTick = date_ticks;
datetick('x','mmm-yyyy','keepticks');


%% CLimatology
plotClimatology(finalTable)
finalTable.dayOfYear=day(finalTable.datetime, 'dayofyear');
legend_entries=["Climatology"];
for i=2014:2023
    if i==2017
        plot(finalTable.dayOfYear(finalTable.datetime.Year==i),finalTable.chla(finalTable.datetime.Year==i),'Color','r','LineWidth',2)

    else
        plot(finalTable.dayOfYear(finalTable.datetime.Year==i),finalTable.chla(finalTable.datetime.Year==i),'Color',rand(1, 3))
    end
    hold on
    legend_entry=num2str(i);
    legend_entries=[legend_entries; legend_entry];
end

legend(legend_entries)
set(gcf,[0 0 1000 800])



%% Tiburon

tiburon=tiburonwatertibc18fbd3b262ba1;
tiburon.time=char(tiburon.time);
year=str2num(tiburon.time(:,1:4));
month=str2num(tiburon.time(:,6:7));
dayz=str2num(tiburon.time(:,9:10));
hour=str2num(tiburon.time(:,12:13));
minute=str2num(tiburon.time(:,15:16));
second=str2num(tiburon.time(:,18:19));
tiburon.datetime=datetime(year,month,dayz,hour,minute,second);
tiburon.datetime.TimeZone = 'America/Los_Angeles';

tiburon_tt=table2timetable(tiburon);
tiburon_tt.dayOfYear = day(tiburon_tt.datetime, 'dayofyear');

%% Plot tibuton chl a
plot(tiburon_tt.datetime,tiburon_tt.mass_concentration_of_chlorophyll_in_sea_water)

close all;
figure()
plotClimatology(tiburon_tt,'mass_concentration_of_chlorophyll_in_sea_water')
hold on
plot(tiburon_tt.dayOfYear(tiburon_tt.datetime.Year==2020),tiburon_tt.mass_concentration_of_chlorophyll_in_sea_water(tiburon_tt.datetime.Year==2020),'r')

%%
function plotClimatology(dataTable,input_var)
    % Extract day of the year for each datetime entry
    dayOfYear = day(dataTable.datetime, 'dayofyear');
    dataTable.dayOfYear = dayOfYear;

    % Group by dayOfYear and compute the mean for each group
    climatology = varfun(@nanmean, dataTable, 'GroupingVariables', 'dayOfYear', 'InputVariables', input_var);

    % Plotting
    figure;
    plot(climatology.dayOfYear, eval(['climatology.nanmean_' input_var]),'Color','k','LineStyle','-','LineWidth',2);
    hold on
    xlabel('Day of Year');
    ylabel('Average of Chl A');
    title('Climatology of Chl A');
    datetick('x')
    grid on;
end


    





function finalTable = concatenateFilesFromFolder(folderPath, filePattern, readFunction)
    % folderPath: path to the directory containing files
    % filePattern: pattern to match files, e.g., '*.csv' for CSV files
    % readFunction: handle to a function that reads a file and returns a table

    % Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 4);
    
    % Specify range and delimiter
    opts.DataLines = [1, Inf];
    opts.Delimiter = ",";
    
    % Specify column names and types
    opts.VariableNames = ["VarName1", "datetime", "VarName3", "chla"];
    opts.VariableTypes = ["double", "datetime", "double", "double"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    % Specify variable properties
    opts = setvaropts(opts, "datetime", "InputFormat", "yyyy-MM-dd HH:mm:ss");

    % List all files in the directory
    fileList = dir(fullfile(folderPath, filePattern));

    % Cell array to store individual tables
    allTables = cell(length(fileList), 1);

    % Loop through each file and read data
    for k = 1:length(fileList)
        filePath = fullfile(folderPath, fileList(k).name);
        allTables{k} = readFunction(filePath, opts);
    end

    % Concatenate all tables into a single table
    finalTable = vertcat(allTables{:});
end








