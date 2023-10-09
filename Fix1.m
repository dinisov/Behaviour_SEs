%Fixation Analysis

close all; clear;
addpath('./Functions');

%% Variables
timeMin = 20;               % minimum time to include experiment
fixWindow = 20;             % +/- fixWindow means fixation
divisions = 10;              % split learning period into this many divisions
timeDiv = 30000/divisions;  % time for each division in ms

%% import data

files = dir('./Fix');
files = files(~[files.isdir],:);
datas = cell(length(files),1);
filenames = cell(length(files),1);
for i = 1:length(files)
    datas{i} = table2array(readtable(fullfile('Fix',files(i).name)));
    disp(files(i).name);
    filenames{i} = files(i).name;
end

nFlies = length(datas);
disp(nFlies);

%% collate and organise data
FLIES = struct;
Results = nan(nFlies+2,divisions);

for fly = 1:nFlies
thisFlyData = datas{fly};
numRows = length(thisFlyData);
%disp(numRows);
FixArray = zeros(3,divisions+1);

% Check for start of flight
timeStart = thisFlyData(1,2);
rowStart = 1;
numZeros = 0;
    for i = 1:numRows
        if(thisFlyData(i,3)==0)
            numZeros = numZeros +1;
        else
            if(numZeros > 100)
                timeStart = thisFlyData(i,2) ;
                rowStart = i;
            end
            numZeros = 0;
        end
    end
timeDelta = timeStart - thisFlyData(1,2);
%disp(timeDelta);
%disp(rowStart);

% Get percentage of loops that pos is within +/-fixWindow for each division
    for i = rowStart:numRows
        timeRow = thisFlyData(i,2);
        pos = thisFlyData(i,4);
        timeDelta = timeRow - timeStart;
        div = floor(1+((timeRow - timeStart)/timeDiv));
        if(abs(pos) < fixWindow)
            FixArray(1,div) = FixArray(1,div) +1;
        else
            FixArray(2,div) = FixArray(2,div) +1;
        end
    end
    
    for i = 1:divisions % put data into results array
        FixArray(3,i) = FixArray(1,i) + FixArray(2,i);
        Results(fly,i) = FixArray(1,i)/(FixArray(1,i) + FixArray(2,i));
    end
FLIES(fly).FixArray = FixArray;     % add to struct
end

%% analyse data (this assumes a sampling rate of 200 Hz)

for fly = 1:nFlies
    meanFix = nanmean(Results);
    
    stdFix = std(Results,[],1,'omitnan');
    nFix = sum(~isnan(Results), 1);
    semFix = stdFix ./ sqrt(nFix);
end

figure
plot(meanFix); title ('Learning Fixation');

Results(nFlies+1,:) = meanFix;
Results(nFlies+2,:) = semFix;

%% Save Data
Folder = 'C:\Experiments\';
baseFilename = [num2str(fixWindow), '_deg_', num2str(divisions), '_div'];
resultsFilenameCSV = [Folder,'Fix_Results_',baseFilename,'.csv'];
writematrix(Results,resultsFilenameCSV);


%myTable = cell2table(ParamArray);
%writetable(myTable, paramFilenameCSV, 'Delimiter', ',');

