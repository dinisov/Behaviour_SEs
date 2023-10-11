%Fixation Analysis

close all; clear;
addpath('./Functions');

%% Variables
timeMin = 20;               % minimum time to include experiment
% fixWindow = 20;             % +/- fixWindow means fixation
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
plotData = zeros(divisions,3);
plotSEM = zeros(divisions,3);

fixWindow = [20, 30, 40];

% different window sizes
for w = 1:3
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
            if(abs(pos) < fixWindow(w))
                FixArray(1,div) = FixArray(1,div) +1;
            else
                FixArray(2,div) = FixArray(2,div) +1;
            end
        end

        for i = 1:divisions % put data into results array
            FixArray(3,i) = FixArray(1,i) + FixArray(2,i);
            Results(fly,i) = FixArray(1,i)/(FixArray(1,i) + FixArray(2,i));
        end
        
        FlightArray = thisFlyData(rowStart:numRows,4);
        FLIES(fly).FlightArray = FlightArray;     % add to struct
    end

    %% analyse data (this assumes a sampling rate of 200 Hz)

    for fly = 1:nFlies
        meanFix = nanmean(Results);

        stdFix = std(Results,[],1,'omitnan');
        nFix = sum(~isnan(Results), 1);
        semFix = stdFix ./ sqrt(nFix);
    end
    
    % figure
    % plot(meanFix); title ('Learning Fixation');

    Results(nFlies+1,:) = meanFix;
    Results(nFlies+2,:) = semFix;
    
    plotData(:,w) = meanFix.';
    plotSEM(:,w) = semFix.';
end

%% Plot Combined
figure; errorbar(repmat((3:3:30).',[1 3]),plotData,plotSEM); xlim([0 33]);
set(gca,'xtick',3:3:30,'ytick',.3:.1:1,'yticklabel',30:10:100); ylabel('Percentage'); xlabel('time (s)'); title('Learning Fixation');
% ylim([0 1]);

%% position Histogram of all data
positionAllFlies = cell(nFlies,1);

for fly = 1:nFlies
    positionAllFlies{fly} = FLIES(fly).FlightArray;
end

positionAllFlies = cell2mat(positionAllFlies);

% RT of all trials
fiXPosArray = nonzeros(positionAllFlies(:,:));
figure;histogram(fiXPosArray,'NumBins', 30); title('Bar Position During Learning');  xlabel('Position (deg)'); xlim([-180 180]);
set(gca,'xtick',-180:30:180); 

%% Save Data
Folder = 'C:\Experiments\';
baseFilename = [num2str(fixWindow), '_deg_', num2str(divisions), '_div'];
resultsFilenameCSV = [Folder,'Fix_Results_',baseFilename,'.csv'];
writematrix(Results,resultsFilenameCSV);


