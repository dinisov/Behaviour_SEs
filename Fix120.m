%Fixation Analysis
% Makes 4 histograms of consecutive 30s for a single 120s fixation data file

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
for k = 1:length(files)
    datas{k} = table2array(readtable(fullfile('Fix',files(k).name)));
    disp(files(k).name);
    filenames{k} = files(k).name;
end

nFlies = length(datas);
disp(nFlies);

%% Remove non-Flight

for fly = 1:nFlies
    thisFlyData = datas{fly};
    thisFlyData(:,2) = thisFlyData(:,2) - thisFlyData(1,2);
    numRows = length(thisFlyData);
    disp(numRows);
    disp(thisFlyData(1,2));
    disp(thisFlyData(numRows-1,2));
    runLength = thisFlyData(numRows-1,2);
    chunkLength = floor(runLength/4);
    
%     indexFlying = plotFlyingData(thisFlyData(:,3),100);
    figure('position',[100 100 1200 200],'color', 'w') ;
 
    for k = 1:4
        chunk_start = thisFlyData(1,2) + chunkLength * (k -1);
        chunk_end = thisFlyData(1,2) + chunkLength * (k);

        chunkIndex = (thisFlyData(:,2) > chunk_start) & (thisFlyData(:,2) < chunk_end);

        WLEChunk = thisFlyData(chunkIndex,3);

        flyingIndexChunk = plotFlyingData(WLEChunk,100);
        
        positionChunk = thisFlyData(chunkIndex,4);
        positionChunk = positionChunk(flyingIndexChunk);

        subplot(1,4,k); 
        %figure;
        %indexFlying = plotFlyingData(thisFlyData(thisFlyData(:,2) > chunk_start & thisFlyData(:,2) < chunk_end,3),50);
        %PosData = thisFlyData(indexFlying,4);
        %IndexChunk = IndexChunk & indexFlying;
        %PosData = thisFlyData(thisFlyData(:,2) > chunk_start & thisFlyData(:,2) < chunk_end,4);

        %figure;
        %histogram(thisFlyData(:,4));
        histogram(positionChunk,'NumBins', 30);
        title(['Bar Position ' num2str((k-1)*30) ' - ' num2str(k*30) ' sec']); xlim([-180 180]); ylim([0 1000]); ylabel('Counts');
        %plotFlyingData(thisFlyData(:,3),100);
    end

end

%% functions
function flyingIdx = plotFlyingData(flyingData, maxRunLength)

    flyingIdx = 1:length(flyingData);
    transitions = diff([0; flyingData == 0; 0]); %find where the array goes from non-zero to zero and vice versa
    runstarts = find(transitions == 1);
    runends = find(transitions == -1); %one past the end
    runlengths = runends - runstarts;
    %keep only those runs of length 4 or more:
    runstarts(runlengths < maxRunLength) = [];
    runends(runlengths < maxRunLength) = [];
    %expand each run into a list indices:
    notFlyingIdx = arrayfun(@(s, e) s:e-1, runstarts, runends, 'UniformOutput', false);
    notFlyingIdx = [notFlyingIdx{:}];  %concatenate the list of indices into one vector
    flyingData(notFlyingIdx) = [];
    flyingIdx(notFlyingIdx) = [];

%     histogram(flyingData);
end