%behaviour SEs

close all; clear;

addpath('./Functions');

%% import data

files = dir('./Data2');

files = files(~[files.isdir],:);

datas = cell(length(files),1);

for i = 1:length(files)
    datas{i} = table2array(readtable(fullfile('Data2',files(i).name)));
end

nFlies = length(datas);

%% collate and organise data

FLIES = struct;

for fly = 1: nFlies

thisFlyData = datas{fly};

%there are some zero rows for some reason
% data = data(logical(sum(data,2)),:);

% data = data(data(:,10) > 0,:);

%starting points of inter-stimulus periods (careful this is actually the index just before the start)
% last trial is usually interrputed so get rid of last four (this may fail is experiment is interrputed during open loop)
interStimulusStart = find(diff(thisFlyData(:,4) == 200) == 1); interStimulusStart = interStimulusStart(1:end-4);
interStimulusEnd = find(diff(thisFlyData(:,4) == 200) == -1); interStimulusEnd = interStimulusEnd(1:end-4);

nTrials = (length(interStimulusStart)/4);

%%  infer trial sequence
stimuli = zeros(nTrials,5);

% figure out first four (open loop) stimuli (gets rid of last four inte)
stimuli(:,1:4) = reshape(sign(thisFlyData(interStimulusStart,4)),[4,nTrials]).';

%figure out last (close loop) stimulus
closedLoopStart = interStimulusEnd+1; 
closedLoopStart = closedLoopStart(4:4:end);%index of start of closed loop
stimuli(:,5) = sign(thisFlyData(closedLoopStart,4));

%% add data to structure for each fly
FLIES(fly).interStimulusStart = interStimulusStart;
FLIES(fly).interStimulusEnd = interStimulusEnd;
FLIES(fly).nTrials = nTrials;
FLIES(fly).stimuli = stimuli;
FLIES(fly).time = thisFlyData(:,2);
FLIES(fly).barPosition = thisFlyData(:,4);
FLIES(fly).closedLoopStart = closedLoopStart;

end

%% analyse data (this assumes a sampling rate of 200 Hz)

% matrix to put results in
R = struct;

% groups (1,32),(2,31),(3,30), etc, as representing the same pattern
% avoids very costly flip() operations later
auxSeq = [1:16 16:-1:1];

for fly = 1:nFlies
    
    thisFly = FLIES(fly);
    
    %matrix to put data separated by trial and sequence
    %800 here assumes 4 seconds at 200 Hz
    tracesSEQ = zeros(800,16,nTrials);
    
    %matrix to store fixation times (calculated for each epoch individually)
    fixTimes = zeros(nTrials,16);
    
    %SCOTT: YOU CAN ADD A MATRIX FOR A NEW MEASURE HERE
    %scottsMatrix = zeros(nTrials,16);
    
    % vector to index bad trials
    badTrials = zeros(nTrials,1);
    
    for trial = 1:thisFly.nTrials
        seq = auxSeq(bin2dec(num2str(thisFly.stimuli(trial,1:5) > 0)) + 1);
        
        thisTrace = thisFly.barPosition(thisFly.closedLoopStart(trial):thisFly.closedLoopStart(trial)+799);
        
        tracesSEQ(:,seq,trial) = thisTrace;
        
        % if the bar close to 0? SCOTT: CHECK WINDOW BOUNDS HERE
        window = [-16 16];
        betweenBounds = thisTrace > window(1) & thisTrace <  window(2);
        
        %calculate fixTime
        if any(betweenBounds)
            fixTimes(trial, seq) = find(betweenBounds,1)*0.005; % again assumes a 200 Hz panel
            if any(diff(thisTrace) > 150)% if it wrapped around before entering window it's a bad trial
                badTrials(trial) = 1;
            end
        else %did not even enter the window
            badTrials(trial) = 1;
        end
        
        % YOU CAN ADD A NEW CALCULATION FOR A DIFFERENT MEASURE HERE
        %scottsMatrix(trial, seq) = *some calculation*
            
    end

    % better reorder here so everything downstream is according to
    % literature
    R(fly).tracesSEQ = tracesSEQ(:,seq_eff_order(5),:);
    R(fly).fixTimes = fixTimes(:,seq_eff_order(5));
    R(fly).badTrials = badTrials;
   
end

%% calculate SEs Scott's way

%SCOTT: THIS SHOULD BE EASY TO MODIFY SO MULTIPLE FLIES ARE CONCATENATED
%ALSO CALCULATIONS SHOULD BE ANALOGOUS FOR ANY OTHER MEASURES
for fly = 1:nFlies
   
    fixTimes = R(fly).fixTimes;
    fixTimes = fixTimes(~badTrials,:);
    fixTimes(fixTimes == 0) = NaN;
    meanFixTimes = nanmean(fixTimes);
    stdFixTimes = std(fixTimes,[],1,'omitnan');
    nFixTimes = sum(~isnan(fixTimes), 1);
    semFixTimes = stdFixTimes ./ sqrt(nFixTimes);
    
    create_seq_eff_plot(meanFixTimes.',[],'errors',semFixTimes.');
    
end