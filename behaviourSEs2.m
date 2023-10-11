%behaviour SEs

close all; clear;

addpath('./Functions');
%% parameters
limitRTB = 150;             %degrees
windowRT = 30;              %degrees
refixQualTime = 200;        %loops
endTrace = 792;             %loops
limitQual = windowRT +10;   %degrees

%% import data
files = dir('./Data2');
files = files(~[files.isdir],:);
datas = cell(length(files),1);

for i = 1:length(files)
    datas{i} = table2array(readtable(fullfile('Data2',files(i).name)));
%     disp(files(i).name);
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
% last trial is usually interrupted so get rid of last four (this may fail is experiment is interrupted during open loop)
interStimulusStart = find(diff(thisFlyData(:,4) == 200) == 1); %interStimulusStart = interStimulusStart(1:end-4);
interStimulusEnd = find(diff(thisFlyData(:,4) == 200) == -1); %interStimulusEnd = interStimulusEnd(1:end-4);

nTrials = (length(interStimulusStart)/4);

% littleBit = interStimulusStart(2) - interStimulusEnd(1);
% 
% openLoopStart = interStimulusStart - littleBit;
% openLoopStart = reshape(openLoopStart,[4,nTrials]);
% openLoopStart = openLoopStart(1,:);
% 
% figure; plot(FLIES(end).barPosition); hold on; plot(openLoopStart,zeros(size(openLoopStart)),'.','markersize',20);

%%  infer trial sequence
stimuli = zeros(nTrials,5);

% figure out first four (open loop) stimuli
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

% groups (1,32),(2,31),(3,30), etc, as representing the same pattern
% avoids very costly flip() operations later
% auxSeq = [1:16 16:-1:1];
for fly = 1:nFlies
    
    thisFly = FLIES(fly);
    nTrials = thisFly.nTrials;
    %matrix to put data separated by trial and sequence
    %800 here assumes 4 seconds at 200 Hz
    tracesSEQ = zeros(endTrace,16,nTrials);
    
    %matrix to store fixation times (calculated for each epoch individually)
<<<<<<< Updated upstream
    fixTimes = zeros(nTrials,16);
    
=======
    fixTimes = zeros(nTrials,32);
    timeouts = zeros(nTrials,32);
    RTB = zeros(nTrials,32);
    Qual = zeros(nTrials,32);
    saveArray = zeros(nTrials,12);    
>>>>>>> Stashed changes
    %SCOTT: YOU CAN ADD A MATRIX FOR A NEW MEASURE HERE
    %scottsMatrix = zeros(nTrials,16);
    
    % vector to index bad trials
    goodTrial = 0;
    badTrials = zeros(nTrials,1);
    
    for trial = 1:nTrials
        seq = bin2dec(num2str(thisFly.stimuli(trial,1:5) > 0)) + 1;
        saveArray(trial,3:7) = thisFly.stimuli(trial,1:5); %put stimuli in
        saveArray(trial,2) = thisFly.closedLoopStart(trial); %put loop number of closed loop start in
        thisTrace = thisFly.barPosition(thisFly.closedLoopStart(trial):thisFly.closedLoopStart(trial)+ endTrace - 1);
        
        tracesSEQ(:,seq,trial) = thisTrace;
        
        % if the bar close to 0? SCOTT: CHECK WINDOW BOUNDS HERE
        window = [-windowRT windowRT];
        betweenBounds = thisTrace > window(1) & thisTrace <  window(2);
        
        %calculate fixTime
        if any(betweenBounds)
            loopRT = find(betweenBounds,1);
            %fixTimes(trial, seq) = loopRT*0.005; % again assumes a 200 Hz panel
            %test = thisFly.time(loopRT)- FLIES(fly).interStimulusEnd(trial);
            fixTimes(trial, seq) = thisFly.time(thisFly.closedLoopStart(trial)+loopRT-1)- thisFly.time(thisFly.closedLoopStart(trial));
            %fprintf('test %.0f,%0.f' , thisFly.time(thisFly.closedLoopStart(trial)),test);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            saveArray(trial,8) = fixTimes(trial, seq);  %put fixtime in
            rtTrace = thisTrace(1:loopRT);
            %disp(rtTrace);
            if(loopRT + refixQualTime > endTrace)  
                qualTrace = thisTrace(loopRT:endTrace);
            else
                qualTrace = thisTrace(loopRT:loopRT + refixQualTime);
            end
            if any(abs(rtTrace) > limitRTB)% if it wrapped around before entering window it's a bad trial
                badTrials(trial) = 1;
                RTB(trial, seq) = 1;
                saveArray(trial,9) = 1;
                saveArray(trial,11) = 1; %put 1 in bad trial and RTB (round the back)
            elseif any(abs(qualTrace) > limitQual)% if it wrapped around before entering window it's a bad trial
                badTrials(trial) = 1;
                Qual(trial, seq) = 1;
                saveArray(trial,9) = 1;
                saveArray(trial,12) = 1; %put 1 in bad trial and RTB (round the back)
            end
        else %did not even enter the window
            badTrials(trial) = 1;
            timeouts(trial, seq) = 1;
            saveArray(trial,9:10) = 1;  %put 1 in bad trial and timeout
        end
        
        if(badTrials(trial) ~= 1)
            goodTrial = goodTrial +1;
        end
        
        % YOU CAN ADD A NEW CALCULATION FOR A DIFFERENT MEASURE HERE
        %scottsMatrix(trial, seq) = *some calculation*
            
    end
    
%     disp('Total Trials');
%     disp(thisFly.nTrials);
%     disp('good Trials');
%     disp(goodTrial);
    
    % get rid of bad trials here so no need to do it downstream
    fixTimes = fixTimes(~badTrials,:);
    tracesSEQ = tracesSEQ(:,:,~badTrials);
    
    FLIES(fly).fixTimes32 = fixTimes;% store original 32-long matrix for isomer calculations
    
    % group by pattern and reduce to 16 sequences
    tracesSEQ = tracesSEQ + flip(tracesSEQ,2);
    tracesSEQ = tracesSEQ(:,1:16,:);
    
    fixTimes = fixTimes + flip(fixTimes,2);
    fixTimes = fixTimes(:,1:16);

    % better reorder here so everything downstream is according to
    % literature
    FLIES(fly).tracesSEQ = tracesSEQ(:,seq_eff_order(5),:);
    FLIES(fly).fixTimes = fixTimes(:,seq_eff_order(5));
    FLIES(fly).badTrials = badTrials;
   
end
%% Analyse Isomers
medianRT = zeros(nFlies,4);
for fly = 1:nFlies
    Hist_Array_L = nonzeros(FLIES(fly).fixTimes32(:,1:16));
    Hist_Array_R = nonzeros(FLIES(fly).fixTimes32(:,17:32));
    Hist_Array_All = nonzeros(FLIES(fly).fixTimes32(:,:));
    numBins = 10;
%         subplot(2,2,1);
%         histogram(Hist_Array_L,'NumBins', numBins);
%         title('Left Isomer RT');
%     subplot(2,2,2);
%         histogram(Hist_Array_R,'NumBins', numBins);
%         title('Right Isomer RT');
%     subplot(2,2,3);
%         histogram(Hist_Array_All,'NumBins', numBins);
%         title('All RT');
    medianL = nanmedian(Hist_Array_L);        medianRT(fly,1) = medianL;%disp(medianL);
    medianR = nanmedian(Hist_Array_R);        medianRT(fly,2) = medianR;%disp(medianR);
    medianAll = nanmedian(Hist_Array_All);  medianRT(fly,3) = medianAll;%disp(medianAll);
    medianRT(fly,4) = medianRT(fly,2) - medianRT(fly,1);
end
% concatenate
figure;
Hist_Array_L = nonzeros(medianRT(:,1));
Hist_Array_R = nonzeros(medianRT(:,2));
Hist_Array_All = nonzeros(medianRT(:,3));
Hist_Array_Diff = nonzeros(medianRT(:,4));
    numBins = 10;
        subplot(2,2,1);
        histogram(Hist_Array_L,'NumBins', numBins);
        title('Left Isomer Median RT'); xlim([0 2000]);
    subplot(2,2,2);
        histogram(Hist_Array_R,'NumBins', numBins);
        title('Right Isomer MedianRT'); xlim([0 2000]);
    subplot(2,2,3);
        histogram(Hist_Array_All,'NumBins', numBins);
        title('All Median RT'); xlim([0 2000]);
    subplot(2,2,4);
        histogram(Hist_Array_Diff,'NumBins', numBins);
        title('Isomer Diff Median RT'); xlim([-600 600]);
        
%% calculate SEs for each fly individually

for fly = 1:nFlies
    
    R(fly) = calculateSEsBehav(FLIES(fly).fixTimes); %#ok<SAGROW>
    
    %figure; create_seq_eff_plot(R(fly).meanFixTimes.',[],'errors',R(fly).semFixTimes.'); title(['Mean Fix Time Fly ' num2str(fly)]);
    %figure; create_seq_eff_plot(R(fly).medianFixTimes.',[],'errors',R(fly).semFixTimes.'); title(['Median Fix Time Fly ' num2str(fly)]);
    
    %plotIsomersBehav(FLIES(fly).fixTimes32);
    
    %SEAnovaBehav(FLIES(fly).fixTimes);
     
end
%% Concatenate
if(nFlies > 1)
    % calculate SEs for all flies using superfly method (concatenate)
    fixTimesAllFlies = cell(nFlies,1);
    fixTimesAllFlies32 = cell(nFlies,1);
    badTrialsAllFlies = cell(nFlies,1);
    for fly = 1:nFlies
        fixTimesAllFlies{fly} = FLIES(fly).fixTimes;
        fixTimesAllFlies32{fly} = FLIES(fly).fixTimes32;
    end

    fixTimesAllFlies = cell2mat(fixTimesAllFlies);
    fixTimesAllFlies32 = cell2mat(fixTimesAllFlies32);
    
    % RT of all trials
    RT_Hist_Array = nonzeros(fixTimesAllFlies(:,:));
    figure;histogram(RT_Hist_Array,'NumBins', 20); title('All RT');  xlabel('Time (ms)');

    RAllFlies = calculateSEsBehav(fixTimesAllFlies);

    %figure; create_seq_eff_plot(RAllFlies.meanFixTimes.',[],'errors',RAllFlies.semFixTimes.');
    %figure; create_seq_eff_plot(RAllFlies.medianFixTimes.',[],'errors',RAllFlies.semFixTimes.');

    %plotIsomersBehav(fixTimesAllFlies32);

    %SEAnovaBehav(fixTimesAllFlies);

    %% calculate SEs for all flies using averaging of profiles method (no weighting)
    allProfilesMean = zeros(nFlies,16);
    allProfilesMedian = zeros(nFlies,16);

    for fly = 1:nFlies
        allProfilesMean(fly,:) = R(fly).meanFixTimes;
        allProfilesMedian(fly,:) = R(fly).medianFixTimes;
    end

    % this does not weigh by the number of data points for each fly
    % taking the mean of the means or medians here as the median of a few
    % points is heavily skewed; error is just std of the means/medians
    %figure; create_seq_eff_plot(nanmean(allProfilesMean).',[],'errors',std(allProfilesMean,[],1,'omitnan').'./sqrt(nFlies));
    %figure; create_seq_eff_plot(nanmean(allProfilesMedian).',[],'errors',std(allProfilesMedian,[],1,'omitnan').'./sqrt(nFlies));

    %SEAnovaBehav(allProfilesMean);

    %% calculate SEs for all flies using averaging of profiles method (weighted, propagated error)
    allProfilesMean = zeros(nFlies,16);
    allProfilesMedian = zeros(nFlies,16);
    semFixTimesAllFlies = zeros(nFlies,16);
    nFixTimesAllFlies = zeros(1,nFlies);
    for fly = 1:nFlies
        nFixTimesAllFlies(fly) = sum(R(fly).nFixTimes);
        allProfilesMean(fly,:) = R(fly).meanFixTimes * nFixTimesAllFlies(fly);
        allProfilesMedian(fly,:) = R(fly).medianFixTimes * nFixTimesAllFlies(fly);

        % part of the propagation error calculation
        semFixTimesAllFlies(fly,:) = R(fly).semFixTimes.^2 * nFixTimesAllFlies(fly).^2;
    end

    % finish propagation of error calculation
    semFixTimesAllFlies = sqrt(sum(semFixTimesAllFlies/(sum(nFixTimesAllFlies)^2),1));

    figure; create_seq_eff_plot(nansum(allProfilesMean).'./sum(nFixTimesAllFlies),[],'errors',semFixTimesAllFlies.'); title('SE weighted Mean');
    figure; create_seq_eff_plot(nansum(allProfilesMedian).'./sum(nFixTimesAllFlies),[],'errors',semFixTimesAllFlies.');title('SE weighted Median');

    %SEAnovaBehav(allProfilesMean);
end
