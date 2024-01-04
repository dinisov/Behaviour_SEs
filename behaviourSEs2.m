%behaviour SEs
close all; clear;
addpath('./Functions');

%% parameters
limitAngle = 90;            %degrees
windowRT = 45;              %degrees
refixQualTime = 300;        %loops
endTrace = 792;             %loops
limitReQual = windowRT +0;   %degrees
timeOut = 4000;             %ms

enRT = 1;                   %Enable RT PLots
enMedRT = 0;                %Enable median isomer RT PLots
enQTplot = 1;               %Enable QT vs trial plot
enPlot4x4indiv =0;          %Enable 4x4 histogram plot individual experiments
enPlot4x4super = 1;         %Enable 4x4 histogram plot of super fly

%% import data
files = dir('./Data2');
files = files(~[files.isdir],:);
datas = cell(length(files),1);

for i = 1:length(files)
    datas{i} = table2array(readtable(fullfile('Data2',files(i).name)));
    disp(files(i).name);
end

nFlies = length(datas);
disp(nFlies);

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

%% work out Qualification Time
qualTime = zeros(nTrials,1);
qualTimeAbs = thisFlyData(closedLoopStart,2) - thisFlyData(1,2) - 667;
qualTime(1) = qualTimeAbs(1);
for i = 1:nTrials-1
    qualTime(i+1) = qualTimeAbs(i+1) - qualTimeAbs(i) -4667;
end
% qualTimeArray = nonzeros(qualTime);
% figure; histogram(qualTimeArray);
% figure; plot(qualTimeArray);

%% add data to structure for each fly
FLIES(fly).interStimulusStart = interStimulusStart;
FLIES(fly).interStimulusEnd = interStimulusEnd;
FLIES(fly).nTrials = nTrials;
FLIES(fly).stimuli = stimuli;
FLIES(fly).time = thisFlyData(:,2);
FLIES(fly).barPosition = thisFlyData(:,4);
FLIES(fly).closedLoopStart = closedLoopStart;
FLIES(fly).qualTime = qualTime;
end

%% plot Qualtimes 
qualTimeArray = zeros(120, nFlies);
for fly = 1:nFlies
    qualTimeArray(1:length(FLIES(fly).qualTime),fly) = FLIES(fly).qualTime;
end
qualTimeArray(~logical(qualTimeArray)) = nan;
%figure; plot(qualTimeArray);
%title('Qualification Time (ms)'); ylim([0 5000]);

aux_x = cell(nFlies,1);
aux_y = cell(nFlies,1);

if(enQTplot)    % switch at top of file
    figure; hold on;
    for i = 1:nFlies
        plot(qualTimeArray(qualTimeArray(:,i) < 60000,i),'b.');
        aux_y{i} = qualTimeArray(qualTimeArray(:,i) < 60000,i);
        aux_x{i} = (1:length(aux_y{i})).';
    end

    x = cell2mat(aux_x); y = cell2mat(aux_y);
    [p,s] = polyfit(x,y,1);
    f = polyval(p,1:120);

    plot(1:120,f,'r','LineWidth',2);
    xlabel('Trial'); ylabel('Qualification Time (ms)');title('Qualification Time vs trial'); ylim([0 5000]); text(60,700,['y = ',num2str(p(1),3),'x + ',num2str(p(2),'%.0f')],'FontSize',14);

% figure; plot(nanmedian(qualTimeArray,2));
end
%% analyse data (this assumes a sampling rate of 200 Hz)

% groups (1,32),(2,31),(3,30), etc, as representing the same pattern
% avoids very costly flip() operations later
% auxSeq = [1:16 16:-1:1];
for fly = 1:nFlies
    
    thisFly = FLIES(fly);
    nTrials = thisFly.nTrials;
    %matrix to put data separated by trial and sequence
    %endTrace here assumes 4 seconds at 200 Hz - needs to be less than 800 as loops a bit over 5ms
    tracesSEQ = zeros(endTrace,16,nTrials);
    
    %matrix to store fixation times (calculated for each epoch individually)
    fixTimes = zeros(nTrials,32);
    timeOuts = zeros(nTrials,32);
    RTB = zeros(nTrials,32);
    reQual = zeros(nTrials,32);
    saveArray = zeros(nTrials,12);    
    
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
        outsideMaxAngle = thisTrace > limitAngle & thisTrace <  -limitAngle;%thisTrace > limitAngle | thisTrace <  -limitAngle; turn on this feature
        
        %calculate fixTime
        if any(betweenBounds)
            loopRT = find(betweenBounds,1); %Loop index of first betweenbounds
            fixTimes(trial, seq) = thisFly.time(thisFly.closedLoopStart(trial)+loopRT-1)- thisFly.time(thisFly.closedLoopStart(trial)); %Time RT in ms
            saveArray(trial,8) = fixTimes(trial, seq);  %put fixtime in
            rtTrace = thisTrace(1:loopRT);
            
            if(loopRT + refixQualTime > endTrace)           % set up Refix Window time
                reQualTrace = thisTrace(loopRT:endTrace);
            else
                reQualTrace = thisTrace(loopRT:loopRT + refixQualTime);
            end
            if any(abs(rtTrace) > limitAngle)% if it wrapped around before entering window it's a bad trial
                badTrials(trial) = 1;
                RTB(trial, seq) = 1;
                saveArray(trial,9) = 1;
                saveArray(trial,11) = 1; %put 1 in bad trial and RTB (round the back)
            elseif any(abs(reQualTrace) > limitReQual)% if it did not stay in the requal window for limitQual it's a bad trial
                badTrials(trial) = 1;
                reQual(trial, seq) = 1;
                saveArray(trial,9) = 1;
                saveArray(trial,12) = 1; %put 1 in bad trial and RTB (round the back)
            end
        else %did not even enter the window
            if any(outsideMaxAngle)
                RTB(trial, seq) = 1;
                saveArray(trial,9) = 1;
                saveArray(trial,11) = 1; %put 1 in bad trial and RTB (round the back)
            else
                timeOuts(trial, seq) = 1;
                saveArray(trial,9:10) = 1;  %put 1 in bad trial and timeout
            end
            badTrials(trial) = 1;
        end
        
        if(fixTimes(trial, seq) > timeOut)
            badTrials(trial) = 1;
            timeOuts(trial, seq) = 1;
            saveArray(trial,9:10) = 1;  %put 1 in bad trial and timeout
        end
        
        if(badTrials(trial) ~= 1)
            goodTrial = goodTrial +1;
        end
        
        % YOU CAN ADD A NEW CALCULATION FOR A DIFFERENT MEASURE HERE
        %scottsMatrix(trial, seq) = *some calculation*
            
    end
    
    RTB = RTB + flip(RTB,2);
    RTB = RTB(:,1:16);
    FLIES(fly).RTB = sum(RTB(:,seq_eff_order(5))).';
    
    reQual = reQual + flip(reQual,2);
    reQual = reQual(:,1:16);
    FLIES(fly).reQual = sum(reQual(:,seq_eff_order(5))).';
    
    timeOuts = timeOuts + flip(timeOuts,2);
    timeOuts = timeOuts(:,1:16);
    FLIES(fly).timeOuts = sum(timeOuts(:,seq_eff_order(5))).';
    
    FLIES(fly).goodTrial = goodTrial;  %add badtrial data to Struct
%     FLIES(fly).RTB = RTB;
%     FLIES(fly).reQual = reQual;
%     FLIES(fly).timeOuts  = timeOuts;
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

%% Plot discarded trials vs Sequence
figure;                 
create_seq_eff_plot([sum([FLIES.RTB],2) sum([FLIES.timeOuts],2) sum([FLIES.reQual],2)],[]);
legend('MaxAngle Fail','Timeout Fails','Requal Fails','Location','northeast','FontSize',8); 
title('Discarded Trials vs Stimulus Sequence');

%% make 4x4 histogram plot for each experiment
if(enPlot4x4indiv)
    load('binomial_x_labels_latex_alt_rep.mat','binomial_x_labels_latex');

    for fly = 1:length(FLIES)
        figure;
        for i = 1:16
           subplot(4,4,i); histogram(FLIES(fly).fixTimes(logical(FLIES(fly).fixTimes(:,i)),i),10); 
           title(binomial_x_labels_latex{seq}(ind_horiz));
        end
    end
end

%% Analyse Badtrials - get stats
totalTrials = zeros(nFlies,5);
for fly = 1:nFlies
    totalTrials(fly,1) = FLIES(fly).nTrials; 
    totalTrials(fly,2) = FLIES(fly).goodTrial; 
    totalTrials(fly,3) = sum(FLIES(fly).RTB,"all");
    totalTrials(fly,4) = sum(FLIES(fly).reQual,"all");
    totalTrials(fly,5) = sum(FLIES(fly).timeOuts,"all");
end
total = sum(totalTrials, 1);
numTotal = sum(totalTrials(fly,1), "all");
numGood = sum(totalTrials(fly,2));
numRTB = sum(totalTrials(fly,3));
numreQual = sum(totalTrials(fly,4));
numTimeout = sum(totalTrials(fly,5));

fprintf('Total Trials %d Good Trials %d RTB %d ReQual Fails %d Timeouts %d\n',total(1),total(2),total(3),total(4),total(5));

%% Analyse Isomers
medianRT = zeros(nFlies,4);
for fly = 1:nFlies
    Hist_Array_L = nonzeros(FLIES(fly).fixTimes32(:,1:16));
    Hist_Array_R = nonzeros(FLIES(fly).fixTimes32(:,17:32));
    Hist_Array_All = nonzeros(FLIES(fly).fixTimes32(:,:));

    medianL = nanmedian(Hist_Array_L);        medianRT(fly,1) = medianL;%disp(medianL);
    medianR = nanmedian(Hist_Array_R);        medianRT(fly,2) = medianR;%disp(medianR);
    medianAll = nanmedian(Hist_Array_All);  medianRT(fly,3) = medianAll;%disp(medianAll);
    medianRT(fly,4) = medianRT(fly,2) - medianRT(fly,1);
end
% plot histograms of medians
if(enMedRT)
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
end    
%% calculate SEs for each fly individually
for fly = 1:nFlies
    
    R(fly) = calculateSEsBehav(FLIES(fly).fixTimes); %#ok<SAGROW>
    
    %figure; create_seq_eff_plot(R(fly).meanFixTimes.',[],'errors',R(fly).semFixTimes.'); title(['Mean Fix Time Fly ' num2str(fly)]);
    %figure; create_seq_eff_plot(R(fly).medianFixTimes.',[],'errors',R(fly).semFixTimes.'); title(['Median Fix Time Fly ' num2str(fly)]);
    
    %plotIsomersBehav(FLIES(fly).fixTimes32);
    
    %SEAnovaBehav(FLIES(fly).fixTimes);
end
%% Concatenate
if(nFlies > 1)      % fails if only one fly
    %% calculate SEs for all flies using superfly method (concatenate)
    fixTimesAllFlies = cell(nFlies,1);
    fixTimesAllFlies32 = cell(nFlies,1);
    badTrialsAllFlies = cell(nFlies,1);
    for fly = 1:nFlies
        fixTimesAllFlies{fly} = FLIES(fly).fixTimes;
        fixTimesAllFlies32{fly} = FLIES(fly).fixTimes32;
    end

    fixTimesAllFlies = cell2mat(fixTimesAllFlies);
    fixTimesAllFlies32 = cell2mat(fixTimesAllFlies32);
    
    RAllFlies = calculateSEsBehav(fixTimesAllFlies);

%    figure; create_seq_eff_plot(RAllFlies.meanFixTimes.',[],'errors',RAllFlies.semFixTimes.');
    figure; create_seq_eff_plot(RAllFlies.medianFixTimes.',[],'errors',RAllFlies.semFixTimes.');title('SE super Median');

    plotIsomersBehav(fixTimesAllFlies32);

%     SEAnovaBehav(fixTimesAllFlies);

%% Plots of Reaction time Left, Right, All and by sequence    
    if(enRT)
        % RT of all trials
        RT_Hist_Array = nonzeros(fixTimesAllFlies);
        meanAll = nanmean(RT_Hist_Array); 
        medianALL = nanmedian(RT_Hist_Array);
        disp([' mean ',num2str(meanAll),' median ',num2str(medianALL)]);
        RT_Hist_Array32 = nonzeros(fixTimesAllFlies32);
        RT_Hist_Array_Left = nonzeros(fixTimesAllFlies32(:,1:16));
        RT_Hist_Array_Right = nonzeros(fixTimesAllFlies32(:,17:32));
        
        RT_Hist_Array_Left_log = log10(RT_Hist_Array_Left);
        RT_Hist_Array_Right_log = log10(RT_Hist_Array_Right);
         [h,p] = ttest2(RT_Hist_Array_Left_log,RT_Hist_Array_Right_log);
         disp(['h ',num2str(h),' p ',num2str(p)]);
        meanLeftLog = nanmean(RT_Hist_Array_Left_log); 
        medianLeftLog = nanmedian(RT_Hist_Array_Left_log);
        meanRightLog = nanmean(RT_Hist_Array_Right_log);
        medianRightLog = nanmedian(RT_Hist_Array_Right_log);
        

         figure; histogram(RT_Hist_Array,'NumBins', 25); title('All RT');  xlabel('Time (ms)'); xlim([0,4200]);
         figure; histogram(RT_Hist_Array_Left,'NumBins', 25,'facecolor',[0.4609 0.7070 0.7695]); title('Left RT');  xlabel('Time (ms)'); xlim([0,4200]);
         figure; histogram(RT_Hist_Array_Right,'NumBins', 25,'facecolor',[0.9140 0.7109 0.4609]); title('Right RT');  xlabel('Time (ms)'); xlim([0,4200]); 
         
         figure('color','w','position',[10 10 1200 400]); 
         subplot(1,2,1);histogram(RT_Hist_Array_Left_log,'NumBins', 25,'facecolor',[0.4609 0.7070 0.7695]); title('Left RT');  xlabel('Time (Log sec)'); ylim([0,90]); ylabel('Count'); xlim([2,4]);
         line([meanLeftLog,meanLeftLog],[0,90],'linewidth',2,'color','k','LineStyle','-'); text(3.3,85,['mean=',num2str(meanLeftLog,3)],'FontSize',14);
         line([medianLeftLog,medianLeftLog],[0,90],'linewidth',2,'color','r','LineStyle','--'); text(3.3,78,['median=',num2str(medianLeftLog,3)],'FontSize',14);
         legend('Left RT','mean','median','Location','northwest'); 
         
         subplot(1,2,2); histogram(RT_Hist_Array_Right_log,'NumBins', 25,'facecolor',[0.9140 0.7109 0.4609]); title('Right RT');  xlabel('Time (Log sec)'); ylim([0,90]); ylabel('Count'); xlim([2,4]);
         line([meanRightLog,meanRightLog],[0,90],'linewidth',2,'color','k','LineStyle','-'); text(3.3,85,['mean=',num2str(meanRightLog,3)],'FontSize',14);
         line([medianRightLog,medianRightLog],[0,90],'linewidth',2,'color','r','LineStyle','--'); text(3.3,78,['median=',num2str(medianRightLog,3)],'FontSize',14);
         legend('Right RT','mean','median','Location','northwest'); 
        
    end
    

    if(enPlot4x4super)
        load('binomial_x_labels_latex_alt_rep.mat','binomial_x_labels_latex');
        ind_horiz = sub2ind(size(binomial_x_labels_latex{1}),1:4,[1 1 1 5]);
        figure;
        for seq = 1:16
           subplot(4,4,seq); histogram(fixTimesAllFlies(logical(fixTimesAllFlies(:,seq)),seq),10,'NumBins', 16); xlim([0,4000]);ylim([0,40]);
           title(binomial_x_labels_latex{seq}(ind_horiz));
        end
    end
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
    semFixTimesAllFlies = sqrt(nansum(semFixTimesAllFlies/(sum(nFixTimesAllFlies)^2),1));

    figure; create_seq_eff_plot(nansum(allProfilesMean).'./sum(nFixTimesAllFlies),[],'errors',semFixTimesAllFlies.'); title('SE weighted Mean'); %ylim([0 1000]);
    figure; create_seq_eff_plot(nansum(allProfilesMedian).'./sum(nFixTimesAllFlies),[],'errors',semFixTimesAllFlies.');title('SE weighted Median');%ylim([0 1000]);

    figure;         % ANOVA
    [P,ANOVATAB,STATS] = SEAnovaBehav(allProfilesMean);     %needs new function - allows looking at the data in the function
    [c,m,h,nms] = multcompare(STATS);
    
    figure;         % Compare Ephys
    load six_hertz.mat
    create_seq_eff_plot(normalize([-six_hertz nansum(allProfilesMean).'./sum(nFixTimesAllFlies)]),[]); ylabel('Z score'); %title('Behaviour vs Ephys Profile'); 
end
