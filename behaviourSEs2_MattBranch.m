%% behaviour SEs

%As of 09/04/26 - 1:06PM

%To do notes:
%   - Z-score profile medians/etc?


close all; 
clear;
addpath('../Functions');

%% Label prep
n_back = 5;
load('binomial_x_labels_latex_alt_rep.mat','binomial_x_labels_latex');
labels = binomial_x_labels_latex;
ind_horiz = sub2ind(size(labels{1}),1:n_back-1,[ones(1,n_back-2) 5]); %Dynamic
exLabels = [];
for s = 1:size(labels,2)
    exLabels{s} = labels{s}(ind_horiz);
end


%% parameters
limitAngle = 90;             %degrees - Default 90
windowRT = [25];              %degrees - Default 45
refixQualTime = 300;        %loops - Default 300
endTrace = 792;             %loops - Default 792 (Note: This functionally caps the max length of a transient to this value)
limitReQual = 180 %windowRT +0;   %degrees - Default windowRT +0
timeOut = 4000;             %ms -Default 4000

enRT = 1;                   %Enable RT PLots
enMedRT = 0;                %Enable median isomer RT PLots
enQTplot = 1;               %Enable QT vs trial plot
enPlot4x4indiv =0;          %Enable 4x4 histogram plot individual experiments
enPlot4x4super = 1;         %Enable 4x4 histogram plot of super fly

%% Matt parameters
addpath('C:\Users\uqmvan13\ANALYSIS\Matt\Scripts');

poolMode = 1; %Whether to pool (1) all trials of a dataset [for certain plots] or not (0)
    %Note: Due to effort involved in coding, only (currently) applicable to certain (labelled) plots

useNTrialExclusion = 0; %Whether to exclude datasets with less than a specified number of nTrials
if useNTrialExclusion
    nTrialThresh = 100; %Empirical; Based on Scott; Was 120
end

plotAllowance = 2; 
%Whether to allow certain plots (ANOVA, etc); 4 - All [Including individual specific plots], 3 - All, 2 - Major, 1 - Only iterative, 0 - None?

fixTimesMedianSubtraction = 0; %Whether to subtract the median fix time at the fly (1) or endgame (2) level [or not at all (0)]
    %Note: Likely that generating 'negative' RTs will cause issues with one or more plots/analyses

doAutomatedSavePlots = 0; %Whether to use automatedSavePlots to save all figs
if doAutomatedSavePlots
    baseFigPath = 'C:\Users\uqmvan13\ANALYSIS\Scott\FIGOUT';
    alsoSaveOverStruct = 1; %Whether to save overStruct alongside
    saveVectors = 1;
    closeFigures = 1;
end

reOrder = [9,10,11,12,13,14,15,16,1,2,3,4,5,6,7,8];

doBoxplots = 1;

%% import data
%files = dir('./DataLara');
%files = dir('D:\group_vanswinderen\Lara\Control\Data');
%files = dir('C:\Users\lab pc\ANALYSIS\Scott\Data\Scott\Data2');
%files = dir('C:\Users\uqmvan13\ScottCode\DATA_RDM_COPY\Good Data\Experiment Data\SE Data 1per');
%files = dir('C:\Users\uqmvan13\ScottCode\LaraRecapitulation\Split\SD')

%dataPath{1} = 'C:\Users\uqmvan13\ANALYSIS\Scott\Data\Scott\Data2';
dataPath{1} = 'C:\Users\uqmvan13\ANALYSIS\Scott\Data\Lara\Control'
dataPath{2} = 'C:\Users\uqmvan13\ANALYSIS\Scott\Data\Lara\SD'
%dataPath{3} = 'C:\Users\uqmvan13\ANALYSIS\Scott\Data\MegaMan'

%% Prepare for folder iteration

if ~exist('dataPath')
    ['## No folders specified ##']
    crash = yes
end

flagParamSaveList = who;
paramStruct = struct;
for i = 1:size(flagParamSaveList,1)
    eval([ 'paramStruct.',flagParamSaveList{i},' = ',flagParamSaveList{i},';' ])
end
%flagParamSaveList = [flagParamSaveList; 'flagParamSaveList'; 'fol'; 'overStruct'];
flagParamSaveList = [flagParamSaveList; 'flagParamSaveList'; 'windi'; 'fol'; 'overStruct'; 'winStruct'; 'folColours'; 'paramStruct'];

winStruct = struct; %Special structure intended to hold iterative window-related data

for windi = 1:size(windowRT,2)
    clearvars('-except', flagParamSaveList{:})

    windowRTActual = windowRT(windi);
    if windowRTActual < 0
        ['## Alert: Sub-zero windowRT requested  ##'];
        crash = yes
    end
    disp([char(10),'~~~ Requested windowRT: ',num2str(windowRTActual),'degs ~~~',char(10)])

    %winStruct(windi).windowRTActual = windowRTActual;
    winStruct(1).windowRTActual(windi) = windowRTActual; %Logi-centric structure
    winStruct(2).windowRTActual(windi) = windowRTActual;
    
    overStruct = struct; %Designed to hold data for comparing data paths; Not explicitly saved into winStruct
    
    folColours = lines(size(dataPath,2));

    for fol = 1:size(dataPath,2)
        %clearvars('-except', flagParamSaveList{:})
        
        disp([char(10),'Now analysing folder ',num2str(fol),' of ',num2str(size(dataPath,2)),char(10)])
        thisDataPath = dataPath{fol};
        disp([thisDataPath])
        temp = split(thisDataPath,filesep);
        datasetName = temp{end};

        nameStr = [datasetName];
        if size(windowRT,2) > 1
            nameStr = [nameStr, ' - windowRT ',num2str(windowRTActual)];
        end
    
        overStruct(fol).datasetName = datasetName;
        overStruct(fol).scott = [];
        overStruct(fol).matt = []; %Note: matt substructure has additional indexing of 'original' (1) and log10-transformed (2)

        winStruct(1).datasetName{fol} = datasetName; %Logi-centric structure
        winStruct(2).datasetName{fol} = datasetName;
    
    %%

        files = dir(thisDataPath);
        files = files(~[files.isdir],:);

        %Old system; Doesn't work well with multiple folders in combination with many iterative windows
        %if exist('datas') & exist('lastFiles') && isequal(lastFiles,files) %Intended to save time by only loading if different
        %    datas = datas; %Unnecessary but symbolic
        %    lastFiles = files;
        %    disp([char(10),'(Data loading skipped)',char(10)])

        if isfield(winStruct(1),'load') && size(winStruct(1).load,2) >= fol && isequal(winStruct(1).load(fol).files, files)
            datas = winStruct(1).load(fol).datas;
            disp(['(Data loaded from winStruct [Folder #',num2str(fol),'])'])
        else
            datas = cell(length(files),1);
            for i = 1:length(files)
                %datas{i} = table2array(readtable(fullfile('./DataLara',files(i).name)));
                %datas{i} = table2array(readtable(fullfile('C:\Users\lab pc\ANALYSIS\Scott\Data\Scott\Data2',files(i).name)));
                %datas{i} = table2array(readtable(fullfile('C:\Users\uqmvan13\ScottCode\DATA_RDM_COPY\Good Data\Experiment Data\SE Data 1per',files(i).name)));
                %datas{i} = table2array(readtable(fullfile('C:\Users\uqmvan13\ScottCode\LaraRecapitulation\Split\SD',files(i).name)));
                datas{i} = table2array(readtable(fullfile(thisDataPath,files(i).name)));
                disp(files(i).name);
            end
            %lastFiles = files;
            winStruct(1).load(fol).files = files;
            winStruct(1).load(fol).datas = datas;
        end
        
        nFlies = length(datas);
        disp(nFlies);
        
        %QA
        if nFlies == 0
            ['-# Alert: No data detected #-']
            crash = yes
        end
        
        %% collate and organise data
        FLIES = struct;
        
        for fly = 1: nFlies
        
        thisFlyData = datas{fly};

        FLIES(fly).filename = files(fly).name;
        
        %there are some zero rows for some reason
        % data = data(logical(sum(data,2)),:);
        
        % data = data(data(:,10) > 0,:);
        
        %starting points of inter-stimulus periods (careful this is actually the index just before the start)
        % last trial is usually interrupted so get rid of last four (this may fail is experiment is interrupted during open loop)
        interStimulusStart = find(diff(thisFlyData(:,4) == 200) == 1); %interStimulusStart = interStimulusStart(1:end-4);
        interStimulusEnd = find(diff(thisFlyData(:,4) == 200) == -1); %interStimulusEnd = interStimulusEnd(1:end-4);
        %QA
        if numel(interStimulusStart) ~= numel(interStimulusEnd)
            if windi == 1 %Only display on first
            disp(['-# Caution: Disparate number of stim starts/ends; Equalising to mod 4 #-'])
            end
            %interStimulusStart = interStimulusStart( 1: );
            %interStimulusEnd = interStimulusEnd( 1:min([numel(interStimulusStart),numel(interStimulusEnd)]) );
                %Assumes no phase-shift or anything weirder
               
            temp = [0:4:min([numel(interStimulusStart),numel(interStimulusEnd)])]; %By inference, curren number of elements is not mod 4, iterating by 4 will give the lowest mod
            interStimulusStart = interStimulusStart( 1:temp(end) );
            interStimulusEnd = interStimulusEnd( 1:temp(end) );
                %Again, assumes no phase-shift
        end
        
        nTrials = (length(interStimulusStart)/4);
        
        if useNTrialExclusion
            if nTrials < nTrialThresh
                disp(['-# Dataset #',num2str(fly),' lacks sufficient nTrials (',num2str(nTrials),' < ',num2str(nTrialThresh),') #-'])
                continue
            end
        end
        
        %fluid
        
        %Testatory for stimulus sync
        %{
        figure
        subplot(2,1,1)
        plot(interStimulusStart)
        hold on
        plot(interStimulusEnd)
        subplot(2,1,2)
        plot(interStimulusEnd - interStimulusStart)
        %}
        
        
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
        
        %%
        
        %Check to see if any FLIES elements empty (Due to exclusion)
        excCount = 0;
        for fly = size(FLIES,2):-1:1
            if isempty( FLIES(fly).nTrials )
                FLIES(fly) = [];
                excCount = excCount + 1;
            end
        end
        nFlies = size(FLIES,2); %This is necessary to prevent issues with later analysis
        if excCount ~= 0
            disp([num2str(excCount),' datasets excluded; ',num2str(size(FLIES,2)),' remain'])
        end
    
        overStruct(fol).N = size(FLIES,2);
        overStruct(fol).files = files;
        
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
        
        if(enQTplot) && plotAllowance >= 3    % switch at top of file
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
            set(gcf,'Name',nameStr);
    
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
                %thisTrace = thisFly.barPosition(thisFly.closedLoopStart(trial):thisFly.closedLoopStart(trial)+ endTrace - 1); %Liable to crash if trial near end
                
                coords = thisFly.closedLoopStart(trial):thisFly.closedLoopStart(trial)+ endTrace - 1;
        
                if nanmax(coords) <= size(thisFly.barPosition,1) %In bounds
                    thisTrace = thisFly.barPosition(coords);
                    tracesSEQ(:,seq,trial) = thisTrace;
                else %Trial off edge of recording
                    thisTrace = nan( size(coords,2), 1 );
                    tracesSEQ(:,seq,trial) = thisTrace;
                    %badTrials(trial) = 1; %Borrow system
                        %Unnecessary, as should be automatically caught below
                end
                
                % if the bar close to 0? SCOTT: CHECK WINDOW BOUNDS HERE
                %window = [-windowRT windowRT]; %Original, single-shot specification
                window = [-windowRTActual windowRTActual];

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

            %stop
            
            % get rid of bad trials here so no need to do it downstream
            fixTimes = fixTimes(~badTrials,:);
            tracesSEQ = tracesSEQ(:,:,~badTrials);

            %Normalise fixation times by individual, if requested
            if fixTimesMedianSubtraction == 1
                %distract
                FLIES(fly).fixTimes32Original = fixTimes; %Since this variable will be modified momentarily, save this pre-subtraction version (Mostly for log sake)

                fixTimes( fixTimes == 0 ) = NaN; %Otherwise multiple critical errors like 0->neg values etc
                fixTimes = fixTimes - nanmedian( fixTimes, [1,2] ); %Subtract median within individual, across sequence
                    %A case could be made for subtracting median within sequence, but that would make comparison harder by definition
                fixTimes( isnan(fixTimes) ) = 0; %Return back to zero-containing state
            end
            
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

            %Some post plots, courtesy of Matt
            if plotAllowance >= 4
                %All transients
                temp = squeeze( nansum(FLIES(fly).tracesSEQ, [2]) ); %Remove sparsity
                figure
                hold on
                for i = 1:size(temp,2)
                    plot(temp(:,i))
                end
                xlabel('Frame')
                ylabel('Position (degs)')
                title(['All transients - Dataset #',num2str(fly),' [n=',num2str(size(temp,2)),']'])
                set(gcf,'Name',nameStr);
            end
           
        end
        
        %% Plot discarded trials vs Sequence
        if plotAllowance >= 3
        figure;                 
        create_seq_eff_plot([sum([FLIES.RTB],2) sum([FLIES.timeOuts],2) sum([FLIES.reQual],2)],[]);
        legend('MaxAngle Fail','Timeout Fails','Requal Fails','Location','northeast','FontSize',8); 
        title('Discarded Trials vs Stimulus Sequence');
        set(gcf,'Name',nameStr);
        end
        
        %% make 4x4 histogram plot for each experiment
        if(enPlot4x4indiv) && plotAllowance >= 2
            load('binomial_x_labels_latex_alt_rep.mat','binomial_x_labels_latex');
        
            for fly = 1:length(FLIES)
                figure;
                for i = 1:16
                   subplot(4,4,i); histogram(FLIES(fly).fixTimes(logical(FLIES(fly).fixTimes(:,i)),i),10); 
                   title(binomial_x_labels_latex{seq}(ind_horiz));
                end
                set(gcf,'Name',nameStr); %Better to use indiv fly names?
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
        if(enMedRT) && plotAllowance >= 2
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
            set(gcf,'Name',nameStr);
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
                if ~fixTimesMedianSubtraction
                fixTimesAllFlies32{fly} = FLIES(fly).fixTimes32;
                else
                fixTimesAllFlies32{fly} = FLIES(fly).fixTimes32Original;
                end
            end
        
            fixTimesAllFlies = cell2mat(fixTimesAllFlies);
            fixTimesAllFlies32 = cell2mat(fixTimesAllFlies32);
            
            RAllFlies = calculateSEsBehav(fixTimesAllFlies);
        
            if plotAllowance >= 2
        %    figure; create_seq_eff_plot(RAllFlies.meanFixTimes.',[],'errors',RAllFlies.semFixTimes.');
            figure; 
            create_seq_eff_plot(RAllFlies.medianFixTimes.',[],'errors',RAllFlies.semFixTimes.');
            
            title('SE super Median');
            set(gcf,'Name',nameStr);

            figure; 
            create_seq_eff_plot(RAllFlies.meanFixTimes.',[],'errors',RAllFlies.semFixTimes.');
            title('SE super Mean');
            set(gcf,'Name',nameStr);

            if exist('reOrder') && ~isempty(reOrder)
            figure; 
            create_seq_eff_plot(RAllFlies.medianFixTimes.',[],'errors',RAllFlies.semFixTimes.', 'reOrder',reOrder); %Note: Relies on updated version of create_seq_eff_plot
            title('SE super Median [reOrder]');
            set(gcf,'Name',nameStr);

            figure; 
            create_seq_eff_plot(RAllFlies.meanFixTimes.',[],'errors',RAllFlies.semFixTimes.', 'reOrder',reOrder); %Note: Relies on updated version of create_seq_eff_plot
            title('SE super Mean [reOrder]');
            set(gcf,'Name',nameStr);
            end

            %to do: reorder isomers
        
            plotIsomersBehav(fixTimesAllFlies32);

            if exist('reOrder') && ~isempty(reOrder)
            tempOrder = 32-flip(reOrder)+1; %Create flipped version of reOrder (Only 75% confident correct)
            temp = fixTimesAllFlies32; %Note: Not yet Dinis ordered
            temp(:,1:16) = temp(:,seq_eff_order(5)); %Order L isomer
            temp(:,17:32) = temp(:,32-flip(seq_eff_order(5))+1); %Be careful to remember to flip seq_eff_order
            temp( temp == 0 ) = NaN;
            tempMedian = nanmedian( temp,1 );
            tempSEM = nanstd(temp,[],1) ./ sqrt( nansum(~isnan(temp),1) );
            %plotIsomersBehav(fixTimesAllFlies32);
            figure; 
            create_seq_eff_plot(tempMedian(1:16).',[],'errors',tempSEM(1:16).', 'reOrder',reOrder); %Plot L(?) isomer
            create_seq_eff_plot(tempMedian(17:32).',[],'errors',tempSEM(17:32).', 'reOrder',flip(reOrder)); %Other isomer
            title('Bootleg isomer plot [reOrder]');
            set(gcf,'Name',nameStr);

            end
            end
        
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
                 disp(['Log L vs R T-test: h ',num2str(h),' p ',num2str(p)]);
                meanLeftLog = nanmean(RT_Hist_Array_Left_log); 
                medianLeftLog = nanmedian(RT_Hist_Array_Left_log);
                meanRightLog = nanmean(RT_Hist_Array_Right_log);
                medianRightLog = nanmedian(RT_Hist_Array_Right_log);
                
        
                 if plotAllowance >= 2
                 figure; histogram(RT_Hist_Array,'NumBins', 25); title('All RT');  xlabel('Time (ms)'); xlim([0,4200]); set(gcf,'Name',nameStr);
                 end
                 if plotAllowance >= 3
                 figure; histogram(RT_Hist_Array_Left,'NumBins', 25,'facecolor',[0.4609 0.7070 0.7695]); title('Left RT');  xlabel('Time (ms)'); xlim([0,4200]); set(gcf,'Name',nameStr);
                 figure; histogram(RT_Hist_Array_Right,'NumBins', 25,'facecolor',[0.9140 0.7109 0.4609]); title('Right RT');  xlabel('Time (ms)'); xlim([0,4200]); set(gcf,'Name',nameStr);
                 end

                 if plotAllowance >= 2
                 figure('color','w','position',[10 10 1200 400]); set(gcf,'Name',nameStr);
                 subplot(1,2,1);histogram(RT_Hist_Array_Left_log,'NumBins', 25,'facecolor',[0.4609 0.7070 0.7695]); title('Left RT');  xlabel('Time (Log sec)'); ylim([0,90]); ylabel('Count'); xlim([2,4]);
                 line([meanLeftLog,meanLeftLog],[0,90],'linewidth',2,'color','k','LineStyle','-'); text(3.3,85,['mean=',num2str(meanLeftLog,3)],'FontSize',14);
                 line([medianLeftLog,medianLeftLog],[0,90],'linewidth',2,'color','r','LineStyle','--'); text(3.3,78,['median=',num2str(medianLeftLog,3)],'FontSize',14);
                 legend('Left RT','mean','median','Location','northwest'); 
                 
                 subplot(1,2,2); histogram(RT_Hist_Array_Right_log,'NumBins', 25,'facecolor',[0.9140 0.7109 0.4609]); title('Right RT');  xlabel('Time (Log sec)'); ylim([0,90]); ylabel('Count'); xlim([2,4]);
                 line([meanRightLog,meanRightLog],[0,90],'linewidth',2,'color','k','LineStyle','-'); text(3.3,85,['mean=',num2str(meanRightLog,3)],'FontSize',14);
                 line([medianRightLog,medianRightLog],[0,90],'linewidth',2,'color','r','LineStyle','--'); text(3.3,78,['median=',num2str(medianRightLog,3)],'FontSize',14);
                 legend('Right RT','mean','median','Location','northwest'); 
                 end
    
                 overStruct(fol).scott.RT_Hist_Array = RT_Hist_Array;
                 overStruct(fol).scott.RT_Hist_Array_Left = RT_Hist_Array_Left;
                 overStruct(fol).scott.RT_Hist_Array_Right = RT_Hist_Array_Right;
                
            end
            
        
            if(enPlot4x4super) && plotAllowance >= 2
                load('binomial_x_labels_latex_alt_rep.mat','binomial_x_labels_latex');
                ind_horiz = sub2ind(size(binomial_x_labels_latex{1}),1:4,[1 1 1 5]);
                figure; 
                for seq = 1:16
                   subplot(4,4,seq); histogram(fixTimesAllFlies(logical(fixTimesAllFlies(:,seq)),seq),10,'NumBins', 16); xlim([0,4000]);ylim([0,40]);
                   title(binomial_x_labels_latex{seq}(ind_horiz));
                end
                set(gcf,'Name',nameStr);
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
        
            if plotAllowance >= 3
            figure; create_seq_eff_plot(nansum(allProfilesMean).'./sum(nFixTimesAllFlies),[],'errors',semFixTimesAllFlies.'); title('SE weighted Mean'); set(gcf,'Name',nameStr); %ylim([0 1000]);
            figure; create_seq_eff_plot(nansum(allProfilesMedian).'./sum(nFixTimesAllFlies),[],'errors',semFixTimesAllFlies.');title('SE weighted Median'); set(gcf,'Name',nameStr); %ylim([0 1000]);
        
            figure;         % ANOVA
            [P,ANOVATAB,STATS] = SEAnovaBehav(allProfilesMean);     %needs new function - allows looking at the data in the function
            [c,m,h,nms] = multcompare(STATS);
            set(gcf,'Name',nameStr);
            
            figure;         % Compare Ephys
            load six_hertz.mat
            create_seq_eff_plot(normalize([-six_hertz nansum(allProfilesMean).'./sum(nFixTimesAllFlies)]),[]); ylabel('Z score'); %title('Behaviour vs Ephys Profile'); 
            set(gcf,'Name',nameStr);
            end
        
            %% Matt specific sequence comparisons
            if plotAllowance >= 3
            %Pre plot with all (medians)
            figure
            for fly = 1:size(allProfilesMedian,1)
                plot(allProfilesMedian(fly,:),'-o')
                hold on
            end
            hold off
            xlim([0,size(allProfilesMedian,2)+1])
            xticks([1:size(allProfilesMedian,2)])
            xticklabels(exLabels)
            xtickangle(270)
            title(['Individual fly medians'])
            set(gcf,'Name',nameStr);
            end
        
            allProfilesMedian_meanCorr = allProfilesMedian - nanmean(allProfilesMedian,2);
        
            %Repeat with mean-corrected medians
            if plotAllowance >= 3
            figure
            for fly = 1:size(allProfilesMedian_meanCorr,1)
                plot(allProfilesMedian_meanCorr(fly,:),'-o')
                hold on
            end
            hold off
            xlim([0,size(allProfilesMedian_meanCorr,2)+1])
            xticks([1:size(allProfilesMedian_meanCorr,2)])
            xticklabels(exLabels)
            xtickangle(270)
            title(['Individual fly medians [Mean-corrected]'])
            set(gcf,'Name',nameStr);
            end
        
        
            if plotAllowance >= 2
            %Comp plot
            useCorrected = 1; %Whether to use mean-corrected median data
            seqCompare = [1,9;...
                          16,8]; %Based on canonical ordering
        
            for seqI = 1:size(seqCompare,1)
                thisComp1 = seqCompare(seqI,1);
                thisComp2 = seqCompare(seqI,2);
                thisCompSeq{1} = exLabels{thisComp1};
                thisCompSeq{2} = exLabels{thisComp2};
                compData = [];
                if ~useCorrected
                    compData = [allProfilesMedian(:,thisComp1),allProfilesMedian(:,thisComp2)];
                else
                    compData = [allProfilesMedian_meanCorr(:,thisComp1),allProfilesMedian_meanCorr(:,thisComp2)]; %Use mean-corrected data
                end
                [h,p] = ttest2(compData(:,1),compData(:,2));
                %compData = nansum(compData,1); %Replicating above
                %compData = compData ./ sum(nFixTimesAllFlies(:,[thisComp1,thisComp2])); %Note: This divides the SEs (median) by the total n ?
                   %The above line only really applicable if using sum        
        
                compMean = nanmean(compData,1);
                compSEM = nanstd(compData,[],1) / sqrt(size(compData,1)); %Mean of median?
                    %These lines applicable if not copying Scott system for preparing SEs plot
        
                %Plot
                figure
                errorbar(compMean,compSEM,'LineWidth',3)
                %plot([1,2],compData,'o')
                xlim([0,3])
                xticks([1,2])
                xticklabels({exLabels{thisComp1},exLabels{thisComp2}})
                xtickangle(270)
                titleStr = [exLabels{thisComp1},' vs ',exLabels{thisComp2},' - mean of per-fly medians [p=',num2str(p,3),']'];
                if useCorrected
                    titleStr = [titleStr,'[Mean-corrected]'];
                end
                title(titleStr)
                set(gcf,'Name',nameStr);
        
            end
            end
        
            %Matt meta stats
            
            %Fixation hists
            binEdges = linspace(-180,180,64);
        
            NColl = nan(1,size(binEdges,2)-1);
            allData = [];
            for fly = 1:size(FLIES,2)
                thisData = FLIES(fly).barPosition;
                thisData( thisData == 60 ) = NaN; %Remove OL periods
                thisData( thisData == -60 ) = NaN;
                allData = [allData; thisData]; %Horrible inefficiency of non-pre-allocating
                [N,~] = histcounts(thisData, binEdges, 'Normalization','pdf');
                NColl(fly,:) = N; %This pseudo-guarantees correct sizing, because it will crash otherwise
            end
        
            %edges = edges(2:end) - (edges(2)-edges(1))/2;
            acEdges = binEdges(2:end);
            NMean = nanmean(NColl,1);
            NMedian = nanmedian(NColl,1);
            NSEM = nanstd(NColl,[],1) / sqrt(size(NColl,1));
        
            if plotAllowance >= 3
            figure
            errorbar(acEdges,NMean,NSEM,'LineWidth',2,'Color','k')
            hold on
            errorbar(acEdges,NMedian,NSEM,'LineWidth',2,'Color','r')
            for fly = 1:size(NColl,1)
                plot(acEdges, NColl(fly,:), ':o');
            end
            xlim([binEdges(1),binEdges(end)])
            xlabel('Bar position (degs)')
            ylabel('Occupancy (prop.)')
            title('SEs bar pos. occupancy hist.')
            set(gcf,'Name',nameStr);
            %Save onion details
            overStruct(fol).plots.onionHist.edges = acEdges;
            overStruct(fol).plots.onionHist.mean = NMean;
            overStruct(fol).plots.onionHist.median = NMedian;
            overStruct(fol).plots.onionHist.SEM = NSEM;
            temp = getframe;
            overStruct(fol).plots.onionHist.frame = temp.cdata;
            overStruct(fol).plots.onionHist.files = files;
            overStruct(fol).plots.onionHist.params = paramStruct;


            %Sum hist
            figure
            [histi] = histogram(allData, 'BinEdges',acEdges);
            xlabel('Bar position (degs)')
            ylabel('Occupancy (sum)')
            title('SEs bar pos. occupancy hist. (Sum)')
            set(gcf,'Name',nameStr);
            %Save sum details
            overStruct(fol).plots.sumHist.edges = histi.BinEdges;
            overStruct(fol).plots.sumHist.values = histi.Values;
            temp = getframe;
            overStruct(fol).plots.sumHist.frame = temp.cdata;
            overStruct(fol).plots.sumHist.files = files;
            overStruct(fol).plots.sumHist.params = paramStruct;

            end
            %shikari
        
            %Meta figure
            if plotAllowance >= 2
            figure
            hold on
            for fly = 1:size(FLIES,2)
                temp = [FLIES(fly).nTrials, nansum(FLIES(fly).RTB), nansum(FLIES(fly).reQual),nansum(FLIES(fly).timeOuts), FLIES(fly).goodTrial];
                plot(temp, '--o')
            end
            %legend
            xticks([1:5])
            xlim([0,6])
            xticklabels([{'nTrials'},{'RTB'},{'reQual'},{'timeOuts'},{'goodTrials'}])
            title('Meta stats')
            set(gcf,'Name',nameStr);
            end
                
            %Matt barplot/s
        
            for logi = 0:1 %Whether to be logMode or not
                if logi == 1 && fixTimesMedianSubtraction
                    continue %Do not do log if median subtracted
                        %Note: This system for altering for loop may not work
                end
    
                overStruct(fol).matt(logi+1).logMode = logi;
            
                %Prepare an appropriate data structure/s
                %omega32Times = nan(nanmax([FLIES.goodTrial]), 32, size(FLIES,2));
                omega32Times = zeros(nanmax([FLIES.goodTrial]), 32, size(FLIES,2));
                    %Note: As assembled, omega32Times will have a row dimension at the max n of events of any one individual fly
                    %An individual fly will only occupy FLIES(fly).goodTrial rows
                for fly = 1:size(FLIES,2)
                    temp = FLIES(fly).fixTimes32;
                    omega32Times( 1:size(temp,1) , : , fly ) = temp;
                end
                %omega32Times( isnan(omega32Times) ) = 0; %STRANGE METHODS (See below secondary note)

                %REORDER STRUCTURE, SINCE THIS NOT DONE AT FLY ACQUISITION LEVEL
                    %See nStructIndexFixing notes for minituae
                    %In essence, og seq 16 -> ordered seq 2, whilst og seq 24 -> ordered seq 29 (mirror of ordered seq 4, if flipped/collapsed)
                omega32Times(:, 1:16, :) = omega32Times(:, seq_eff_order(5), :); %Stage 1
                omega32Times( : , 17:32, : ) = omega32Times( :, 32-flip(seq_eff_order(5))+1, : ); %Stage 2
                    %Note: As reordered, omega32Times is now cols 1:16 = Dinis classic (5back), 17:32 = mirror (i.e. 16==17, 15==18, etc)
                    %Write nStruct inds accordingly
                %Test cases, more to understand logic than actually to do anything
                %{
                om1 = omega32Times(:,:,1); %Grab first layer (Fly) of omega32Times
                om1 = om1 + flip(om1,2); %Classic Dinis flip addition
                om1 = om1(:,1:16); %Classic Dinis subselection
                nRows = FLIES(1).goodTrial; %Need to find, because omega32Times pre-assembled for max potential event #
                isequal( om1(1:nRows,:) , fixTimesAllFlies(1:nRows,:) ) %If not 1, then error has occurred
                    %Note: Obviously won't work if logi = log
                %}
                %grant

                omega32Times( omega32Times == 0 ) = NaN; %Do not allow zeroes
                    %These are present because fixTimes32 is a sparse array
                    %Note: Mathematically confirmed to be more or less identical to fixTimesAllFlies32, which is used for Dinis SEs
                    %(following rearrangement to tile 3rd dimension as rows)
                    %Secondary note: Must be done after reordering, because NaNs do not addition nicely with values
        
                if logi == 1
                    omega32Times = log10(omega32Times);
                end
    
                overStruct(fol).matt(logi+1).omega32Times = omega32Times;
                overStruct(fol).matt(logi+1).omega32Times_isOrdered = 1; %"Never Forget"
            
                %flyColours = parula(size(omega32Times,3));
            
                %Overall median
                %jude
                if ~poolMode %Calculate median/SEM of each individual
                    flyColours = parula(size(omega32Times,3));
                    overMedian = squeeze( nanmedian(omega32Times , [1,2]) );
                    overSEM = squeeze( nanstd( omega32Times(:,:,:), [], [1,2] ) )  ./ sqrt( squeeze( nansum( ~isnan(omega32Times(:,:,:)) , [1,2] ) ) );
                    %winStruct(windi).logi(logi+1).thisMedian(fol) = nanmedian(overMedian); %Use of fol allows cross-dataset comparison later
                    %winStruct(windi).logi(logi+1).thisSEM(fol) = nanstd(overMedian,[],1) / sqrt(size(overMedian,1));
                else %Calculate median/SEM of *all* events
                    %flyColours = parula(1);
                    flyColours = winter(1);
                    overMedian = squeeze( nanmedian(omega32Times , 'all') ); %Empirically confirmed same as reshaping all to one col (or row) and applying median
                    overSEM = squeeze( nanstd( omega32Times, [], 'all' ) )  ./ sqrt( squeeze( nansum( ~isnan(omega32Times) , 'all' ) ) ); %Note: n ballpark somewhere in the ~2000s for a good-sized dataset
                end
                winStruct(logi+1).thisMedian(windi,fol) = nanmedian(overMedian); %Use of fol allows cross-dataset comparison later
                winStruct(logi+1).thisSEM(windi,fol) = nanstd(overMedian,[],1) / sqrt(size(overMedian,1));
            
                if plotAllowance >= 2
                %Overall jenga plot
                figure
                scatter(overMedian,[1:size(overMedian,1)], [], flyColours)
                hold on
                for fly = 1:size(overMedian,1)
                    line([overMedian(fly)-0.5*overSEM(fly),overMedian(fly)+0.5*overSEM(fly)], [fly,fly],'LineStyle','--','Color',flyColours(fly,:))
                end
                xlabel('RT (ms)')
                if ~poolMode
                    ylabel('Dataset #')
                    yticks([1:size(overMedian,1)])
                else
                    ylabel('Pooled data')
                    yticks([])
                end
                
                title('Overall median RT')
                %And barplot, for some reason
                figure
                thisMedian = nanmedian(overMedian);
                thisSEM = nanstd(overMedian,[],1) / sqrt(size(overMedian,1));
                bar([1],[thisMedian])
                hold on
                errorbar([thisMedian],[thisSEM])
                xlabel('Overall median RT')
                %xticklabels({'Left','Right'})
                xticklabels([])
                ylabel('Median of median RTs (ms)')
                %title('Overall median RT barplot')
                if ~poolMode
                    titleStr = ['Overall median RT barplot'];
                else
                    titleStr = ['Pooled RT barplot'];
                end
                if logi
                    titleStr = [titleStr,' (Log)'];
                end
                title(titleStr)
                set(gcf,'Name',nameStr);
                end
            
                
                %L vs R
                if ~poolMode
                    LDataMean = squeeze( nanmean( omega32Times(:,1:16,:), [1,2] ) ); %Average across trial and seq to give per-fly mean RT for L
                    LDataMedian = squeeze( nanmedian( omega32Times(:,1:16,:), [1,2] ) ); %Ditto, median
                    LDataSEM = squeeze( nanstd( omega32Times(:,1:16,:), [], [1,2] ) )  ./ sqrt( squeeze( nansum( ~isnan(omega32Times(:,1:16,:)) , [1,2] ) ) ); %Per-dataset SEM for L ending fix times
                
                    RDataMean = squeeze( nanmean( omega32Times(:,17:32,:), [1,2] ) ); 
                    RDataMedian = squeeze( nanmedian( omega32Times(:,17:32,:), [1,2] ) );
                    RDataSEM = squeeze( nanstd( omega32Times(:,17:32,:), [], [1,2] ) )  ./ sqrt( squeeze( nansum( ~isnan(omega32Times(:,17:32,:)) , [1,2] ) ) ); 

                else
                    LDataMean = squeeze( nanmean( omega32Times(:,1:16,:), 'all' ) ); %Average across all
                    LDataMedian = squeeze( nanmedian( omega32Times(:,1:16,:), 'all' ) ); %Ditto, median
                    LDataSEM = squeeze( nanstd( omega32Times(:,1:16,:), [], 'all' ) )  ./ sqrt( squeeze( nansum( ~isnan(omega32Times(:,1:16,:)) , 'all' ) ) ); %Per-dataset SEM for L ending fix times
                    LDataAll = reshape( omega32Times(:,1:16,:), [size(omega32Times,1)*16*size(omega32Times,3)], 1 );
                    LDataAll( isnan(LDataAll) ) = [];

                    RDataMean = squeeze( nanmean( omega32Times(:,17:32,:), 'all' ) ); 
                    RDataMedian = squeeze( nanmedian( omega32Times(:,17:32,:), 'all' ) );
                    RDataSEM = squeeze( nanstd( omega32Times(:,17:32,:), [], 'all' ) )  ./ sqrt( squeeze( nansum( ~isnan(omega32Times(:,17:32,:)) , 'all' ) ) ); 
                    RDataAll = reshape( omega32Times(:,17:32,:), [size(omega32Times,1)*16*size(omega32Times,3)], 1 );
                    RDataAll( isnan(RDataAll) ) = [];
                
                end

                %stronger

                %'Jenga' plot
                if plotAllowance >= 3
                figure
            
                for sidi = 1:2
                    subplot(1,2,sidi)
                        hold on
                    if sidi == 1
                        thisMean = LDataMean;
                        thisMedian = LDataMedian;
                        thisSEM = LDataSEM;
                        titleStr = ['L ending seqs median (and mean)'];
                        if logi
                            titleStr = [titleStr,' (Log)'];
                        end
                        title(titleStr)
                    else
                        thisMean = RDataMean;
                        thisMedian = RDataMedian;
                        thisSEM = RDataSEM;
                        titleStr = ['R ending seqs median (and mean)'];
                        if logi
                            titleStr = [titleStr,' (Log)'];
                        end
                        title(titleStr)
                    end
            
                    %for fly = 1:size(thisMean,1)
                    scatter(thisMean,[1:size(thisMean,1)], [], flyColours)
                    %end
                    scatter(thisMedian,[1:size(thisMedian,1)], [], flyColours, 'Marker', '+')
                    for fly = 1:size(thisMedian,1)
                        line([thisMedian(fly)-0.5*thisSEM(fly),thisMedian(fly)+0.5*thisSEM(fly)], [fly,fly],'LineStyle','--','Color',flyColours(fly,:))
                    end
            
                    xlabel('RT (ms)')
                    if ~poolMode
                        ylabel('Dataset #')
                        yticks([1:size(overMedian,1)])
                    else
                        ylabel('Pooled data')
                        yticks([])
                    end
                    %yticks([1:size(overMedian,1)])
            
                end
                set(gcf,'Name',nameStr);
                end
            
                %Supermedian (L vs R)
                superMedian = [];
                superSEM = [];
                if ~poolMode
                    superMedian{1} = nanmedian(LDataMedian,1);
                    superSEM{1} = nanstd(LDataMedian,[],1) / sqrt(size(LDataMedian,1));
                    superMedian{2} = nanmedian(RDataMedian,1);
                    superSEM{2} = nanstd(RDataMedian,[],1) / sqrt(size(RDataMedian,1));
                    superData{1} = LDataMedian;
                    superData{2} = RDataMedian;

                else
                    superMedian{1} = LDataMedian; %Use pooled data from above, since functionally identical already to this operation
                    superSEM{1} = LDataSEM;
                    superMedian{2} = RDataMedian;
                    superSEM{2} = RDataSEM;
                    superData{1} = LDataAll;
                    superData{2} = RDataAll;

                end
            
                if plotAllowance >= 2
                [H,p] = ttest2(superData{1},superData{2});

                figure
                %bar([1,2],[superMedian{1},superMedian{2}])
                hold on
                errorbar([superMedian{1},superMedian{2}],[superSEM{1},superSEM{2}],'Color',folColours(fol,:))
                xlim([0,3])
                xticks([1,2])
                xlabel('Ending side')
                xticklabels({'Left','Right'})
                if ~poolMode
                    ylabel('Median of median RTs (ms)')
                    titleStr = ['L vs R ending overall medians'];
                else
                    ylabel('Median pooled RTs (ms)')
                    titleStr = ['L vs R ending pooled medians'];
                end
                titleStr = [titleStr,' p=',num2str(p), ' [L n: ',num2str(numel(superData{1})),', R n: ',num2str(numel(superData{1})),']'];
                if logi
                    titleStr = [titleStr,' (Log)'];
                end
                title(titleStr)
                set(gcf,'Name',nameStr);
                end
    
                %overStruct(fol).matt(logi+1).LDataMedian = LDataMedian;
                %overStruct(fol).matt(logi+1).RDataMedian = RDataMedian;
                overStruct(fol).matt(logi+1).superMedian = superMedian;
                overStruct(fol).matt(logi+1).superData = superData;
            
            
                %Various nBack plots (Bootleg style)
                nStruct = struct;
                nStruct(2).inds = [1:8;9:16];
                %nStruct(2).inds = [nStruct(2).inds, nStruct(2).inds+16]; %Lazy easy way of mirroring indices across 32 full-isomer seqs; PATENTLY INCORRECT
                nStruct(2).inds = [nStruct(2).inds,32-nStruct(2).inds+1]; %Actually (hopefully) correct way of *flipping* (not mirroring) indices (e.g. 1==32, 2==31, ..., 16==17, etc)
                nStruct(2).name = {'R','A'};

                %{{
                nStruct(3).inds = [1:4;  5:8;  9:12;  13:16]; 
                %nStruct(3).inds = [nStruct(3).inds, nStruct(3).inds+16]; 
                nStruct(3).inds = [nStruct(3).inds,32-nStruct(3).inds+1];
                nStruct(3).name = {'RR','AR', 'RA', 'AA'};

                nStruct(4).inds = [1,2;  3,4;  5,6;  7,8; 9,10; 11,12; 13,14; 15,16]; 
                %nStruct(4).inds = [nStruct(4).inds, nStruct(4).inds+16]; 
                nStruct(4).inds = [nStruct(4).inds,32-nStruct(4).inds+1];
                nStruct(4).name = {'RRR','ARR', 'RAR', 'AAR','RRA','ARA','RAA','AAA'};

                %nStruct(5).inds = [1:16]'; 
                %nStruct(5).inds = [nStruct(5).inds, nStruct(5).inds+16]; 
                %nStruct(5).name = exLabels;
                %temp = seq_eff_order(6); %Note: No 6back sequence does *not* seem to be same as 1:16 w/ isomers (Unless natural vs canonical taken into account, plus interleaving?)
                %                nStruct(5).inds = temp(1:16);  
                nStruct(5).inds = [1:16]'; 
                %nStruct(5).inds = [nStruct(5).inds,temp(17:32)]; 
                nStruct(5).inds = [nStruct(5).inds,32-nStruct(5).inds+1];
                nStruct(5).name = exLabels;
                %}

                %haywyre

                %Duplicate nStruct, but separated into isomers (even bootylegger)
                nIsoStruct = struct; %1 - Isomer 1, etc
                for enny = 1:size(nStruct,2)
                    if ~isempty(nStruct(enny).inds)
                        indWidth = size(nStruct(enny).inds,2);
                        if rem(indWidth,2) ~= 0
                            ['-# Alert: Non-symmetric indices potentially present #-']
                            continue
                        end
                        if numel(nStruct(enny).inds) ~= 32 %Assumes both no overlap in indices, and capture of all
                            ['-# Alert: Not all sequences present in indices potentially #-']
                            continue
                        end

                        nIsoStruct.nBack(enny).nStruct(1).name = nStruct(enny).name;
                        nIsoStruct.nBack(enny).nStruct(2).name = nStruct(enny).name;
                        
                        nIsoStruct.nBack(enny).nStruct(1).inds = nStruct(enny).inds(:,1:0.5*indWidth);
                        nIsoStruct.nBack(enny).nStruct(2).inds = nStruct(enny).inds(:,0.5*indWidth+1:end);
                            %Note: **Critically** relies on flip present in indices
                    end
                end

                %Reporter
                if plotAllowance >= 4

                    for enny = 1:size(nStruct,2)
                        disp([char(10),'nBack: ',num2str(enny)])   
                        disp(['Indices:'])
                        disp([num2str(nStruct(enny).inds)]) 
                        disp(['Names:'])
                        disp([nStruct(enny).name]) 
                    end

                    for naisu = 1:2 %Note: Different usage compared to below
                        disp(['Isomer ',num2str(naisu)])
                        for enny = 1:size(nStruct,2)
                            disp([char(10),'nBack: ',num2str(enny)])
                            if ~isempty(nIsoStruct.nBack(enny).nStruct)
                            disp(['Indices:'])
                            disp([num2str(nIsoStruct.nBack(enny).nStruct(naisu).inds)]) 
                            disp(['Names:'])
                            disp([nIsoStruct.nBack(enny).nStruct(naisu).name]) 
                            end
                        end
                    end
                    
                end

            
                for naisu = 1:3 %1 - Merged isomers, 2 - Isomer 1, 3 - Isomer 2
                    for enny = 1:size(nStruct,2) %Shared across naisu
                        if ~isempty(nStruct(enny).inds)
                            %Choose what doing
                            if naisu == 1
                                thisInds = nStruct(enny).inds;
                                thisNames = nStruct(enny).name;
                            elseif naisu == 2
                                thisInds = nIsoStruct.nBack(enny).nStruct(1).inds;
                                thisNames = nIsoStruct.nBack(enny).nStruct(1).name;
                            elseif naisu == 3
                                thisInds = flip( nIsoStruct.nBack(enny).nStruct(2).inds ,2 ); %Flip, to maintain correctness with labels/etc; NOTE CRITICALITY OF DIMENSION SPECIFICATION
                                thisNames = nIsoStruct.nBack(enny).nStruct(2).name;
                            end

                            %thisInds = nStruct(enny).inds; %Old
                            thisData = []; %'seq' x fly
                            for subEn = 1:size(thisInds,1)
                                thisSubInds = thisInds(subEn,:);
                                if ~poolMode
                                    thisData(subEn,:) = squeeze( nanmedian( omega32Times(:,thisSubInds,:) , [1,2] ) )';
                                else
                                    thisData(subEn,:) = reshape( omega32Times(:,thisSubInds,:), [size(omega32Times,1)*size(thisSubInds,2)*size(omega32Times,3)], 1 );
                                    %Unlike above, keep NaNs, for squareness sake (Note: Assumes that omega32Times is square/rectangular)
                                end
                            end
                            thisMedianMedian = nanmedian(thisData,2);
                            %nansum(thisData,2) %Testatory, to check correctness of data acquisition (A pseudo-hash of sort)
                            %thisSEM = nanstd(thisData,[], 2) / sqrt(size(thisData,2)); %Size-based n count; Weak to pooling
                            thisSEM = nanstd(thisData,[], 2) ./ sqrt( nansum(~isnan(thisData),2) ); %Non-NaN based n count
                
                            %Stats
                            %signofthetimes
                            %[P,ANOVATAB,STATS]=anova1([statsData],[statsGroups], 'off');
                            [P,ANOVATAB,STATS]=anova1(thisData', [], 'off');
                            statCompare = multcompare(STATS, 'CType', 'bonferroni', 'display', 'off');
                
                            if plotAllowance >= 2
                            figure
                            errorbar(thisMedianMedian,thisSEM,'.','Color',folColours(fol,:))
                            xlim([0,size(thisInds,1)+1])
                            xticks([1:size(thisInds,1)])
                            %xticklabels(nStruct(enny).name) %Note: Currently only as guaranteed as correct indices specified
                            xticklabels(thisNames)
                            %title([nStruct(enny).name,' comparison'])
                
                            %to add: stats, correct indexing checks, etc
                
                            %Stats display, Mk 2
                            yLim = get(gca,'YLim');
                            ylim([yLim(1),yLim(2)+0.25*range(yLim)]) %Add 25% on top
                            statPos = linspace(yLim(2),yLim(2)+0.2*range(yLim),size(statCompare,1)); %Defines possible places for putting sig bar
                            for statI = 1:size(statCompare,1)
                                thisAlt = statPos(statI);
                                if statCompare(statI,6) < 0.05 %Hardcoded alpha
                                    line([statCompare(statI,1),statCompare(statI,2)],[thisAlt,thisAlt],'Color','k', 'LineWidth',1.5)
                                else
                                    line([statCompare(statI,1),statCompare(statI,2)],[thisAlt,thisAlt],'Color',[0.5,0.5,0.5],'LineStyle',':')
                                end
                            end
                
                            titleStr = [];
                            for i = 1:size(thisNames,2)%size(nStruct(enny).name,2)
                                %titleStr = [titleStr,nStruct(enny).name{i},' vs '];
                                titleStr = [titleStr,thisNames{i},' vs '];
                            end
                            %title([nStruct(enny).name,' comparison (p=',num2str(P),')'])
                            titleStr = [titleStr,' comparison (p=',num2str(P),')'];
                            if poolMode
                                titleStr = [titleStr,' [Pooled]'];
                            end
                            if logi
                                titleStr = [titleStr,' (Log)'];
                            end
                            if naisu == 1
                                titleStr = [titleStr,char(10),'Merged isom.'];
                            elseif naisu == 2
                                titleStr = [titleStr,char(10),'Isomer 1'];
                            elseif naisu == 3
                                titleStr = [titleStr,char(10),'Isomer 2'];
                            end
                            title(titleStr)
                            set(gcf,'Name',nameStr);
                            end

                            %Pair histograms
                            %paradigm
                            if plotAllowance >= 2 %Starting to get a lot of plots at this tier...
                                
                                figure
                                for subEn = 1:size(thisData,1)
                                    subplot( ceil( sqrt( size(thisData,1) ) ) , ceil( sqrt( size(thisData,1) ) ), subEn )
                                    histogram( thisData(subEn,:) , ceil(128/enny) )
                                    xlabel('Time')
                                    title(thisNames{subEn})
                                    ylabel('Count')
                                end
                                set(gcf,'Name',nameStr);

                            end


        
                            if naisu == 1
                            nStruct(enny).data = thisData;
                            nStruct(enny).SEM = thisSEM; %Note: Relates to data
                            nStruct(enny).medianMedian = thisMedianMedian;
                            nStruct(enny).stats = statCompare;
                            elseif naisu == 2 || naisu == 3
                            nIsoStruct.nBack(enny).nStruct(naisu-1).data = thisData;
                            nIsoStruct.nBack(enny).nStruct(naisu-1).SEM = thisSEM; %Note: Relates to data
                            nIsoStruct.nBack(enny).nStruct(naisu-1).medianMedian = thisMedianMedian;
                            nIsoStruct.nBack(enny).nStruct(naisu-1).stats = statCompare; 
                            %{
                            elseif naisu == 3
                            nIsoStruct.nBack(enny).nStruct(2).data = thisData;
                            nIsoStruct.nBack(enny).nStruct(2).SEM = thisSEM; %Note: Relates to data
                            nIsoStruct.nBack(enny).nStruct(2).medianMedian = thisMedianMedian;
                            nIsoStruct.nBack(enny).nStruct(2).stats = statCompare;
                            %}
                            end
                            
                
                        end
                    end
                end
    
                overStruct(fol).matt(logi+1).nStruct = nStruct;
                overStruct(fol).matt(logi+1).nIsoStruct = nIsoStruct;

                %thomas

                %Super-testatory comparison with Scott SEs
                %{
                figure; create_seq_eff_plot(RAllFlies.medianFixTimes.',[],'errors',RAllFlies.semFixTimes.');title('test Median');
                hold on
                errorbar( overStruct.matt(1).nStruct(5).medianMedian, overStruct.matt(1).nStruct(5).SEM, 'LineWidth', 2 ,'Marker','X')
                legend({'Scott SEs','Matt SEs'})
                %}
            
            
            end
        
            %clash
        
        end
    
    %fol end
    end
    
    %% Cross-folder comparisons
    
    if size(dataPath,2) > 1
    
        clearvars('-except', flagParamSaveList{:}) %NO BLEED
    
        %Prepare some things
        exLabels = [];
        for fol = 1:size(overStruct,2)
            exLabels{fol} = overStruct(fol).datasetName;
        end
    
        %All RT (Scott) comparison
            %Note: Functionally 'pooled' nature
        if isfield(overStruct(1).scott,'RT_Hist_Array') %Check first element for existence
            collAllRT = [];
            temp = [];
            for fol = 1:size(overStruct,2)
                collAllRT{fol} = overStruct(fol).scott.RT_Hist_Array;
                temp = [temp,size(overStruct(fol).scott.RT_Hist_Array,1)];
            end
            collAllRTFlat = nan(nanmax(temp),size(overStruct,2));
            for fol = 1:size(overStruct,2)
                collAllRTFlat(1:temp(fol),fol) = collAllRT{fol};
            end
    
            if plotAllowance >= 2
            %Orig
            figure
            for fol = 1:size(collAllRTFlat,2)
                subplot(size(collAllRTFlat,2),1,fol)
                histogram(collAllRTFlat(:,fol),'NumBins', 25, 'FaceColor',folColours(fol,:)) %No guarantee same binning
                title([overStruct(fol).datasetName,' All RT'])
                xlabel('RT (ms)')
                ylabel('Count')
            end
    
            %Log
            if ~fixTimesMedianSubtraction
            figure
            for fol = 1:size(collAllRTFlat,2)
                subplot(size(collAllRTFlat,2),1,fol)
                thisData = log10( collAllRTFlat(:,fol) );
                thisMean = nanmean(thisData);
                thisMedian = nanmedian(thisData);
                histogram(thisData,'NumBins', 25,'FaceColor',folColours(fol,:)) %No guarantee same binning
    
                yLim = get(gca,'YLim');
                xLim = get(gca,'XLim');
                line([thisMean,thisMean],[0,nanmax(yLim)],'linewidth',2,'color','k','LineStyle','-'); 
                text(nanmax(xLim)*0.75,nanmax(yLim)*0.9,['mean=',num2str(thisMean,3)],'FontSize',14);
                line([thisMedian,thisMedian],[0,nanmax(yLim)],'linewidth',2,'color','r','LineStyle','--'); 
                text(nanmax(xLim)*0.75,nanmax(yLim)*0.85,['median=',num2str(thisMedian,3)],'FontSize',14);            
                title([overStruct(fol).datasetName,' All RT Log'])
    
                xlabel('log10RT')
                ylabel('Count')
            end
            end
            end
    
        end
    
        if plotAllowance >= 2
        %Matt [super]Median calcs
        for logi = 1:2
            if logi == 2 && fixTimesMedianSubtraction
                continue %Do not do log if median subtracted
            end
        
            if logi == 2
                logiLabel = ' [Log]';
            else
                logiLabel = ' [Original]';
            end

            plotData = [];
            plotN = []; %Only used (currently) by pooled mode
            for fol = 1:size(overStruct,2)
                if ~poolMode
                    plotData{fol} = squeeze( nanmedian( overStruct(fol).matt(logi).omega32Times, [1,2] ) ); 
                else
                    plotData{fol} = reshape( overStruct(fol).matt(logi).omega32Times , ...
                        [size(overStruct(fol).matt(logi).omega32Times,1)*size(overStruct(fol).matt(logi).omega32Times,2)*size(overStruct(fol).matt(logi).omega32Times,3)], 1 ); 
                    plotData{fol}( isnan(plotData{fol}) ) = [];
                    plotN(fol) = nansum(~isnan(plotData{fol}));
                end
                plotMedian(fol) = nanmedian(plotData{fol});
                plotSEM(fol) = nanstd(plotData{fol}) / sqrt(size(plotData{fol},1)); %Note: Unlike overSEM, this is not based on underlying total n
            
                %S
                overStruct(1).plots.medianRTs.overall(logi).datasetName{fol} = overStruct(fol).datasetName;
                overStruct(1).plots.medianRTs.overall(logi).data{fol} = plotData{fol}; %Reminder that terminal architecture is for separation by fol (not logi)
                
            end
    
            statsData = [];
            statsGroups = [];
            for fol = 1:size(overStruct,2)
                temp = plotData{fol};
                temp( isnan(temp) ) = [];
                statsData = [statsData; temp];
                statsGroups = [statsGroups; repmat(fol,numel(temp),1)]; %Note: Will fail hard if data not n rows x 1 col
            end
            [P,ANOVATAB,STATS]=anova1([statsData], [statsGroups], 'off');
            statCompare = multcompare(STATS, 'CType', 'bonferroni', 'display', 'off');
    
            figure
            %S
            overStruct(1).plots.medianRTs.overall(logi).logiState = logiLabel;
            overStruct(1).plots.medianRTs.overall(logi).pooledState = poolMode;

            %errorbar(plotMedian,plotSEM,'+')
            hold on
            iqDeets = [];
            for fol = 1:size(plotMedian,2)
                if ~doBoxplots
                errorbar(fol, plotMedian(fol),plotSEM(fol),'+','Color',folColours(fol,:))
                else
                %boxplot( plotData{fol}, 'Positions', fol,'Color',folColours(fol,:) )
                boxplot( plotData{fol}, 'Positions', fol,'Color',folColours(fol,:) , 'Whisker', 0, 'Symbol', '') %No whiskers, no outliers
                iqDeets = [iqDeets; nanmedian(plotData{fol}) - iqr(plotData{fol})*1.2 , nanmedian(plotData{fol}) + iqr(plotData{fol})*1.2 ];
                end
            end
            hold off
            xlim([0,size(overStruct,2)+1])
            xticks([1:size(overStruct,2)])
            xticklabels(exLabels)
            yLabelStr = ['Median RT (ms)'];
            if logi == 2 %Different from above, for indexing simplicity
                yLabelStr = [yLabelStr,' [Log]'];
            end
            ylabel(yLabelStr)
            titleStr = ['Comparative median RT (p=',num2str(P),')'];
            if ~poolMode
                titleStr = [titleStr, ' N= ' num2str([overStruct.N])];
            else
                titleStr = [titleStr, ' N= ' num2str([plotN]), ' [Pooled]'];
            end
            if logi == 2
                titleStr = [titleStr,' [Log]'];
            end
            title(titleStr)
            %S
            overStruct(1).plots.medianRTs.overall(logi).median = plotMedian;
            overStruct(1).plots.medianRTs.overall(logi).SEM = plotSEM;

            if doBoxplots
                ylim([nanmin(iqDeets(:,1)),nanmax(iqDeets(:,2))]) %Use more bespoke yLim for boxplots
            end
    
            yLim = get(gca,'YLim');
            if ~doBoxplots
            ylim([yLim(1),yLim(2)+0.25*range(yLim)]) %Add 25% on top
            statPos = linspace(yLim(2),yLim(2)+0.2*range(yLim),size(statCompare,1)); %Defines possible places for putting sig bar
            else
            ylim([yLim(1),yLim(2)+0.1*range(yLim)])
            statPos = linspace(yLim(2),yLim(2)+0.05*range(yLim),size(statCompare,1)); %Defines possible places for putting sig bar
            end
            for statI = 1:size(statCompare,1)
                thisAlt = statPos(statI);
                if statCompare(statI,6) < 0.05 %Hardcoded alpha
                    line([statCompare(statI,1),statCompare(statI,2)],[thisAlt,thisAlt],'Color','k', 'LineWidth',1.5)
                else
                    line([statCompare(statI,1),statCompare(statI,2)],[thisAlt,thisAlt],'Color',[0.5,0.5,0.5],'LineStyle',':')
                end
            end
            %rhcp

            %S
            overStruct(1).plots.medianRTs.overall(logi).stats = statCompare;

            %laraControlVsSD_medianRTs_overall = overStruct(1).plots.medianRTs.overall;

    
        end
        end
    
    
        %Matt nBack calcs
        nStruct = overStruct(1).matt(1).nStruct;
        %rmfield(nStruct,'data') %Attempt to remove 'spurious' data
    
        %Borrow above architecture, with mods
        xMod = linspace(-0.125,0.125,size(overStruct,2));
        for enny = 1:size(nStruct,2) 
            if ~isempty(nStruct(enny).inds)
                thisInds = nStruct(enny).inds;

                %S
                overStruct(1).plots.nBackOverall(enny).seqNames = nStruct(enny).name;

                thisData = [];
                thisMedianMedian = [];
                thisSEM = [];
                for fol = 1:size(overStruct,2)
                    %thisData = [];
                    for subEn = 1:size(thisInds,1)
                        thisSubInds = thisInds(subEn,:);
                        if ~poolMode
                            thisData{fol}(subEn,:) = squeeze( nanmedian( overStruct(fol).matt(1).omega32Times(:,thisSubInds,:) , [1,2] ) )'; %Hardcoded orig
                                %Row = Seq sub/superset, Col = Fly
                        else
                            thisData{fol}(subEn,:) = reshape( overStruct(fol).matt(1).omega32Times(:,thisSubInds,:), ...
                                1, [size(overStruct(fol).matt(1).omega32Times,1)*size(thisSubInds,2)*size(overStruct(fol).matt(1).omega32Times,3)] );
                                %Note: As above, rectangularity is probably guaranteed by rectangularity of omega32Times
                        end
                    end
                    %QA
                    if any(thisData{fol} == 0, 'all')
                        ['## Alert: Zeroes detected in pooled data; Potential rectangularity failure ##'];
                        crash = yes
                    end
                    thisMedianMedian{fol} = nanmedian(thisData{fol},2);
                    %thisSEM{fol} = nanstd(thisData{fol},[], 2) / sqrt(size(thisData{fol},2)); %Size-based
                    thisSEM{fol} = nanstd(thisData{fol},[], 2) ./ sqrt( nansum(~isnan(thisData{fol}),2) ); %Non-NaN based
                
                    %S
                    overStruct(1).plots.nBackOverall(enny).datasetName{fol} = overStruct(fol).datasetName;

                end
    
                %Stats
                    %Note: Performed more bespokely than above (For within datasets)
                %{{
                PVals = [];
                statCompare = {};
                bootlegPostHocBonf = 1; %Whether to multiply P vals by num comparisons to maintain a sort of Bonferroni nature
                %haydn
                for subEn = 1:size(thisInds,1)
                    statsData = [];
                    statsGroups = [];
                    for fol = 1:size(overStruct,2)
                        temp = thisData{fol}(subEn,:)';
                        temp( isnan(temp) ) = [];
                        statsData = [statsData; temp]; %Unlike above, this groups data of one seq across dataset 
                        statsGroups = [statsGroups; repmat(fol,numel(temp),1)]; %Note: Will fail hard if data not n rows x 1 col
                    end
                    [preVal,ANOVATAB,STATS]=anova1([statsData], [statsGroups], 'off');
                    preStat = multcompare(STATS, 'CType', 'bonferroni', 'display', 'off');
                    if bootlegPostHocBonf
                        preVal = preVal * size(thisInds,1);
                        preStat(:,6) = preStat(:,6) * size(thisInds,1); %Multiply P value by number of tests
                    end
                    PVals(subEn) = preVal;
                    statCompare{subEn} = preStat;
                end
                %}
    
    
                if plotAllowance >= 2
                %Plot
                figure
                hold on
                iqDeets = [];
                for fol = 1:size(overStruct,2)
                    if ~doBoxplots
                    errorbar( [1:size(thisMedianMedian{fol},1)]+xMod(fol), thisMedianMedian{fol}, thisSEM{fol}, '+', 'Color', folColours(fol,:) )
                    else
                    boxplot( thisData{fol}' , 'Positions', [1:size(thisMedianMedian{fol},1)]+xMod(fol), 'Color',folColours(fol,:) , 'Whisker', 0, 'Symbol', '')
                    iqDeets = [iqDeets; nanmin(nanmedian(thisData{fol}',1)) - nanmax(iqr(thisData{fol}'))*1.2 , nanmax(nanmedian(thisData{fol}',1)) + nanmax(iqr(thisData{fol}'))*1.2 ];
                    end
                end
                xlim([0,size(thisInds,1)+1])
                xticks([1:size(thisInds,1)])
                xticklabels(nStruct(enny).name)
                if ~doBoxplots %Doesn't seem to work with boxplots
                legend(exLabels)
                end

                if doBoxplots
                    ylim([nanmin(iqDeets(:,1)),nanmax(iqDeets(:,2))]) %Use more bespoke yLim for boxplots
                end
    
                titleStr = [];
                for i = 1:size(nStruct(enny).name,2)
                    titleStr = [titleStr,nStruct(enny).name{i},' vs '];
                end
                %title([nStruct(enny).name,' comparison (p=',num2str(P),')'])
                %P = NaN;
                titleStr = [titleStr,' comparison (p=',num2str(PVals),')']; %Note multiple P vals
                if poolMode
                    titleStr = [titleStr,char(10),' [Pooled]'];
                end
                if bootlegPostHocBonf
                    titleStr = [titleStr,char(10),' [Bootleg post-hoc Bonff. corr.]'];
                end
                %if logi
                %    titleStr = [titleStr,' (Log)'];
                %end
                title(titleStr)
                %set(gcf,'Name',nameStr);
                legend('AutoUpdate','off')
                legend('boxoff') 
            
                yLim = get(gca,'YLim');
                for subEn = 1:size(statCompare,2)
                    if ~doBoxplots
                    ylim([yLim(1),yLim(2)+0.25*range(yLim)]) %Add 25% on top
                    statPos = linspace(yLim(2),yLim(2)+0.2*range(yLim),size(statCompare,1)); %Defines possible places for putting sig bar
                    else
                    ylim([yLim(1),yLim(2)+0.1*range(yLim)])
                    statPos = linspace(yLim(2),yLim(2)+0.05*range(yLim),size(statCompare,1)); %Defines possible places for putting sig bar
                    end
                    for statI = 1:size(statCompare{subEn},1)
                        thisAlt = statPos(statI);
                        if statCompare{subEn}(statI,6) < 0.05 %Hardcoded alpha
                            
                            %line([statCompare{subEn}(statI,1),statCompare{subEn}(statI,2)],[thisAlt,thisAlt],'Color','k', 'LineWidth',1.5)
                            line([subEn + xMod(statCompare{subEn}(statI,1)),subEn + xMod(statCompare{subEn}(statI,2))],...
                                [thisAlt,thisAlt],'Color','k', 'LineWidth',1.5)
                        else
                            %line([statCompare{subEn}(statI,1),statCompare{subEn}(statI,2)],[thisAlt,thisAlt],'Color',[0.5,0.5,0.5],'LineStyle',':')
                            line([subEn + xMod(statCompare{subEn}(statI,1)),subEn + xMod(statCompare{subEn}(statI,2))],[thisAlt,thisAlt],'Color',[0.5,0.5,0.5],'LineStyle',':')
                        end
                    end
                end


                %Isomer plot of same
                %televisor
                for naisu = 1:2
                    %Stats as done above
                    PVals = [];
                    statCompare = {};
                    %bootlegPostHocBonf = 1; %Use above inline specification
                    for subEn = 1:size(thisInds,1)
                        statsData = [];
                        statsGroups = [];
                        for fol = 1:size(overStruct,2)
                            temp = overStruct(fol).matt(1).nIsoStruct.nBack(enny).nStruct(naisu).data(subEn,:)';
                            temp( isnan(temp) ) = [];
                            statsData = [statsData; temp]; %Unlike above, this groups data of one seq across dataset 
                            statsGroups = [statsGroups; repmat(fol,numel(temp),1)]; %Note: Will fail hard if data not n rows x 1 col
                        end
                        [preVal,ANOVATAB,STATS]=anova1([statsData], [statsGroups], 'off');
                        preStat = multcompare(STATS, 'CType', 'bonferroni', 'display', 'off');
                        if bootlegPostHocBonf
                            preVal = preVal * size(thisInds,1);
                            preStat(:,6) = preStat(:,6) * size(thisInds,1); %Multiply P value by number of tests
                        end
                        PVals(subEn) = preVal;
                        statCompare{subEn} = preStat;
                    end

                    %Plot
                    figure
                    hold on
                    thisData = [];
                    for fol = 1:size(overStruct,2)
                        if ~doBoxplots
                        thisIsoMedian = overStruct(fol).matt(1).nIsoStruct.nBack(enny).nStruct(naisu).medianMedian;
                        errorbar( [1:size(thisIsoMedian,1)]+xMod(fol), thisIsoMedian, thisSEM{fol}, '+', 'Color', folColours(fol,:) )
                        else
                        thisData{fol} = overStruct(fol).matt(1).nIsoStruct.nBack(enny).nStruct(naisu).data;
                        boxplot( thisData{fol}' , 'Positions', [1:size(thisMedianMedian{fol},1)]+xMod(fol), 'Color',folColours(fol,:) , 'Whisker', 0, 'Symbol', '')
                        iqDeets = [iqDeets; nanmin(nanmedian(thisData{fol}',1)) - nanmax(iqr(thisData{fol}'))*1.2 , nanmax(nanmedian(thisData{fol}',1)) + nanmax(iqr(thisData{fol}'))*1.2 ];
                        end
                    end
                    xlim([0,size(thisInds,1)+1])
                    ylabel('Reaction time (ms)')
                    xticks([1:size(thisInds,1)])
                    xticklabels(nStruct(enny).name)
                    if ~doBoxplots
                    legend(exLabels)
                    legend('AutoUpdate','off')
                    legend('boxoff') 
                    end
                    %title(['Isomer ',num2str(naisu)])
                    titleStr = [];
                    %for i = 1:size(nStruct(enny).name,2)
                    %    titleStr = [titleStr,nStruct(enny).name{i},' vs '];
                    %end
                    if enny <= 3
                    titleStr = ['Isomer ',num2str(naisu),' nStruct ',num2str(enny),' comparison (p=',num2str(PVals),')']; %Note multiple P vals
                    else
                    titleStr = ['Isomer ',num2str(naisu),' nStruct ',num2str(enny),' comparison (p=<multicular>)'];
                    end
                    if poolMode
                        titleStr = [titleStr,char(10),' [Pooled]'];
                    end
                    if bootlegPostHocBonf
                        titleStr = [titleStr,char(10),' [Bootleg post-hoc Bonff. corr.]'];
                    end
                    title(titleStr)                    
                
                    yLim = get(gca,'YLim');
                    for subEn = 1:size(statCompare,2)
                        if ~doBoxplots
                        ylim([yLim(1),yLim(2)+0.25*range(yLim)]) %Add 25% on top
                        statPos = linspace(yLim(2),yLim(2)+0.2*range(yLim),size(statCompare,1)); %Defines possible places for putting sig bar
                        else
                        ylim([yLim(1),yLim(2)+0.1*range(yLim)])
                        statPos = linspace(yLim(2),yLim(2)+0.05*range(yLim),size(statCompare,1)); %Defines possible places for putting sig bar
                        end
                        for statI = 1:size(statCompare{subEn},1)
                            thisAlt = statPos(statI);
                            if statCompare{subEn}(statI,6) < 0.05 %Hardcoded alpha
                                line([subEn + xMod(statCompare{subEn}(statI,1)),subEn + xMod(statCompare{subEn}(statI,2))],...
                                    [thisAlt,thisAlt],'Color','k', 'LineWidth',1.5)
                            else
                                line([subEn + xMod(statCompare{subEn}(statI,1)),subEn + xMod(statCompare{subEn}(statI,2))],[thisAlt,thisAlt],'Color',[0.5,0.5,0.5],'LineStyle',':')
                            end
                        end
                    end


                end %naisu end

                end %plotAllowance end

                %pylot

                %S
                overStruct(1).plots.nBackOverall(enny).medianData = thisMedianMedian;
                overStruct(1).plots.nBackOverall(enny).semData = thisSEM;

                
            end
        end
    %cross-folder end
    end

%window loop end
end

%Cross-window analyses
    %Note: Currently as written cannot be done for N=1 dataset

if size(winStruct(1).windowRTActual,2) > 1 && isfield(winStruct(1),'thisMedian') 
    disp(['Commencing cross-window analyses'])

    for logi = 1:2
        if logi == 2 && fixTimesMedianSubtraction
            continue %Do not do log if median subtracted
        end
        if logi == 2
            logiLabel = ' [Log]';
        else
            logiLabel = ' [Original]';
        end
        %S
        overStruct(1).plots.iterative(logi).state = logiLabel;
        subWin = winStruct(logi);

        %Median of median medians (Median)
        figure
        hold on
        for fol = 1:size(subWin(1).thisMedian,2) %Use first as guide
            %temp = reshape( [subWin.thisMedian], size(winStruct(1).logi(logi).thisMedian,2), size(winStruct,2) ); %Necessary because MATLAB bad at structure -> matrices
            %thisIsMedian = temp(fol,:);
            temp = [subWin.thisMedian];
            thisIsMedian = temp(:,fol);
            %temp = reshape( [winStruct.thisSEM], size(winStruct(1).thisSEM,2), size(winStruct,2) ); %Necessary because MATLAB bad at structure -> matrices
            %thisIsSEM = temp(fol,:);
            temp = [subWin.thisSEM];
            thisIsSEM = temp(:,fol);
    
            errorbar([subWin.windowRTActual],[thisIsMedian],[thisIsSEM])
            %Save-related stuff
            overStruct(1).plots.iterative(logi).median(fol).data = thisIsMedian; %As implied by indexing, elements under median field are individual folders (not logi state)
            overStruct(1).plots.iterative(logi).SEM(fol).data = thisIsSEM;
        end
        xlim([nanmin([subWin.windowRTActual])-1, nanmax([subWin.windowRTActual])+1])
        xticks([subWin.windowRTActual])
        xlabel(['windowRT (degs)'])
        ylabel(['Median RT',logiLabel]);
        legend(subWin.datasetName)
        titleStr = ['Iterative window plot'];
        if poolMode
            titleStr = [titleStr,' [Pooled]'];
            %Sav
            overStruct(1).plots.iterative(logi).pooled(1).state = ' [Pooled]'; %Note funny indexing
        end
        if logi == 2
            titleStr = [titleStr,' (Log)'];
            %Sav
            %overStruct(1).plots.iterative(logi).logi(1).state = ' (Log)';
        end
        title(titleStr)
        %Save iterative details
        %overStruct(1).plots.iterative.edges = histi.BinEdges;
        %overStruct(1).plots.iterative.values = histi.Values;
        temp = getframe;
        overStruct(1).plots.iterative(logi).frame = temp.cdata;
        overStruct(1).plots.iterative(logi).files(1).data = overStruct(1).files;
        overStruct(1).plots.iterative(logi).files(2).data = overStruct(2).files;
        overStruct(1).plots.iterative(logi).params = paramStruct;
            %Note: Since iterative plot operates comparatively, it is only saved to index 1 of overStruct, for simplicity

    end


end

if doAutomatedSavePlots
    figPath = baseFigPath; %Maybe in future make each dataset save to own separate folder?
    automatedSavePlots(figPath,'clearOldFigures',1,'automatedSaveVectors',saveVectors,'closeFiguresAfterSaving',closeFigures)
    if alsoSaveOverStruct
        save([figPath,filesep,'overStruct.mat'],'overStruct')
        disp(['(overStruct also saved)'])
    end
end

disp([char(10),'Fin'])

