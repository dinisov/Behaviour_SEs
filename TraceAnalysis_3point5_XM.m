% Supplementary/Additional analysis of Scott/Lara data, by Matt
%   "Tracer"
%Portions of this code based on BehaviourSEs2_MattBranch/Fix1
%Mk 1 - Core functionality
%   2 - Scott exclusion, general improvements
%   3 - Multi-dataset capabilities
%    .5 - Plot tweaks etc, Arbitrary trace shortening

%  ML_data_SD* datemodified:05/04/2024 .. 01/06/2024
%ML_Fix_SD* datemodified:05/04/2024 .. 01/06/2024
%       Windows explorer code to find Lara data

clear
close all
addpath('../Functions');
addpath('C:\Users\uqmvan13\ANALYSIS\Matt\Scripts'); %Supplies automatedSavePlots

%% Parameters

utiliseScottExclusion = 1;

%Stock Scott parameters
% "Yurei"
limitAngle = 90;             %degrees
windowRT = 45;              %degrees
refixQualTime = 300;        %loops
endTrace = 792;             %loops; If Scott Exclusion is in effect, this will be used to 'cap' the max length of a transient, as per his code's operation
limitReQual = 180 %windowRT +0;   %degrees
timeOut = 4000;             %ms

%% Matt exclusion parameters

utiliseMattExclusion = 0; %Whether to implement (Matt) exclusion
    %Note that currently, exclusion criteria will always be calculated
posRotationThreshold = 720; %Degs, how far to allow a bar position to ever be (Intended to counter infinispins)
wingbeatSDThreshold = 2; %How many SD to tolerate the majority of a wingbeat value being

%% Matt other parameters

arbitraryTraceShortening = 12; %If non-zero, arbitrarily remove this many frames from the start of every 'transient'
    %This is in reaction to 1 in 5 of Scott's data showing ~0.6s delays at the OL->CL period
    %Thus, it seems safest to unilaterally apply a very short (~20ms) slicing off of the start
        %Technically this means that the transients are now not perfectly identical between Scott/Matt analysis, but that is probably forgivable
    %Secondary note: For processing reasons, this will also be applied to endTrace to shorten it (If applicable)

%interpolateMeanTransients = 1; %Whether to interpolate mean transients to a more accurate time reference
    %Note: Abandoned

useNTrialExclusion = 0; %Whether to exclude datasets with less than a specified number of nTrials
if useNTrialExclusion
    nTrialThresh = 120; %Empirical; Based on Scott
end

calculateFirstEffectualSaccade = 0;
if calculateFirstEffectualSaccade
    effectualSaccadeSizeLimit = 6000; %Poly intersection calcs become v slow after a point, thus use this to limit which transients are analysed
end

%doIndividualPlots = 0; %Might make granular later?
plotAllowance = 4; %4 - All, 3 - Minor, 2 - Minor, 1 - Just cross-dataset, 0 - ???
    %Note: Unlike the previous doIndividualPlots system, (most) math is done anyway, so as to facilitate overStruct-related things
    %Secondary note: Currently only 4 does anything (Everything else is full allowed)
rawDataException = 0; %Whether to show raw data plots anway

%Plot stuff
numBins = 16; %Defines the number of bins for histograms
maxTime = 16; %Similar, defines the X limit for RT plots

doAutomatedSavePlots = 0; %Whether to use automatedSavePlots to save all figs
if doAutomatedSavePlots
    baseFigPath = 'C:\Users\uqmvan13\ANALYSIS\Scott\FIGOUT';
    alsoSaveOverStruct = 1; %Whether to save overStruct alongside
end

%%

%Lara Control vs SD
%{
dataPaths = ...
    {'C:\Users\uqmvan13\ANALYSIS\Scott\Data\Lara\Control',...
     'C:\Users\uqmvan13\ANALYSIS\Scott\Data\Lara\SD'};
%}
%Scott vs Lara control
%{
dataPaths = ...
    {'C:\Users\uqmvan13\ANALYSIS\Scott\Data\Scott\Data2',...
    'C:\Users\uqmvan13\ANALYSIS\Scott\Data\Scott\Data2'};
%}
%Singular
dataPaths = {'C:\Users\uqmvan13\ANALYSIS\Scott\Data\Scott\Data2'};

%%
%Some quick interim calcs
endTraceActual = endTrace - arbitraryTraceShortening; %If zero, no effect, otherwise shorten the max length of a trace by the arbitrary amount imposed above

%%

%Prevent potential data bleed
flagParamSaveList = who;
%flagParamSaveStruct = struct;
%for i = 1:size(flagParamSaveList,1)
%    eval(['flagParamSaveStruct.',flagParamSaveList{i},' = ', flagParamSaveList{i},';'])
%end
flagParamSaveList = [flagParamSaveList; 'flagParamSaveList'; 'fol'; 'overStruct'];

overStruct = struct;
%Dataset loop
for fol = 1:size(dataPaths,2)
    clearvars('-except', flagParamSaveList{:})

    
    %% Additional params
    
    %dataPath = 'C:\Users\uqmvan13\ScottCode\Matlab\Data24'; %Scott testatory set?
    %dataPath = 'C:\Users\uqmvan13\ScottCode\Data2'; %Scott thesis data?
    %dataPath = "C:\Users\uqmvan13\ScottCode\LaraRecapitulation\Split\SD";
    
    dataPath = dataPaths{fol};

    %% Pre-report
    temp = split(dataPath,filesep);
    datasetName = temp{end};
    disp([char(10),'--------------------------------',char(10),'Processing dataset ',datasetName])

    overStruct(fol).datasetName = datasetName;
    
    
    %% Import
    
    disp(['Searching ',char(dataPath),' for valid data'])

    overStruct(fol).ancillary.dataPath = dataPath;
    
    files = dir(dataPath);
    files = files(~[files.isdir],:);
    datas = cell(length(files),1);
    
    disp(['Detected files:'])
    for i = 1:length(files)
        %datas{i} = table2array(readtable(fullfile('D:\group_vanswinderen\Lara\Control\Data',files(i).name)));
        datas{i} = table2array(readtable(fullfile(dataPath,files(i).name)));
        disp(files(i).name);
    end
    
    nFlies = length(datas);
    disp([char(10),'Total # of detected files: ',num2str(nFlies)]);
    overStruct(fol).ancillary.nFlies = nFlies;
    if nFlies == 0
        ['## Alert: No detected flies ##']
        crash = yes %Probably warranted
    end
    
    %% Collate
    
    FLIES = struct;
    
    for fly = 1:nFlies
        disp([char(10),'Now analysing ',files(fly).name])
    
        %% Preprocess
    
        %Name/Identity
        disp(['Attempting to derive unique fly identity'])
        thisName = files(fly).name;
        thisName = strrep( thisName, ' ', ''); %Remove whitespace
        thisNameSep = strsplit(thisName,'_');
        prosDate = [];
        prosFly = [];
        prosExp = [];
        for i = 1:size(thisNameSep,2)
            temp = thisNameSep{i};
            %disp([num2str(i),' - ',temp])
            %Date
            %str2num(temp)
            if size(temp,2) == 6 && ~isempty(str2num(temp)) %DDMMYY presumably && Can be converted to number
                try
                    datetime(temp,'InputFormat','ddMMyy','TimeZone','+10'); %Attempt to convert to datetime, mostly as a check
                    prosDate = temp;
                catch
                    ['Prospective date cannot be converted to datetime']
                end
            end
            %Fly
            if contains(temp,'fly') || contains(temp,'f')
                if i ~= size(thisNameSep,2)
                    %try
                    if ~isempty( str2num(thisNameSep{i+1}) )
                        prosFly = thisNameSep{i+1};
                    else
                        ['fly keyword found but next string element not apparent #']
                    end
                    %catch
                    %end
                end
            end
            %Experiment
            if contains(temp,'exp')
                if i ~= size(thisNameSep,2)
                    temp2 = thisNameSep{i+1};
                    temp2 = strrep(temp2,'.csv','');
                    if ~isempty( str2num(temp2) )
                        prosExp = temp2;
                    else
                        ['Experiment keyword found but next string element not apparent #']
                        crosh
                    end
                end
            end
        end
    
        thisDate = [];
        thisFly = [];
        thisExp = [];
    
        thisDate = prosDate;
        try
            thisFly = str2num(prosFly);
        end
        try
            thisExp = str2num(prosExp);
        end
    
        FLIES(fly).date = thisDate;
        FLIES(fly).fly = thisFly;
        FLIES(fly).exp = thisExp;
    
        thisDatasetName = strcat('File #',num2str(fly),' - Date ',thisDate,'-f',num2str(thisFly),'-e',num2str(thisExp));

        overStruct(fol).individuals(fly).thisDatasetName = thisDatasetName;
    
        %% Process
    
    
        %Data
        data = datas{fly}; %All columns of data file
    
        timePoints = data(:,1); %Posix with milliseconds, presumably
        barPos = data(:,4); %Pixels[?], panel reference frame[?]
        wingBeat = data(:,5); %AU[?]; Only moderately confident this is indeed wingbeat amplitude/etc
        state = data(:,7); %Believe to be the scott-defined 'state' of stimulus
    
        timePointsRel = [timePoints-timePoints(1)]/1000; %Convert to relative seconds (from start)
    
    
        %Testatory figure
        %if doIndividualPlots || rawDataException
        if plotAllowance >= 4 || rawDataException
            figure
            plot(timePointsRel,barPos) %Di
            xlabel(['Time (s)'])
            ylabel(['Bar pos. (px)'])
            xlim([0,nanmax(timePointsRel)])
            legend('Bar pos')
            %xlim([0,0.25*nanmax(timePointsRel)])
            hold on
            plot(timePointsRel,wingBeat*20 - 200)
            %legend([{'Bar pos'},{'Wingbeat (norm.)'}])
            %legend('AutoUpdate','off')
            %hold on
            plot(timePointsRel,state*10 - 150,'Color','k')
            legend([{'Bar pos'},{'Wingbeat (norm.)'},{'State'}])
            legend('AutoUpdate','off')
            line([0,timePointsRel(end)],[15,15],'LineStyle','--','Color','r')
            line([0,timePointsRel(end)],[-15,-15],'LineStyle','--','Color','r')
            %line([0,timePointsRel(end)],[30,30],'LineStyle','--','Color','g')
            %line([0,timePointsRel(end)],[-30,-30],'LineStyle','--','Color','g')
            line([0,timePointsRel(end)],[45,45],'LineStyle','--','Color','b')
            line([0,timePointsRel(end)],[-45,-45],'LineStyle','--','Color','b')
            line([0,timePointsRel(end)],[60,60],'LineStyle','--','Color','m')
            line([0,timePointsRel(end)],[-60,-60],'LineStyle','--','Color','m')
            line([0,timePointsRel(end)],[0,0],'LineStyle',':','Color','k','LineWidth',1.5)
            temp = barPos; temp( temp ~= -60 ) = NaN;
            plot(timePointsRel,temp,'Color','c','LineWidth',3)
            temp = barPos; temp( temp ~= 60 ) = NaN;
            plot(timePointsRel,temp,'Color','m','LineWidth',3)
            title([thisDatasetName,' : Raw data'])
            set(gcf,'Name',[thisDatasetName,' : Raw data'])
        end
    
        %Processing
        %'Transients'
        stateTemp = state;
        stateTemp(stateTemp == 3) = 1; %Force state 3 ("Post-OL, not yet fixating") to be boolean true
        stateTemp(stateTemp == 2) = 0; %Disregard OL portions
        stateBW = bwlabel(stateTemp);

        nTrials = nanmax(stateBW); %More or less

        disp(['There were ',num2str(nTrials),' apparent OL->Fixation periods'])

        if useNTrialExclusion
            if nTrials < nTrialThresh
                disp(['-# Dataset lacks sufficient nTrials (',num2str(nTrials),' < ',num2str(nTrialThresh),') #-'])
                continue
            end
        end

        FLIES(fly).nTrials = nTrials; %Needs to be done after potential exclusion
    
        transStruct = struct;
        transColours = jet(nanmax(stateBW)); %Make a unique colour for each transient
        for i = 1:nanmax(stateBW)
            transStruct(i).ID = i; %Useful for post-exclusion checks possibly
    
            theseInds = find( stateBW == i );
            if utiliseScottExclusion
                if length(theseInds) > endTrace
                    theseInds = theseInds(1:endTrace); %'Cap' out max length
                    %NOTE: IT IS INTENTIONAL TO USE endTrace HERE, OTHERWISE THE ARBITRARY SHORTENING WOULD BE APPLIED FUNCTIONALLY TWICE
                end
            end
            %Derive 'actual' theseInds, for use with bar pos/time/wingbeat
                %It is unwise to mess with theseInds, since it is critical to sequence derivation/etc
            if arbitraryTraceShortening ~= 0
                theseIndsActual = theseInds( arbitraryTraceShortening+1:end ); %"Eliminate first <arbitraryTraceShortening> elements from start of indices"
            else
                theseIndsActual = theseInds;
            end
            %QA whether to not process this transient because too short
            if isempty(theseIndsActual)
                disp(['-# Alert: Transient #',num2str(i),' empty indices #-'])
                if arbitraryTraceShortening ~= 0
                    disp(['(Original indices length: ',num2str(numel(theseInds)),'; Post-shortening: ',num2str(numel(theseIndsActual)),')'])
                end
                continue %Note: Will need to account for empty transients if doing this
            end
    
            %thisTrans = barPos( theseInds ); %Original
            %thisTime = timePointsRel( theseInds );
            %thisWing = wingBeat( theseInds );
            thisTrans = barPos( theseIndsActual );
            thisTime = timePointsRel( theseIndsActual );
            thisWing = wingBeat( theseIndsActual );
    
            %blirg = unwrap(deg2rad(collTrans),[],2);
            thisTrans = unwrap(deg2rad(thisTrans),[],1); %Unwrap
            thisTrans = rad2deg(thisTrans); %Return to deg
            thisRelTime = thisTime - thisTime(1);
    
            transStruct(i).trans = thisTrans;
            transStruct(i).time = thisTime;
            transStruct(i).wing = thisWing;
    
            transStruct(i).size = size(thisTrans,1);
            transStruct(i).relTime = thisRelTime;
            transStruct(i).mean = nanmean(thisTrans);
            transStruct(i).wingMean = nanmean(thisWing);
    
            transStruct(i).colour = transColours(i,:);
    
    
            %Sequence derivation
                %Note: As written by Scott, an OL period is a cycling between Pos=200 (Null) -> 60 (Left?) -> 200 -> -60 (Right?) (Actual seq randomised ofc)
                %Secondary note: The CL period appears to actually start concommitant with the final stimulus, as in the last stim is presented and immediately control is given back
            seq = [];
            seqIdent = [];
            lastStimSide = [];
            %Derive inter-stimulus period
            if theseInds(1) == 1 %No sequence for first period (Based on N=1 data)
                transStruct(i).seq = [];
                transStruct(i).seqIdent = [];
                transStruct(i).lastStimSide = [];
                continue %Hopefully should skip rest of loop for first trial
            end
            isiEnd = theseInds(1)-1; %Note: Assumption
            isiStart = find( state(1:theseInds(1)-1) == 1, 1, 'last' )+1;
            %QA
            if isempty(isiStart) || isempty(isiEnd) || state(isiEnd) ~= 2 || numel( unique(barPos(isiStart:isiEnd)) ) > 3
                ['## Alert: Fatal error in pre-trial (#',num2str(i),') stimulus collection ##']
                ['isiStart: ',num2str(isiStart),', isiEnd: ',num2str(isiEnd),', pos. uniques: ',num2str(unique(barPos(isiStart:isiEnd)'))]
                crash = yes
            end
            temp = barPos( isiStart : isiEnd );
            temp2 = temp;
            temp2(temp2 == 200) = 0;
            temp2(temp2 == 60) = 128;
            temp2(temp2 == -60) = 128;
            temp2 = bwlabel(temp2);
            seq = [];
            for k = 1:nanmax(temp2)
                seq = [seq,temp(find(temp2 == k,1,'first'))]; %Find position at each stim cycle onset                  
            end
            seq = [seq,round(barPos(theseInds(1)),0)]; %Append final stim
            %QA
            if numel(seq) ~= 5 %Assumes nBack of 5
                ['## Alert: Critical failure to identify all (apparent) stim. elements ##']
                crash = yes
            end
            %if doIndividualPlots
            if plotAllowance >= 4
                disp(['Trial #',num2str(i),' sequence (pos): ',num2str(seq)])
            end
            seq( seq == -60 ) = 0;
            seq( seq == 60 ) = 1;
            %if doIndividualPlots
            if plotAllowance >= 4
                disp(['Sequence (binary): ',num2str(seq)])
            end
            seqIdent = bin2dec(num2str(seq)) + 1; %Functionally equivalent to the variable named "seq" in behaviourSEs2
    
            transStruct(i).seq = seq;
            transStruct(i).seqIdent = seqIdent;
            transStruct(i).lastStimSide = seq(end);
            transStruct(i).lastStimPos = round(barPos(theseInds(1)),0); %Borrow above math
    
            %Other ancillary metrics
            %story
            phaseCrossIndex = [-360,0,360]; %Not as dynamic as it seems, given hardcoded booleans below wrt > <
            phaseCrosses = [0,0,0]; %Crosses of -360, 0, and 360 respectively
            phaseCrossTime = [NaN,NaN,NaN]; %Time of above crosses
            %-360
            if any(thisTrans < -360)
                phaseCrosses(1) = 1;
                phaseCrossTime(1) = transStruct(i).relTime( find(thisTrans < -360, 1, 'first') );
            end
            %0
            if ( transStruct(i).lastStimPos == -60 && any(thisTrans > 0) ) %Longer bool to account for direction that fixating fly needs to move towards zero
                phaseCrosses(2) = 1;
                phaseCrossTime(2) = transStruct(i).relTime( find(thisTrans > 0, 1, 'first') );
            elseif ( transStruct(i).lastStimPos == 60 && any(thisTrans < 0) )
                phaseCrosses(2) = 1;
                phaseCrossTime(2) = transStruct(i).relTime( find(thisTrans < 0, 1, 'first') );
            end
            %360
            if any(thisTrans > 360)
                phaseCrosses(3) = 1;
                phaseCrossTime(3) = transStruct(i).relTime( find(thisTrans > 360, 1, 'first') );
            end
    
            transStruct(i).phaseCrosses = phaseCrosses;
            transStruct(i).phaseCrossTime = phaseCrossTime;
    
            %Sum (simple)
            thisSum = sum(abs(thisTrans)); %This uses the abs the deg value from 0 as a metric
                %It does not account for a fixation on a phase+1 or phase-1 
                %It is however agnostic of the starting position (-60 or 60), as long as possible starting positions are equidistant of 0
            transStruct(i).sumSimple = thisSum;
            %Sum (phase-associated)
                %Note: Could also easily calculate this always for all phases (even if not crossed)
            thisPhaseSum = nan(1,numel(phaseCrossIndex));
            if phaseCrosses(1) == 1
                thisPhaseSum(1) = sum( abs( thisTrans - phaseCrossIndex(1) ) ); %Calculate distance from -360
            end
            if phaseCrosses(2) == 1
                thisPhaseSum(2) = sum( abs( thisTrans ) ); %Ditto, 0
            end
            if phaseCrosses(3) == 1
                thisPhaseSum(3) = sum( abs( thisTrans - phaseCrossIndex(3) ) ); %Ditto, 0
            end
            transStruct(i).sumPhase = thisPhaseSum;
    
            %Half-'rise' to first inflection
            if calculateFirstEffectualSaccade
                firstEffectualSaccadeEnd = [];
                firstEffectualSaccadeHalftime = [];
                firstEffectualSaccadeTotalSpeed = [];
                %waves
                if transStruct(i).size < effectualSaccadeSizeLimit 
                    temp = diff(thisTrans);
                    tempx = [1:size(temp,1)];
                    temp2 = zeros(size(temp,1),1);
                    [exxy, ~] = polyxpoly([tempx], [temp], [tempx], [temp2] ); %Use polynomial solving to find places where derivative intersects ~0
                    exxyRound = round(exxy,0); %Minor risk of very close values causing issues here (e.g. 221.2 and 221.45 both round to 221 etc)
                    
                    firstZeroReturnIntent = find( abs(thisTrans) < abs(transStruct(i).lastStimPos), 1, 'first' );
                        %Note: Currently this only signals a move to zero, not to full phase (in the inefficient direction)
                    if ~isempty(firstZeroReturnIntent) && ~isempty(exxy) %"If moved towards zero, and non-zero number of 'saccades' (defined by stabilisation periods)
                        firstInd = find( exxy > firstZeroReturnIntent, 1, 'first' );
                        firstEffectualSaccadeEnd = exxyRound(firstInd); %Find ind of first point after intent where stabilised
                        firstEffectualSaccadeHalftime = transStruct(i).relTime( round( firstEffectualSaccadeEnd/2 ,0) ); %Rounding a rounded number but eh...
                        if firstInd ~= 1
                            firstEffectualSaccadeTotalSpeed = ...
                                abs( ( thisTrans( exxyRound(firstInd) ) - thisTrans( exxyRound(firstInd-1) ) ) / (thisRelTime( exxyRound(firstInd) ) - thisRelTime( exxyRound(firstInd-1) )) );
                                %Find total degs between end of first effectual saccade and previous stability point, then divide by time taken to give deg/s for first effectual saccadde
                        else
                            firstEffectualSaccadeTotalSpeed = ...
                                abs( ( thisTrans( exxyRound(firstInd) ) - transStruct(i).lastStimPos) / thisRelTime( exxyRound(firstInd) ) ); 
                                %As above, but simply for time from start to current degrees value of saccade end
                        end
                    end
                else
                    disp(['-# Transient too long (',num2str(transStruct(i).size),'); First effectual saccade pos. not calculated #-'])
                end
                
                %Testatory figure
                %{
                figure
                temp3 = normalize(thisTrans)*0.5;
        
                %{
                %Incomplete multi-axis attempt
                t=tiledlayout(1,1);
                ax1=axes(t);
                set(ax1,'XLim',[min(tempx) max(tempx)]);
                %plot(ax1,x1,y)
                %ax1.XTick=x1;
                %ax1.YTick=[];
                %ax1.XTickLabel=x1;
                %xlabel('Frame')
                %}
        
                plot(temp3)
                hold on
                plot(temp)
                plot(temp2,'LineStyle',':','Color','k')
                scatter(exxy,temp(round(exxy,0)));
                xlabel('Frame')
                if ~isempty(firstEffectualSaccadeEnd)
                    scatter(firstEffectualSaccadeEnd,temp(firstEffectualSaccadeEnd),'filled')
                    line([firstEffectualSaccadeEnd/2,firstEffectualSaccadeEnd/2],[0,1],'Color','g','LineStyle',':','LineWidth',1.5)
                    legend({['This transient'],['Trans deriv.'],['Zero diff line'],['Intersection/s'],['First saccade end'],['First saccade end halftime']})
                else
                    legend({['This transient'],['Trans deriv.'],['Zero diff line'],['Intersection/s']})
                end
        
                %{
                ax2=axes(t);
                set(ax2,'XLim',[min(transStruct(i).relTime) max(transStruct(i).relTime)]);
                %plot(ax2,x1,y)
                ax2.YTick=[];
                ax2.XAxisLocation='top';
                %ax2.XTick=x1;
                %ax2.XTickLabel=x2;
                xlabel('Time (s)')
                %}
        
                %legend({['This transient'],['Trans deriv.'],['Zero diff line'],['Intersection/s']})
                %}
        
                transStruct(i).firstSaccadeEnd = firstEffectualSaccadeEnd; %i.e. Time to first effective movement towards zero
                transStruct(i).firstSaccadeHalfTime = firstEffectualSaccadeHalftime; %i.e. Halved time from zero to this movement
                transStruct(i).firstSaccadeTotalSpeed = firstEffectualSaccadeTotalSpeed; %i.e. Total speed (deg/s) *during* this saccade (Probs superior to halftime)
            end
    
                    
            clear thisTrans thisWing thisTime
        end
    
        %----------

        %Quick pre-removal of transients that hard-failed for one reason or another
        contExcCount = 0;
        for i = size(transStruct,2):-1:1
            if isempty(transStruct(i).trans)
                transStruct(i) = [];
                contExcCount = contExcCount + 1;
            end
        end
        if contExcCount ~= 0
            disp([num2str(contExcCount),' transient/s excluded due to probably continue reasons (i.e. Too short post-shortening)'])
        end
    
        %Exclude if applicable
        excMattCount = 0;
        excMattInds = [];
        excScottCount = 0;
        excScottInds = [];
        for i = 1:size(transStruct,2)
            exClusion = [0,0,0]; %Abs. pos exceeding threshold, Wingbeat SD threshold, Hybrid Scott Too-Short
            exe = 0;
            exScotClusion = [0,0,0,0]; %>limitAngle, leaves windowRT during refix, no refix, refix >timeout
            exeS = 0;
    
            thisTrans = transStruct(i).trans;
            thisWing = transStruct(i).wing;
            thisTime = transStruct(i).time;
    
            %Matt
            %Absolute position
            if any(abs(thisTrans) > posRotationThreshold) %"Any bar position > specified threshold"
                exClusion(1) = 1;
                exe = 1;
            end
    
            %Wingbeat SD
            if nansum( abs([ thisWing + nanmedian([transStruct.wingMean]') ]) > nanstd([transStruct.wingMean]')*wingbeatSDThreshold ) > 0.5*size(thisWing,1) %"More than half wingbeats exceed SD threshold"
                    %Note: Uses median of individual wingbeat means, to cut down on length bias of potentially aberrant transients
                    %Note: High risk of sign issues causing this to not work as intended
                exClusion(2) = 1;
                exe = 1;
            end
    
            %Hybrid Scott Too-Short detection
                %Not valid to do here because not yet calculated loopRT
            %if size(thisTrans,1) <  ( loopRT + refixQualTime )
            %    exClusion(3) = 1;
            %    exe = 1;
            %end
    
            %Moved below to support Scott hybrid exclusion
            %transStruct(i).exClusion = exClusion;
            %excMattCount = excMattCount + exe;
            %if exe
            %    excMattInds = [excMattInds,i];
            %end
    
            %Scott
            %Preparatory (directly copied from behaviourSEs2)
                %NOTE: MAY NEED TO IMPLEMENT ANTI-FIRST TRIAL HERE, TO REMOVE SPURIOUS INSTA-FIXATION ETC
            rtTrace = [];
            if i > 1
                window = [-windowRT windowRT];
                %betweenBounds = thisTrace > window(1) & thisTrace <  window(2);
                betweenBounds = thisTrans > window(1) & thisTrans <  window(2); %Note: Not 100% guaranteed that thisTrans == thisTrace
                %outsideMaxAngle = thisTrace > limitAngle & thisTrace <  -limitAngle;
                outsideMaxAngle = thisTrans > limitAngle & thisTrans <  -limitAngle;
                if any(betweenBounds)
                    loopRT = find(betweenBounds,1); %Loop index of first betweenbounds
                    %neptr
                    %fixTimes(trial, seq) = thisFly.time(thisFly.closedLoopStart(trial)+loopRT-1)- thisFly.time(thisFly.closedLoopStart(trial)); %Time RT in ms
                    thisFixTime = thisTime(loopRT) - thisTime(1); %Functionally same as Scott calcs above...hopefully
                    %rtTrace = thisTrace(1:loopRT);
                    rtTrace = thisTrans(1:loopRT);
                    %if(loopRT + refixQualTime > endTrace)           % set up Refix Window time
                    if(loopRT + refixQualTime > endTraceActual)           % set up Refix Window time
                        %reQualTrace = thisTrace(loopRT:endTrace);
                        reQualTrace = thisTrans(loopRT:endTraceActual);
                    else
                        %reQualTrace = thisTrace(loopRT:loopRT + refixQualTime);
                        if loopRT + refixQualTime <= size(thisTrans,1) %Normal
                            reQualTrace = thisTrans(loopRT:loopRT + refixQualTime);
                        else %Aberrantly short?
                            reQualTrace = thisTrans(loopRT:end);
                            exClusion(3) = 1;
                            exe = 1;
                            disp(['Trans. ',num2str(i),' duration (',num2str(size(thisTrans,1)),') < loopRT + refixQualTime (',num2str(loopRT + refixQualTime),')'])
                        end
                    end
                        %Note: This appears to be a cutdown portion of the transient/trial at a 'fixed' duration from when refixation occurred
                            %'fixed' because if it is long (Bool result 1) then it is from refixation:endTrace(e.g.792), otherwise (Bool 2) it is refixation:refixQualTime(e.g.300)
                    if any(abs(rtTrace) > limitAngle)% if it wrapped around before entering window it's a bad trial
                        %badTrials(trial) = 1;
                        exScotClusion(1) = 1; %"Probably looped a phase before refixating" (More accurately, "Departed a nominal (e.g. 90 deg) distance from 0 before fixating" (e.g. Inefficient return) 
                    elseif any(abs(reQualTrace) > limitReQual)% if it did not stay in the requal window for limitQual it's a bad trial
                        %badTrials(trial) = 1;
                        exScotClusion(2) = 1; %"Left window during refixation"
                    end
                else
                    exScotClusion(3) = 1; %"Did not enter window" 
                end
        
                %if(fixTimes(trial, seq) > timeOut)
                if thisFixTime > timeOut/1000 %Note Scott time in ms
                    %badTrials(trial) = 1;
                    exScotClusion(4) = 1;
                end
        
                if any(exScotClusion)
                    exeS = 1;
                end
            end
    
            transStruct(i).exScotClusion = exScotClusion;
            transStruct(i).scottRTTrace = rtTrace; %Scott-specific refixation trace
                %Note: Some of Scott exclusion based on rtTrace, thus the 'remainder' of the trace, while odd looking may be 'valid' as defined by Scott
            excScottCount = excScottCount + exeS;
            if exeS
                excScottInds = [excScottInds,i];
            end
    
            transStruct(i).exClusion = exClusion;
            excMattCount = excMattCount + exe;
            if exe
                excMattInds = [excMattInds,i];
            end
    
        end
        %loathing
        %Report
        if excMattCount > 0
            disp([num2str(excMattCount),' transient/s failed one or more exclusion criteria'])
    
            %if doIndividualPlots
            if plotAllowance >= 4
                figure
                for i = 1:excMattCount
                    subplot( ceil(sqrt(excMattCount)), ceil(sqrt(excMattCount)), i )
                    hold on
                    thisI = excMattInds(i);
                    %Pos
                    yyaxis left
                    plot( transStruct(thisI).relTime, transStruct(thisI).trans )
                    line([0,transStruct(thisI).relTime(end)],[posRotationThreshold,posRotationThreshold],'Color','b','LineStyle',':')
                    line([0,transStruct(thisI).relTime(end)],[-posRotationThreshold,-posRotationThreshold],'Color','b','LineStyle',':')
                    ylabel('Bar pos. (deg)')
                    xlabel('Relative time (s)')
        
                    %Wingbeat SD
                    yyaxis right
                    plot( transStruct(thisI).relTime, transStruct(thisI).wing + nanmedian([transStruct.wingMean]') )
                    line([0,transStruct(thisI).relTime(end)],[nanstd([transStruct.wingMean]')*wingbeatSDThreshold,nanstd([transStruct.wingMean]')*wingbeatSDThreshold],'Color','r','LineStyle',':')
                    line([0,transStruct(thisI).relTime(end)],-[nanstd([transStruct.wingMean]')*wingbeatSDThreshold,nanstd([transStruct.wingMean]')*wingbeatSDThreshold],'Color','r','LineStyle',':')
                    ylabel('Adjusted wingbeat (a.u.)')
        
                    %Title
                    title(['Transient #',num2str(thisI),' - Exc. ',num2str(transStruct(thisI).exClusion)])
        
                end
                %falling
                %title([thisDatasetName,'   Exclusion'])
                set(gcf,'Name',[thisDatasetName,' Matt exclusion'])
            end
        end
        if excScottCount > 0
            disp([num2str(excScottCount),' transient/s failed one or more Scott exclusion criteria'])
    
            %if doIndividualPlots
            if plotAllowance >= 4
                figure
                for i = 1:excScottCount
                    subplot( ceil(sqrt(excScottCount)), ceil(sqrt(excScottCount)), i )
                    hold on
                    thisI = excScottInds(i);
                    %xieuha
                    %Pos
                    plot( transStruct(thisI).relTime, transStruct(thisI).trans )
                    %line([0,transStruct(thisI).relTime(end)],[posRotationThreshold,posRotationThreshold],'Color','b','LineStyle',':')
                    %line([0,transStruct(thisI).relTime(end)],[-posRotationThreshold,-posRotationThreshold],'Color','b','LineStyle',':')
                    line([0,transStruct(thisI).relTime(end)],[windowRT,windowRT],'Color','g','LineStyle','--')      
                    line([0,transStruct(thisI).relTime(end)],[-windowRT,-windowRT],'Color','g','LineStyle','--')
                    if transStruct(thisI).exScotClusion(1)
                        line([0,transStruct(thisI).relTime(end)],[limitAngle,limitAngle],'Color','r','LineStyle',':')      
                        line([0,transStruct(thisI).relTime(end)],[-limitAngle,-limitAngle],'Color','r','LineStyle',':')
                        title(['Trans. #',num2str(thisI),' - LimitAngle (',num2str(limitAngle),'degs)'])
                    elseif transStruct(thisI).exScotClusion(2) && limitReQual ~= Inf
                        line([0,transStruct(thisI).relTime(end)],[limitReQual,limitReQual],'Color','r','LineStyle',':')      
                        line([0,transStruct(thisI).relTime(end)],[-limitReQual,-limitReQual],'Color','r','LineStyle',':')
                        title(['Trans. #',num2str(thisI),' - LimitReQual (',num2str(limitReQual),'degs)'])
                    elseif transStruct(thisI).exScotClusion(3)
                        line([0,transStruct(thisI).relTime(end)],[windowRT,windowRT],'Color','r','LineStyle',':')      
                        line([0,transStruct(thisI).relTime(end)],[-windowRT,-windowRT],'Color','r','LineStyle',':')
                        title(['Trans. #',num2str(thisI),' - NeverWindow (',num2str(windowRT),'degs)'])            
                    elseif transStruct(thisI).exScotClusion(4)
                        yLim = get(gca,'YLim');
                        line([timeOut/1000,timeOut/1000],[yLim(1),yLim(2)],'Color','r','LineStyle',':')      
                        line([timeOut/1000,timeOut/1000],[yLim(1),yLim(2)],'Color','r','LineStyle',':')
                        title(['Trans. #',num2str(thisI),' - TimeOut (',num2str(timeOut/1000),'s)'])  
                    end
                    ylabel('Bar pos. (deg)')
                    xlabel('Relative time (s)')
                    %title([thisDatasetName,'   Exclusion'])
                    set(gcf,'Name',[thisDatasetName,' Scott exclusion'])
        
                    %Title
                    %title(['Transient #',num2str(thisI),' - Exc. ',num2str(transStruct(thisI).exClusion)])
        
                end
                %falling
            end
        end
        %And actually exclude, if requested
        if utiliseMattExclusion && ~isempty( excMattInds )
            %Old, risky system
            %{
            for i = size(transStruct,2):-1:1
                if ~isempty( find(excMattInds == i) )
                    if ~exist('excStruct') || isempty(excStruct)
                        excStruct = transStruct(i);
                    else
                        excStruct(i) = transStruct(i);
                    end
                    transStruct(i) = [];
                end
            end
            %}
            %New, index/ID-based system
            excRemov = 0;
            for i = 1:numel(excMattInds)
                temp = [transStruct.ID]; %Note: Critical to recapitulate this before each loop, since indices will change with exclusion
                thisExc = excMattInds(i);
                thisInd = find(temp == thisExc);
                if ~isempty(thisInd)
                    if ~exist('excStruct') || isempty(excStruct)
                        excStruct = transStruct(thisInd);
                    else
                        excStruct(thisInd) = transStruct(thisInd);
                    end
                    transStruct(thisInd) = [];     
                    excRemov = excRemov + 1;
                end
            end
            disp([num2str(excRemov),' transients excluded based on Matt criteria (',num2str(excRemov),' of ',num2str(excMattCount),' originally designated)'])
                %Note: Deviations from intended vs actual excluded may come from pre-exclusion, etc
        end
        if utiliseScottExclusion && ~isempty( excScottInds )
            %New, index/ID-based system
            excRemov = 0;
            for i = 1:numel(excScottInds)
                temp = [transStruct.ID]; %Note: Critical to recapitulate this before each loop, since indices will change with exclusion
                thisExc = excScottInds(i);
                thisInd = find(temp == thisExc);
                if ~isempty(thisInd)
                    if ~exist('excStruct') || isempty(excStruct)
                        excStruct = transStruct(thisInd);
                    else
                        excStruct(thisInd) = transStruct(thisInd);
                    end
                    transStruct(thisInd) = [];     
                    excRemov = excRemov + 1;
                end
            end
            disp([num2str(excRemov),' transients excluded based on Scott criteria (',num2str(excRemov),' of ',num2str(excScottCount),' originally designated)'])
                %Note: Deviations from intended vs actual excluded may come from pre-exclusion, etc
        end
        %Quick interim exclusion of first transient if no sequence data
        if ~isempty(transStruct) && isempty(transStruct(1).seq) %"Only do if first transient indeed contains no sequence"
            if  ~isempty( find(excMattInds == i) )
                excStruct = transStruct(1);
            else
                excStruct(1) = transStruct(1);
            end
            transStruct(1) = [];
            disp(['First transient excluded'])
        end
        %sun
    
        if utiliseMattExclusion | utiliseScottExclusion
            disp([num2str(size(transStruct,2)),' transients remain of ',num2str(nanmax(stateBW)),' original'])
        end
    
        %Abort if transStruct empty (all excluded, most likely)
        if isempty(transStruct)
            ['-# No valid transients detected #-']
            continue
        end
    
        FLIES(fly).transStruct = transStruct;
    
        %----------
    
        %Collate transient data
        maxSize = nanmax( [transStruct.size]' );
        medSize = nanmedian( [transStruct.size]' );
        %Quick clean
        if rem(medSize,1) ~= 0
            medSize = round(medSize,0); %Round to 0
            disp(['-# Caution: Median size had to be rounded to whole number (',num2str(medSize),') #-'])
        end
        if utiliseScottExclusion ~= 1 %Only applicable properly if not capping transients as per Scott
            disp(['Max transient size/duration: ',num2str(maxSize),'samples /',...
                num2str(transStruct( find([transStruct.size]' == maxSize) ).time(end) - transStruct( find([transStruct.size]' == maxSize) ).time(1)),'s'])
            disp(['Median transient size/duration: ',num2str(medSize),'samples /',...
                num2str(transStruct( find([transStruct.size]' == maxSize) ).time(medSize) - transStruct( find([transStruct.size]' == maxSize) ).time(1)),'s']) %Cheat and pull median timepoint from longest
        else
            %disp(['Max transient size/duration capped at ',num2str(endTrace),' samples (Scott exclusion)'])
            disp(['Max transient size/duration capped at ',num2str(endTraceActual),' samples (Scott exclusion)'])
            if arbitraryTraceShortening ~= 0
                disp(['(',num2str(endTrace),' before arbitrary shortening [',num2str(arbitraryTraceShortening),' frames])'])
            end
            disp(['(Min dur: ',num2str(nanmin([transStruct.size])),'; Med. dur: ',num2str(nanmedian([transStruct.size])),')'])
        end
        
        collTrans = nan(size(transStruct,2),maxSize); %Bar pos
        collWing = nan(size(transStruct,2),maxSize); %Wingbeat
        collTime = nan(size(transStruct,2),maxSize); %Rel time
        collLinear = []; %Basically ligated transients
        collWingLinear = []; %Basically ligated transients
        %collWingMeanSD = nan(size(transStruct,2),2);
        for i = 1:size(transStruct,2)
            thisTrans = transStruct(i).trans;
            thisWing = transStruct(i).wing;
            thisTime = transStruct(i).relTime;
    
            collTrans(i, 1:size(thisTrans,1) ) = thisTrans;
            collWing(i, 1:size(thisWing,1) ) = thisWing;
            collTime(i, 1:size(thisTime,1) ) = thisTime;
            collLinear = [collLinear;thisTrans]; %Note: Inefficient memory allocation
            collWingLinear = [collWingLinear;thisWing]; %Note: Inefficient memory allocation
            %collWingMeanSD(i,:) = [ nanmean(thisWing) , nanstd(thisWing) ];
        end

        overStruct(fol).individuals(fly).nTrans = size(collTrans,1);
        overStruct(fol).individuals(fly).transSize = size(collTrans,2);

        %And some other stuff
        phaseColl = reshape([transStruct.phaseCrosses], size(phaseCrosses,2), size(transStruct,2) )'; %Collated phase crossings
        phaseTimeColl = reshape([transStruct.phaseCrossTime], size(phaseCrosses,2), size(transStruct,2) )'; %Ditto, for time
        phaseSumColl = reshape([transStruct.sumPhase], size(phaseCrosses,2), size(transStruct,2) )'; %Ditto, for sum
        phaseSumNormColl = phaseSumColl ./ repmat( [transStruct.size]', 1, size(phaseSumColl,2) ); %Normalised sum by size of transient
        phaseSideColl = [transStruct.lastStimSide]'; %Mostly unnecessary to collate, but useful
        phaseSideDegsColl = [transStruct.lastStimPos]'; %Mostly unnecessary to collate, but useful

        %timeBase = nanmean(collTime,1);
        timeBase = nanmedian(collTime,1); %Use median, to reduce potential issues with mean where timing inconsistencies exist (See below QA)

        %Quick check of variability in times
        %if fly == 5
        %truecolours
        %end
        temp = collTime - timeBase;
        temp2 = nanstd(collTime,[],1);
        if any(any(abs(temp) > 8*temp2)) %"Are any elements of collTime more than 8SD deviant from mean collTime" (Note: 8SD is empirical, based on Scott data)
            ['-# Caution: Dataset #',num2str(fly),' One or more transients have timing outside of +- 8SD of all #-']
            overStruct(fol).individuals(fly).timeIssue = 1; %Note that this may not be critical, but it's worth knowing
        else
            overStruct(fol).individuals(fly).timeIssue = 0;
        end

        if plotAllowance >= 4
            %Testatory time plot
            figure
            hold on
            for i = 1:size(collTime,1)
                plot(collTime(i,:))
            end
            plot(timeBase,'LineWidth',2,'Color','k')
            title([datasetName,' - Dataset #',num2str(fly)])
        end

        %if fly == 5
        %alice
        %end

    
        %Sidedness
        unSides = unique([transStruct.lastStimSide]);
        unPos = []; %Parallel of above, but needs extra processing to make sure synchronised
        for unp = 1:numel(unSides)
            temp = [transStruct.lastStimPos];
            temp2 = temp( [transStruct.lastStimSide] == unSides(unp) );
            %Quick QA
            if numel(unique(temp2)) > 1
                ['## Alert: Potential contamination of lastStimPos detected ##']
                crash = yes %This implies an error in lastStimPos calcs and/or lastStimSide
            end
            unPos(unp) = unique(temp2);
        end
    
        %% Plots
        %if doIndividualPlots
    
            if plotAllowance >= 4
            %Occupancy histogram
            figure
            hist(collLinear,32)
            %xlim([-2000,2000])
            xlabel('Bar position')
            title('Pos. occupancy hist.')
            %title([thisDatasetName,'   Exclusion'])
            set(gcf,'Name',[thisDatasetName,' Occupancy hist'])
        
            %--------------------------
        
            %Plot all transients
            figure
            hold on
            for i = 1:size(collTrans,1)
                plot(collTime(i,:),collTrans(i,:),'Color',transStruct(i).colour)
                plot(collTime(i,:),collWing(i,:)*25 - 200,'Color',transStruct(i).colour,'LineStyle','--')
            end
            xlabel('Relative time (s)')
            ylabel('Bar pos. (deg)')
            title('All transients')
            %title([thisDatasetName,'   Exclusion'])
            set(gcf,'Name',[thisDatasetName,' All trans'])
            %xlim([0    9.9496])
            %ylim([ -546.1594  473.3915])
        
            %Plot histogram of transient durations
            figure
            hist([transStruct.size]',16)
            xlabel('Transient duration (samples)')
            %title([thisDatasetName,'   Exclusion'])
            set(gcf,'Name',[thisDatasetName,' Trans dur'])
        
            %--------------------------
        
            %Plot collated transients
                %Note: This plot significantly suffers from both sidedness and rolling decreases in N
            figure
            %timeBase = nanmean(collTime,1); %Moved above
            plot( timeBase , nanmean(collTrans,1), 'LineWidth', 2, 'Color', 'k' )
            hold on
            %errorbar( nanmean(collTime,1) , nanmean(collTrans,1), nanstd(collTrans,1) )
            %groupShadeCoordsX = [1:1:size(collTime,2),size(collTime,2):-1:1];
            groupShadeCoordsX = [timeBase(1:1:size(collTime,2)),timeBase(size(collTime,2):-1:1)];
            groupShadeCoordsY = [nanmean(collTrans,1)+(nanstd(collTrans,1)/sqrt(size(collTrans,1))),flip(nanmean(collTrans,1)-(nanstd(collTrans,1)/sqrt(size(collTrans,1))))];
            fill(groupShadeCoordsX, groupShadeCoordsY,'r') %Error shading
            alpha(0.15)
            xlim([0,nanmean(timeBase(medSize))]) %Artificially limit X axis to median duration
            xlabel(['Relative time (s)'])
            ylabel('Bar pos. (degs)')
            title(['Mean transient (Median X limit [Max true: ',num2str(timeBase(end)),'s])'])
            %title([thisDatasetName,'   Exclusion'])
            set(gcf,'Name',[thisDatasetName,' Mean trans'])
        
            %-----------------------------------------
            %plotAllowance
            end
        
            %Plot separated transients
            if plotAllowance >= 4
            figure
            end
            sideColours = jet(numel(unSides));
            for sidi = 1:numel(unSides)
                thisSide = unSides(sidi);
                if plotAllowance >= 4
                subplot(1,numel(unSides),sidi)
                end
                matchTrans = find( [transStruct.lastStimSide] == thisSide );
                %laine
                thisTime = timeBase;
                thisData = collTrans(matchTrans,:);
                thisMean = nanmean(thisData,1);
                thisSEM = nanstd(thisData,[],1) / sqrt(size(thisData,1));
                %thisSEM(isnan(thisSEM)) = 0;
                if any(isnan(thisMean))
                    %QA
                    if find( ~isnan(thisMean) == 1 , 1 , 'last' ) > find( isnan(thisMean) == 1 , 1 , 'first' )
                        ['## Alert: Potential mix of NaN and non-NaN elements ##']
                        %crash = yes %Not sure what case might cause this; Not critical
                    end
                    thisData = thisData(:, 1:find(isnan(thisMean),1,'first')-1); %Note slightly different indexing, because not collapsed
                    thisTime = thisTime(1:find(isnan(thisMean),1,'first')-1);
                    thisSEM = thisSEM(1:find(isnan(thisMean),1,'first')-1);
                    thisMean = thisMean(1:find(isnan(thisMean),1,'first')-1);
                end
                if plotAllowance >= 4
                plot( thisTime , thisMean, 'LineWidth', 2, 'Color', sideColours(sidi,:) )
                hold on
                groupShadeCoordsX = [thisTime(1:1:size(thisTime,2)),thisTime(size(thisTime,2):-1:1)];
                groupShadeCoordsY = [thisMean+thisSEM,flip(thisMean-thisSEM)];
                fill(groupShadeCoordsX, groupShadeCoordsY, sideColours(sidi,:)) %Error shading
                alpha(0.15)
                for i = 1:size(thisData,1)
                    plot(thisTime,thisData(i,:),'Color',transStruct(matchTrans(i)).colour)
                    alpha(0.15)
                end
                line([0,thisTime(end)],[0,0],'LineStyle','--','Color','k')
                if nanmax(thisData,[],'all') > 360
                    line([0,thisTime(end)],[360,360],'LineStyle','--','Color','k')            
                end
                if nanmin(thisData,[],'all') < -360
                    line([0,thisTime(end)],[-360,-360],'LineStyle','--','Color','k')  
                end
                xlabel('Relative time (s)')
                ylabel('Bar position (degs)')
                title(['Transient for stim ending ',num2str(thisSide)])
                end

                %Save this data for later use
                overStruct(fol).individuals(fly).sepTrans(sidi).thisTime = thisTime;
                overStruct(fol).individuals(fly).sepTrans(sidi).thisMean = thisMean;
                overStruct(fol).individuals(fly).sepTrans(sidi).thisSEM = thisSEM;
                    %Architecture: sepTrans - 2 sides
        
            end
            if plotAllowance >= 4
            %title([thisDatasetName,'   Exclusion'])
            set(gcf,'Name',[thisDatasetName,' Sep trans'])
            end
        
            %----------------------
        
            %Plot separated transients, additionally separated by ending position
            if plotAllowance >= 4
            figure
            end
            unSides = unique([transStruct.lastStimSide]);
            sideColours = jet(numel(unSides));
            for sidi = 1:numel(unSides)
                thisSide = unSides(sidi);
                %subplot(1,numel(unSides),sidi)
                for ph = 1:size(phaseCrosses,2) %Use definition from above
                    if plotAllowance >= 4
                    subplot(size(phaseCrosses,2),numel(unSides), sub2ind( [numel(unSides),size(phaseCrosses,2)], sidi, ph ) ) %Caved and used sub2ind
                    end
                    matchTrans = intersect( find( [transStruct.lastStimSide] == thisSide ), find( phaseColl(:,ph) == 1 ) );
                        %Return transients both ending on the correct side AND corresponding to this phase cross
                    %laine
                    thisMedianPhaseSums = NaN; %Prepare in case empty
                    thisSEMPhaseSums = NaN;
                    thisMedianPhaseTime = NaN;
                    thisSEMPhaseTime = NaN;
        
                    if ~isempty(matchTrans)
                        %Main data calcs
                        thisTime = timeBase;
                        thisData = collTrans(matchTrans,:);
                        thisMean = nanmean(thisData,1);
                        thisSEM = nanstd(thisData,[],1) / sqrt(size(thisData,1));
                        %thisSEM(isnan(thisSEM)) = 0;
                        if any(isnan(thisMean))
                            %No QA because done above
                            thisData = thisData(:, 1:find(isnan(thisMean),1,'first')-1); %Note slightly different indexing, because not collapsed
                            thisTime = thisTime(1:find(isnan(thisMean),1,'first')-1);
                            thisSEM = thisSEM(1:find(isnan(thisMean),1,'first')-1);
                            thisMean = thisMean(1:find(isnan(thisMean),1,'first')-1);
                        end
                        %Tack on phase cross sum stuff
                        temp = reshape([transStruct(matchTrans).sumPhase], size(phaseCrosses,2), size(matchTrans,1) )'; %Returns a collated list of phase sums for this matchTrans
                            %Note: It is possible that a match will have phase results for more than one crossing, so the indexing next is also necessary
                        thisMedianPhaseSums = nanmedian( temp(:,ph) );
                        thisSEMPhaseSums = nanstd( temp(:,ph) ) / sqrt( size(matchTrans,1) );
                        %Repeat for time
                        temp2 = reshape([transStruct(matchTrans).phaseCrossTime], size(phaseCrosses,2), size(matchTrans,1) )'; %Returns a collated list of phase sums for this matchTrans
                        thisMedianPhaseTime = nanmedian( temp2(:,ph) );
                        thisSEMPhaseTime = nanstd( temp2(:,ph) ) / sqrt( size(matchTrans,1) );
                        %ghost
        
                        if plotAllowance >= 4
                        plot( thisTime , thisMean, 'LineWidth', 2, 'Color', sideColours(sidi,:) )
                        hold on
                        groupShadeCoordsX = [thisTime(1:1:size(thisTime,2)),thisTime(size(thisTime,2):-1:1)];
                        groupShadeCoordsY = [thisMean+thisSEM,flip(thisMean-thisSEM)];
                        fill(groupShadeCoordsX, groupShadeCoordsY, sideColours(sidi,:)) %Error shading
                        alpha(0.15)
                        for i = 1:size(thisData,1)
                            plot(thisTime,thisData(i,:),'Color',transStruct(matchTrans(i)).colour)
                            alpha(0.15)
                            %Note: Could also use temp data above to report individual transient phase sums
                        end
                        end
                    end
                    if plotAllowance >= 4
                    line([0,thisTime(end)],[0,0],'LineStyle','--','Color',[0.25,0.25,0.25])
                    line([0,thisTime(end)*0.1],[unPos(sidi),unPos(sidi)],'LineStyle',':','Color','r','LineWidth',2)
                    line([0,thisTime(end)],[phaseCrossIndex(ph),phaseCrossIndex(ph)],'LineStyle',':','Color','g','LineWidth',2)
                    %if nanmax(thisData,[],'all') > 360
                    %    line([0,thisTime(end)],[360,360],'LineStyle','--','Color','k')            
                    %end
                    %if nanmin(thisData,[],'all') < -360
                    %    line([0,thisTime(end)],[-360,-360],'LineStyle','--','Color','k')  
                    %end
                    xlabel('Relative time (s)')
                    ylabel('Bar position (degs)')
                    %ylim([phaseCrossIndex(1)*1.25,phaseCrossIndex(3)*1.25,])
                    %%end
                    %title(['Stim ending ',num2str(thisSide),' - Phase cross: ',num2str(phaseCrossIndex(ph))])
                    title(['Stim ending ',num2str(thisSide),' (',num2str(unPos(sidi)),' degs) - Phase cross: ',num2str(phaseCrossIndex(ph)),...
                        ' [Time: ',num2str(thisMedianPhaseTime),'+- ',num2str(thisSEMPhaseTime),'s, Sum: ',num2str(thisMedianPhaseSums),'+- ',num2str(thisSEMPhaseSums),' degs]'])
                    end
                end
        
            end
            %title([thisDatasetName,'   Exclusion'])
            if plotAllowance >= 4
            set(gcf,'Name',[thisDatasetName,' Sep trans by side'])
            end
        
            %--------------------------
        
            %Barplot of phase reaction times
                %Note: Not separated by side?
            phaseTimeCollMedian = nanmedian(phaseTimeColl,1);
            phaseTimeCollSEM = nanstd(phaseTimeColl,[],1) ./ sqrt( nansum( ~isnan( phaseTimeColl ) , 1 ) );
            if plotAllowance >= 4
            figure
            errorbar([1:size(phaseTimeColl,2)],phaseTimeCollMedian,phaseTimeCollSEM, '.')
            xlim([0,size(phaseTimeColl,2)+1])
            xticks([1:size(phaseTimeColl,2)])
            xticklabels(phaseCrossIndex)
            xtickangle(90)
            xlabel('Phase crossed')
            ylabel('Reaction Time (s)')
            title(['Phase cross median RTs (n=',num2str(nansum( ~isnan( phaseTimeColl ) , 1 )),')'])
            grid on
            %legend(legNameList)
            %title([thisDatasetName,'   Exclusion'])
            set(gcf,'Name',[thisDatasetName,' Phase cross RTs'])
            end
        
            %--------------------------
        
            %Barplot of phase sum
                %Note: Not separated by side?
            phaseSumCollMedian = nanmedian(phaseSumColl,1);
            phaseSumCollSEM = nanstd(phaseSumColl,[],1) ./ sqrt( nansum( ~isnan( phaseSumColl ) , 1 ) );
            if plotAllowance >= 4
            figure
            errorbar([1:size(phaseSumColl,2)],phaseSumCollMedian,phaseSumCollSEM, '.')
            xlim([0,size(phaseSumColl,2)+1])
            xticks([1:size(phaseSumColl,2)])
            xticklabels(phaseCrossIndex)
            xtickangle(90)
            xlabel('Phase crossed')
            ylabel('Sum (delta degs)')
            title(['Phase cross median sums (n=',num2str(nansum( ~isnan( phaseSumColl ) , 1 )),')'])
            grid on
            set(gcf,'Name',[thisDatasetName,' Phase sums'])
            end
    
        %doIndividualPlots end
        %end
    
    
    
        %%
    
        %Save overdata
        FLIES(fly).phaseColl = phaseColl;
        FLIES(fly).phaseTimeColl = phaseTimeColl;
        FLIES(fly).phaseSumColl = phaseSumColl;
        FLIES(fly).phaseSumNormColl = phaseSumNormColl;
        FLIES(fly).phaseSideColl = phaseSideColl;
        FLIES(fly).phaseSideDegsColl = phaseSideDegsColl;

        FLIES(fly).collLinear = collLinear;
    
    end

    %casca
    %Check to see if any FLIES elements empty (Due to exclusion)
        %Borrowed from behaviourSEs2_MattBranch
    excCount = 0;
    for fly = size(FLIES,2):-1:1
        if isempty( FLIES(fly).nTrials )
            FLIES(fly) = [];
            excCount = excCount + 1;
        end
    end
    %nFlies = size(FLIES,2); %This is necessary to prevent issues with later analysis; Not actually used here
    if excCount ~= 0
        disp([num2str(excCount),' datasets excluded; ',num2str(size(FLIES,2)),' remain'])
    end

    %overStruct(fol).FLIES = FLIES; %Unnecessary to save (probably)
    
    
    %%
    
    %Cross-fly analyses
    %Note: Currently cannot handle all data of one side and/or plot type being NaN (e.g. N=1 fly analysed)
    
    %Prepare data
    superFlyColl = []; %Ancillary, to store fly identity
    superPhaseColl = [];
    superPhaseTimeColl = []; %Remin
    superPhaseSumColl = [];
    superPhaseSumNormColl = [];
    superPhaseSideColl = [];
    superPhaseSideDegsColl = [];
    superCollLinear = [];
    for fly = 1:size(FLIES,2)
        superPhaseColl = [superPhaseColl;FLIES(fly).phaseColl]; %Inefficient
        superFlyColl = [superFlyColl;repmat(fly,size(FLIES(fly).phaseColl,1),1)];
        superPhaseTimeColl = [superPhaseTimeColl;FLIES(fly).phaseTimeColl]; 
        superPhaseSumColl = [superPhaseSumColl;FLIES(fly).phaseSumColl];
        superPhaseSumNormColl = [superPhaseSumNormColl;FLIES(fly).phaseSumNormColl]; 
        superPhaseSideColl = [superPhaseSideColl;FLIES(fly).phaseSideColl]; 
        superPhaseSideDegsColl = [superPhaseSideDegsColl;FLIES(fly).phaseSideDegsColl]; 

        superCollLinear = [superCollLinear;FLIES(fly).collLinear]; 
    end
    
    %Some plot stuff
    %%numBins = 10; %Moved above
    %maxSum = nanmax(superPhaseSumColl,[],'all');
    medSum = nanmedian(superPhaseSumColl,'all');
    medSumNorm = nanmedian(superPhaseSumNormColl,'all');
    %%maxTime = 12; %Arbitrary, based on Scott paradigm settings vaguely
    
    %'Simple' plots
        %Now 50% more generalised!
    
    sourceData{1} = superPhaseTimeColl;
    sourceDataIndex{1} = 'Phase times'; %Not exactly correct here to specify zero-cross, since that is done in plotting
    sourceData{2} = superPhaseSumColl;
    sourceDataIndex{2} = 'Phase sums';
    sourceData{3} = superPhaseSumNormColl;
    sourceDataIndex{3} = 'Phase sums norm.';
    
    for plotI = 1:size(sourceData,2)
        
        %L vs R
        plotData = [];
        plotMean = [];
        plotMedian = [];
        %plotSEM = [];
        unSides = unique(superPhaseSideColl); %Recalculate in case of differential expression across individuals
        unColours = summer(numel(unSides));
    
        figure
        for unI = 1:numel(unSides)
            thisUn = unSides(unI);
            thisUnDeg = unique( superPhaseSideDegsColl( superPhaseSideColl == thisUn ) );
            %QA
            if numel(thisUnDeg) > 1
                ['## Alert: Differing side (degs) values obtained for same sides ##']
                crash = yes
                %Not quite sure what would cause this
            end
            %plotData{unI} = superPhaseTimeColl( superPhaseSideColl == thisUn, 2 ); %Collect phase=zero crossings for this side
            plotData{unI} = sourceData{plotI}( superPhaseSideColl == thisUn, 2 );
            %QA for no data
            if nansum(~isnan(plotData{unI})) == 0
                ['-# Caution: No non-NaN data detected #-']
                continue
            end
            plotMean(unI) = nanmean(plotData{unI}); %Not armoured against empty data
            plotMedian(unI) = nanmedian(plotData{unI});
            plotN(unI) = nansum( ~isnan(plotData{unI}) );
        
            subplot(1,numel(unSides),unI)
            if plotI == 1
                xLimVal = maxTime;
            elseif plotI == 2
                xLimVal = medSum*2;
            elseif plotI == 3
                xLimVal = medSumNorm*2;
            end
            theseBins = linspace(0,xLimVal,numBins); %Leave as hardcoded
            %histogram(plotData{unI},'NumBins', numBins,'BinEdges',theseBins);
            histogram(plotData{unI},'NumBins', numBins,'BinEdges',theseBins,'FaceColor',unColours(unI,:));
            yLim = get(gca,'YLim');
            xLim = get(gca,'XLim');
        
            %Tack on mean/median calcs
            line([plotMean(unI),plotMean(unI)],[0,nanmax(yLim)],'linewidth',2,'color','k','LineStyle','-'); 
            text(nanmax(xLim)*0.75,nanmax(yLim)*0.9,['mean=',num2str(plotMean(unI),3)],'FontSize',14);
            line([plotMedian(unI),plotMedian(unI)],[0,nanmax(yLim)],'linewidth',2,'color','r','LineStyle','--'); 
            text(nanmax(xLim)*0.75,nanmax(yLim)*0.85,['median=',num2str(plotMedian(unI),3)],'FontSize',14);
        
            %xlabel('Time to phase cross (s)')
            xlabel([sourceDataIndex{plotI}])
            ylabel('Count')
            %title(['Stim ending ',num2str(thisUn),' [',num2str(thisUnDeg),' degs] - n=',num2str(plotN(unI)),' zero phase cross times'])
            title([datasetName,' - Stim ending ',num2str(thisUn),' [',num2str(thisUnDeg),' degs] - n=',num2str(plotN(unI)),' zero phase ',sourceDataIndex{plotI}])
            set(gcf,'Name',[datasetName,' - Cross-fly Zero phase cross ',sourceDataIndex{plotI}])
        end
        %Follow on log plot reusing data
        %bodies
        %if nansum(~isnan(plotData{unI})) ~= 0
        minVal = nanmin( cell2mat(plotData') , [], 'all' ); %Min of both sides, used for below
        minLogVal = log10(minVal); %Log version
        if minLogVal == -Inf
            ['-# Alert: At least one 0 value existed in ',sourceDataIndex{plotI},' #-']
            %This will be ignored for log hist X limit purposes (Because log10(0) = Inf)
            minLogVal = 0; 
        end
        %end
        figure
        logData = [];
        logMean = [];
        logMedian = [];
        for unI = 1:numel(unSides)
            thisUn = unSides(unI);
            thisUnDeg = unique( superPhaseSideDegsColl( superPhaseSideColl == thisUn ) ); %Don't redo QA
        
            %QA for no data
            if nansum(~isnan(plotData{unI})) == 0
                ['-# Caution: No non-NaN data detected #-']
                continue
            end
            logData{unI} = log10(plotData{unI});
            logMean(unI) = nanmean(logData{unI});
            logMedian(unI) = nanmedian(logData{unI});
        
            subplot(1,numel(unSides),unI)
            if plotI == 1
                xLimVal = log10(maxTime);
            elseif plotI == 2
                xLimVal = log10(medSum*8);
            elseif plotI == 3
                xLimVal = log10(medSumNorm*8);
            end
            %theseBins = linspace(0,xLimVal,numBins);
            theseBins = linspace(minLogVal,xLimVal,numBins); %Will crash minVal exceeds xLimVal probs
            histogram(logData{unI},'NumBins', numBins,'BinEdges',theseBins,'FaceColor',unColours(unI,:));
            yLim = get(gca,'YLim');
            xLim = get(gca,'XLim');
        
            %Tack on mean/median calcs
            line([logMean(unI),logMean(unI)],[0,nanmax(yLim)],'linewidth',2,'color','k','LineStyle','-'); 
            text(nanmax(xLim)*0.75,nanmax(yLim)*0.9,['mean=',num2str(logMean(unI),3)],'FontSize',14);
            line([logMedian(unI),logMedian(unI)],[0,nanmax(yLim)],'linewidth',2,'color','r','LineStyle','--'); 
            text(nanmax(xLim)*0.75,nanmax(yLim)*0.85,['median=',num2str(logMedian(unI),3)],'FontSize',14);
        
            %xlabel('Time to phase cross (log s)')
            xlabel(['Log ',sourceDataIndex{plotI}])
            ylabel('Count')
            %title(['Stim ending ',num2str(thisUn),' [',num2str(thisUnDeg),' degs] - n=',num2str(plotN(unI)),' zero phase cross log times'])
            title([datasetName,' - Stim ending ',num2str(thisUn),' [',num2str(thisUnDeg),' degs] - n=',num2str(plotN(unI)),' zero phase log ',sourceDataIndex{plotI}])
            set(gcf,'Name',[datasetName,' - Cross-fly Zero phase cross log ',sourceDataIndex{plotI}])
        end
    
    
        %Repeat, but for All
        allData = cell2mat(plotData');
        allMean = nanmean(allData);
        allMedian = nanmedian(allData);
        allSEM = nanstd(allData) / sqrt(nansum(~isnan(allData))); %Use ~isnan for better count of non-NaN values
    
        %Original
        figure
        if plotI == 1
            xLimVal = maxTime;
        elseif plotI == 2
            xLimVal = medSum*2;
        elseif plotI == 3
            xLimVal = medSumNorm*2;
        end
        %theseBins = linspace(0,xLimVal,numBins);
        theseBins = linspace(minVal,xLimVal,numBins); %Switch to minVal rather than hardcoded 0
        histogram(allData,'NumBins', numBins,'BinEdges',theseBins,'FaceColor','k');
        hold on
        yLim = get(gca,'YLim');
        xLim = get(gca,'XLim');
    
        %Tack on mean/median calcs
        line([allMean,allMean],[0,nanmax(yLim)],'linewidth',2,'color','k','LineStyle','-'); 
        text(nanmax(xLim)*0.75,nanmax(yLim)*0.9,['mean=',num2str(allMean,3)],'FontSize',14);
        line([allMedian,allMedian],[0,nanmax(yLim)],'linewidth',2,'color','r','LineStyle','--'); 
        text(nanmax(xLim)*0.75,nanmax(yLim)*0.85,['median=',num2str(allMedian,3)],'FontSize',14);
    
        xlabel([sourceDataIndex{plotI}])
        ylabel('Count')
        title([datasetName,' - All stim - n=',num2str(nansum(~isnan(allData))),' zero phase ',sourceDataIndex{plotI}])
        set(gcf,'Name',[datasetName,' - Cross-fly Zero phase cross ',sourceDataIndex{plotI}])
    
    
        %Log
        allLogData = log10(allData);
        allLogData( allLogData == Inf ) = NaN; %Cleaning required to prevent mean from being contaminated by Inf
        allLogMean = nanmean(allLogData);
        allLogMedian = nanmedian(allLogData);
    
        figure
        if plotI == 1
            xLimLogVal = log10(maxTime);
        elseif plotI == 2
            xLimLogVal = log10(medSum*2);
        elseif plotI == 3
            xLimLogVal = log10(medSumNorm*2);
        end
        theseBins = linspace(minLogVal,xLimLogVal,numBins);
        histogram(allLogData,'NumBins', numBins,'BinEdges',theseBins,'FaceColor','k');
        hold on
        yLim = get(gca,'YLim');
        xLim = get(gca,'XLim');
    
        %Tack on mean/median calcs
        line([allLogMean,allLogMean],[0,nanmax(yLim)],'linewidth',2,'color','k','LineStyle','-'); 
        text(nanmax(xLim)*0.75,nanmax(yLim)*0.9,['mean=',num2str(allLogMean,3)],'FontSize',14);
        line([allLogMedian,allLogMedian],[0,nanmax(yLim)],'linewidth',2,'color','r','LineStyle','--'); 
        text(nanmax(xLim)*0.75,nanmax(yLim)*0.85,['median=',num2str(allLogMedian,3)],'FontSize',14);
    
        xlabel(['Log ',sourceDataIndex{plotI}])
            ylabel('Count')
        title([datasetName,' - All stim - n=',num2str(nansum(~isnan(allData))),' zero phase log ',sourceDataIndex{plotI}])
        %set(gcf,'Name',['Cross-fly Zero phase cross log ',sourceDataIndex{plotI}])
        set(gcf,'Name',[datasetName,' - Cross-fly Zero phase cross log ',sourceDataIndex{plotI}])

        %nagare

        thisPlotNameClean = matlab.lang.makeValidName(sourceDataIndex{plotI});
        thisPlotNameClean = strrep(thisPlotNameClean,'_','');
        overStruct(fol).plots.(thisPlotNameClean).plotData = plotData;
        %overStruct(fol).plots.(thisPlotNameClean).plotMean = plotMean; %Don't bother saving this stuff
        %overStruct(fol).plots.(thisPlotNameClean).plotMedian = plotMedian;
        overStruct(fol).plots.(thisPlotNameClean).logData = logData;

        overStruct(fol).plots.(thisPlotNameClean).allData = allData;
        %overStruct(fol).plots.(thisPlotNameClean).allMean = allMean;
        %overStruct(fol).plots.(thisPlotNameClean).allMedian = allMedian;
        overStruct(fol).plots.(thisPlotNameClean).allLogData = allLogData;
        %overStruct(fol).plots.(thisPlotNameClean).allLogMean = allLogMean;
        %overStruct(fol).plots.(thisPlotNameClean).allLogMedian = allLogMedian;
        overStruct(fol).plots.(thisPlotNameClean).xLimVal = xLimVal;
        overStruct(fol).plots.(thisPlotNameClean).xLimLogVal = xLimLogVal;
        overStruct(fol).plots.(thisPlotNameClean).minVal = minVal;
        overStruct(fol).plots.(thisPlotNameClean).minLogVal = minLogVal;


        %Overall position occupancy hist
        figure
        hist(superCollLinear,64)
        xlabel('Bar position')
        title([datasetName,' Pos. occupancy hist. (Med.=',num2str(nanmedian(superCollLinear),3),'; SD=',num2str(nanstd(superCollLinear),3),')'])
        set(gcf,'Name',[datasetName,' - Cross-fly Occupancy hist'])
        yLim = get(gca,'YLim');
        xLim = get(gca,'XLim');

        overStruct(fol).specialPlots.posOccupancyHist.allData = superCollLinear;
        overStruct(fol).specialPlots.posOccupancyHist.xLim = xLim;
        overStruct(fol).specialPlots.posOccupancyHist.yLim = yLim;
    
    
    end


%fol end
end


%% Cross-dataset comparisons
%------------------------------------------------------------------------------------------------------------------------------------

if size(overStruct,2) < 2
    disp(['-- Only one dataset specified for analysis; Ending --'])
    if doAutomatedSavePlots %Do this here if about to end
        figPath = baseFigPath; %Maybe in future make each dataset save to own separate folder?
        automatedSavePlots(figPath,'clearOldFigures',1)
        if alsoSaveOverStruct
            save([figPath,filesep,'overStruct.mat'],'overStruct')
            disp(['(overStruct also saved)'])
        end
    end
    return
end

clearvars('-except', flagParamSaveList{:})

disp([char(10),'-- Performing cross-dataset analysis --'])

disp([num2str(size(overStruct,2)),' datasets exist for comparison'])
for fol = 1:size(overStruct,2)
    disp([overStruct(fol).datasetName])
end

%disp(['Now comparing plots'])
folColours = jet(size(overStruct,2));

%%

plotTypes = fieldnames(overStruct(1).plots); %Note: Assumption of all plots done for all datasets
for plotI = 1:size(plotTypes,1)
    disp(['Now comparing ',plotTypes{plotI}])

    %Original
    %{
    figure
    for fol = 1:size(overStruct,2)
        thisPlotInd = sub2ind( [size(overStruct,2),2] , 1, fol ); %Why this seems to be in col, row order i do not know
        subplot(size(overStruct,2),2,thisPlotInd)

        allData = overStruct(fol).plots.(plotTypes{plotI}).allData;
        allMean = nanmean(allData);
        allMedian = nanmedian(allData);

        %onourside
        xLimVal = overStruct(fol).plots.(plotTypes{plotI}).xLimVal;
        theseBins = linspace(0,xLimVal,numBins);
        histogram(allData,'NumBins', numBins,'BinEdges',theseBins,'FaceColor','k');
        hold on
        yLim = get(gca,'YLim');
        xLim = get(gca,'XLim');
    
        %Tack on mean/median calcs
        line([allMean,allMean],[0,nanmax(yLim)],'linewidth',2,'color','k','LineStyle','-'); 
        text(nanmax(xLim)*0.75,nanmax(yLim)*0.9,['mean=',num2str(allMean,3)],'FontSize',14);
        line([allMedian,allMedian],[0,nanmax(yLim)],'linewidth',2,'color','r','LineStyle','--'); 
        text(nanmax(xLim)*0.75,nanmax(yLim)*0.85,['median=',num2str(allMedian,3)],'FontSize',14);
    
        xlabel([plotTypes{plotI}])
        ylabel('Count')
        title([overStruct(fol).datasetName,' - All stim - n=',num2str(nansum(~isnan(allData))),' zero phase ',plotTypes{plotI}])
        %set(gcf,'Name',[datasetName,' - Cross-fly Zero phase cross ',plotTypes{plotI}])
        set(gcf,'Name',['Cross-dataset Zero phase cross ',plotTypes{plotI}])
    end
    %}

    %Generalised across og and log
    subPlotTypes = {'allData','allLogData'};
    subPlotXLims = {'xLimVal','xLimLogVal'}; %Paired with above ofc
    subPlotMins = {'minVal','minLogVal'}; %Paired with above ofc
    subPlotNames = {'original','log'};

    for subPlotI = 1:size(subPlotTypes,2)
        figure
        
        %Quick pre-calculate of cross-dataset values
        thisMinVal = [];
        thisXLimVal = [];
        temp = [];
        temp2 = [];
        for fol = 1:size(overStruct,2)
            temp = [temp,overStruct(fol).plots.(plotTypes{plotI}).(subPlotMins{subPlotI})];
            temp2 = [temp2,overStruct(fol).plots.(plotTypes{plotI}).(subPlotXLims{subPlotI})];
        end
        thisMinVal = nanmin(temp);
        thisXLimVal = nanmax(temp2);

        for fol = 1:size(overStruct,2)
            thisPlotInd = sub2ind( [2,size(overStruct,2)] , 1, fol ); %Note: MATLAB arrays are linearly indicised vertically, but plots are indicised horizontally
            subplot(size(overStruct,2),2,thisPlotInd)
                %Note: The desired subplot behaviour is the datasets in col 1, arranged by number, then col 2 fullheight occupied by cross-dataset errorbarplot (Below)
    
            %thisData = overStruct(fol).plots.(plotTypes{plotI}).allData;
            thisData = overStruct(fol).plots.(plotTypes{plotI}).(subPlotTypes{subPlotI});
            allMean = nanmean(thisData);
            allMedian = nanmedian(thisData);
    
            %onourside
            %xLimVal = overStruct(fol).plots.(plotTypes{plotI}).xLimVal;
            %thisXLimVal = overStruct(fol).plots.(plotTypes{plotI}).(subPlotXLims{subPlotI}); %As below
            %theseBins = linspace(0,thisXLimVal,numBins);
            %thisMinVal = overStruct(fol).plots.(plotTypes{plotI}).(subPlotMins{subPlotI}); %Pre-defined above now
            theseBins = linspace(thisMinVal,thisXLimVal,numBins);
            histogram(thisData,'NumBins', numBins,'BinEdges',theseBins,'FaceColor',folColours(fol,:));
            hold on
            yLim = get(gca,'YLim');
            xLim = get(gca,'XLim');
        
            %Tack on mean/median calcs
            line([allMean,allMean],[0,nanmax(yLim)],'linewidth',2,'color','k','LineStyle','-'); 
            text(nanmax(xLim)*0.75,nanmax(yLim)*0.9,['mean=',num2str(allMean,3)],'FontSize',14);
            line([allMedian,allMedian],[0,nanmax(yLim)],'linewidth',2,'color','r','LineStyle','--'); 
            text(nanmax(xLim)*0.75,nanmax(yLim)*0.85,['median=',num2str(allMedian,3)],'FontSize',14);
        
            xlabel([plotTypes{plotI}])
            ylabel('Count')
            title([overStruct(fol).datasetName,' - All stim - n=',num2str(nansum(~isnan(thisData))),' zero phase ',plotTypes{plotI},' ',subPlotNames{subPlotI}])
            %set(gcf,'Name',[datasetName,' - Cross-fly Zero phase cross ',plotTypes{plotI}])
            set(gcf,'Name',['Cross-dataset Zero phase cross ',plotTypes{plotI},' ',subPlotNames{subPlotI}])
        end

        %Stats
        statsData = [];
        statsGroups = [];
        for fol = 1:size(overStruct,2)
            temp = overStruct(fol).plots.(plotTypes{plotI}).(subPlotTypes{subPlotI});
            temp( isnan(temp) ) = [];
            statsData = [statsData; temp];
            statsGroups = [statsGroups; repmat(fol,numel(temp),1)]; %Note: Will fail hard if data not n rows x 1 col
        end
        [P,ANOVATAB,STATS]=anova1([statsData],[statsGroups], 'off');
        statCompare = multcompare(STATS, 'CType', 'bonferroni', 'display', 'off');
        %Consider using barStats custom function here?

        %figure
        thisPlotInd = sub2ind( [2, size(overStruct,2)] , repmat(2,1,size(overStruct,2)), 1:size(overStruct,2) ); %Why this seems to be in col, row order i do not know
            %Note: Calcs not tested above 2x2 plot
        subplot(size(overStruct,2),2,thisPlotInd)
        hold on
        exLabels = {};
        for fol = 1:size(overStruct,2)
            temp = overStruct(fol).plots.(plotTypes{plotI}).(subPlotTypes{subPlotI});
            temp( isnan(temp) ) = [];
            thisMean = nanmean(temp);
            thisSEM = nanstd(temp) / sqrt(size(temp,1));
            errorbar(fol,thisMean,thisSEM,'LineWidth',2,'Color',folColours(fol,:))
            exLabels{fol} = overStruct(fol).datasetName;
        end
        xlim([0,size(overStruct,2)+1])
        ylim([thisMinVal,thisXLimVal]) %Inherit from histograms
        xticks([1:size(overStruct,2)])
        xticklabels(exLabels)
        titleStr = [overStruct(fol).datasetName,' - All stim - n=',num2str(nansum(~isnan(thisData))),' zero phase ',plotTypes{plotI},' ',subPlotNames{subPlotI},' stats'];
        titleStr = [titleStr,char(10),'p=',num2str(P)];
        title(titleStr)
        if size(overStruct,2) > 2
            disp(['-# Caution: One p-value reported on plot but number of comparisons was greater #-'])
            %I forget which ANOVA p is reported by default
        end

    %subPlotI end
    end

%plotI end
end

%e

%Comparative pos occupancy hist
tempX = [];
tempY = [];
for i=1:size(overStruct,2)
    tempX(i,:) = overStruct(i).specialPlots.posOccupancyHist.xLim;
    tempY(i,:) = overStruct(i).specialPlots.posOccupancyHist.yLim;
end


theseBins = linspace(nanmin(tempX(:,1)),nanmax(tempX(:,2)),64);

figure
for fol = 1:size(overStruct,2)
    subplot( size(overStruct,2) , 1, fol )

    thisData = overStruct(fol).specialPlots.posOccupancyHist.allData;
    allMean = nanmean(thisData);
    allMedian = nanmedian(thisData);

    histogram(thisData,'NumBins', 64,'BinEdges',theseBins,'FaceColor',folColours(fol,:));
    hold on
    yLim = get(gca,'YLim');
    xLim = get(gca,'XLim');

    %Tack on mean/median calcs
    line([allMean,allMean],[0,nanmax(yLim)],'linewidth',2,'color','k','LineStyle','-'); 
    text(nanmax(xLim)*0.75,nanmax(yLim)*0.9,['mean=',num2str(allMean,3)],'FontSize',14);
    line([allMedian,allMedian],[0,nanmax(yLim)],'linewidth',2,'color','r','LineStyle','--'); 
    text(nanmax(xLim)*0.75,nanmax(yLim)*0.85,['median=',num2str(allMedian,3)],'FontSize',14);

    xlabel(['Pos.'])
    ylabel('Count')
    title([overStruct(fol).datasetName,' - Pos occupancy hist - All stim '])
    set(gcf,'Name',['Cross-dataset pos occupancy hist'])
end


%--------------------------------------

%Mean transients
useMean = 0; %Whether to use mean (1) or median (0)

%Per-dataset
for fol = 1:size(overStruct,2)
    superTrans = nan(size(overStruct(fol).individuals,2), nanmax([overStruct(fol).individuals.transSize]), size(overStruct(fol).individuals(1).sepTrans,2));
    superTime = nan(size(overStruct(fol).individuals,2), nanmax([overStruct(fol).individuals.transSize]), size(overStruct(fol).individuals(1).sepTrans,2));
        %Row - Fly, Col - Time, Layer - Side
    for fly = 1:size(overStruct(fol).individuals,2)
        for sidi = 1:size(overStruct(fol).individuals(fly).sepTrans,2)
            thisMean = overStruct(fol).individuals(fly).sepTrans(sidi).thisMean;
            superTrans( fly, [1:size(thisMean,2)], sidi ) = thisMean;      
            thisTime = overStruct(fol).individuals(fly).sepTrans(sidi).thisTime;
            superTime( fly, [1:size(thisTime,2)], sidi ) = thisTime;     
        end
    end
    %hinkik

    %Interpolate, if necessary
        %Abandoned (For now)
    %{
    if interpolateMeanTransients
        %temp = reshape( superTime, size(superTime,1)*size(superTime,3), size(superTime,2) ); %Flatten to 2D matrix (Ignore sides); INCORRECT
        temp = [];
        for sidi = 1:size(superTime,3)
            temp = [temp; superTime(:,:,sidi)]; %Inefficient, but eh
        end
        medianTime = nanmedian(temp, 1);

        for i = 1:size(temp,1)
            norg= interp1(blirg, temp(i,:), medianTime);
        end
    end
    %}

    %Testatory time plot
    figure
    for sidi = 1:size(superTime,3)
        subplot(1,size(superTime,3),sidi)
        hold on
        for fly = 1:size(superTime,1)
            plot(superTime(fly,:,sidi))
        end
        xlabel('Frame')
        ylabel('Time')
        title(['Stim ending ',num2str(sidi)])
        set(gcf,'Name',[overStruct(fol).datasetName,' - Time comparison'])
    end

    %Plot mean trans per side for this dataset
    sideColours = jet(size(superTime,3));
    figure
    for sidi = 1:size(superTrans,3)
        subplot(1,size(superTrans,3),sidi)
        sidMean = squeeze( nanmean( superTrans(:,:,sidi), [1,3] ) );
        sidMedian = squeeze( nanmedian( superTrans(:,:,sidi), [1,3] ) );
        sidSEM = nanstd( superTrans(:,:,sidi), [], 1 ) / sqrt( size(superTrans,1) );
        sidTime = squeeze( nanmedian( superTime(:,:,sidi), [1,3] ) );
            %Note: Time is calculated as median, so theoretically may not be truly representative
        
        %Mean, Frame time ref
        if useMean
        plot( sidTime , sidMean, 'LineWidth', 2, 'Color', sideColours(sidi,:) )
        hold on
        groupShadeCoordsX = [sidTime,flip(sidTime)];
        groupShadeCoordsY = [sidMean+sidSEM,flip(sidMean-sidSEM)];
        else
        %Median, rel time ref
        plot( sidTime , sidMedian, 'LineWidth', 2, 'Color', sideColours(sidi,:) )
        hold on
        groupShadeCoordsX = [sidTime,flip(sidTime)]; %Note: Frame reference, not time
        groupShadeCoordsY = [sidMedian+sidSEM,flip(sidMedian-sidSEM)];
        end

        fill(groupShadeCoordsX, groupShadeCoordsY, sideColours(sidi,:)) %Error shading
        alpha(0.15)
        %for i = 1:size(thisData,1)
        %    plot(thisTime,thisData(i,:),'Color',transStruct(matchTrans(i)).colour)
        %    alpha(0.15)
        %end
        if useMean
        title(['Stim ending ',num2str(sidi),' (Mean)'])
        else
        title(['Stim ending ',num2str(sidi),' (Median)'])
        end
        %xlabel(['Frame no. [TIME NOT CONSISTENT ACROSS INDIVIDUALS]'])
        ylabel(['Pos. (degs)'])
        xlim([0,nanmax(sidTime)])
        xlabel(['Time (s)'])
    end
    if useMean
    set(gcf,'Name',[overStruct(fol).datasetName,' - Mean transient side sep'])
    else
    set(gcf,'Name',[overStruct(fol).datasetName,' - Median transient side sep'])
    end

    overStruct(fol).plots.meanTrans.superTrans = superTrans;
    overStruct(fol).plots.meanTrans.superTime = superTime;

end

clear sidMean sidMedian sidtime

%Comparative
if size(overStruct,2) > 1
    nSides = size(overStruct(1).plots.meanTrans.superTrans,3); %Use first as example

    folColours = lines(size(overStruct,2));

    %Comparative median transient
    for plotMode = 1:2 %1 - Sep by sides, 2 - Combined
        figure
        for sidi = 1:nSides
            if plotMode == 1
                subplot( 1, nSides, sidi )
            end
            legLabels = {};
            a = 1;
            sidTimeMax = [];
            for fol = 1:size(overStruct,2)
                sidMean = squeeze( nanmean( overStruct(fol).plots.meanTrans.superTrans(:,:,sidi), [1,3] ) );
                sidMedian = squeeze( nanmedian( overStruct(fol).plots.meanTrans.superTrans(:,:,sidi), [1,3] ) );
                sidSEM = nanstd( overStruct(fol).plots.meanTrans.superTrans(:,:,sidi), [], 1 ) / sqrt( size(overStruct(fol).plots.meanTrans.superTrans,1) );
                sidTime = squeeze( nanmedian( overStruct(fol).plots.meanTrans.superTime(:,:,sidi), [1,3] ) );
                sidTimeMax = [sidTimeMax,nanmax(sidTime)];
        
                if useMean
                plot( sidTime , sidMean, 'LineWidth', 2, 'Color', folColours(fol,:) )
                else
                plot( sidTime , sidMedian, 'LineWidth', 2, 'Color', folColours(fol,:) )
                end
                hold on
                groupShadeCoordsX = [sidTime,flip(sidTime)];
                if useMean
                groupShadeCoordsY = [sidMean+sidSEM,flip(sidMean-sidSEM)];
                else
                groupShadeCoordsY = [sidMedian+sidSEM,flip(sidMedian-sidSEM)];
                end
                            
                fill(groupShadeCoordsX, groupShadeCoordsY, folColours(fol,:) ) %Error shading
                alpha(0.15)
                legLabels{a} = overStruct(fol).datasetName;
                legLabels{a+1} = overStruct(fol).datasetName;
                a = a + 2;
            end
            ylabel('Position (degs)')
            %xlabel('Frame [TIME NOT CONSISTENT ACROSS INDIVIDUALS]')
            xlabel('Time (s)')
            legend(legLabels)
            legend('AutoUpdate','off')
            %line([0,size(sidMean,2)],[0,0],'LineStyle','--','Color','k')
            %line([sidTime(1),sidTime(end)],[0,0],'LineStyle','--','Color','k')
            line([0,4],[0,0],'LineStyle','--','Color','k') %4s
            %yLim = get(gca,'YLim');
            if plotMode == 1
            if sidi == 1
                ylim([-60,10]) %Hardcoded, because too lazy to make an overengineered system
            elseif sidi == 2
                ylim([-10,60])
            end
            else
            ylim([-70,70])
            end
    
            xlim([0,4]); %HARDCODED 4s X LIMIT
            if nanmax(sidTimeMax) > 4 && sidi == 1
                disp(['-# Caution: X-lim is 4s, but data exceeds this scope (',num2str(nanmax(sidTimeMax)),'s) #-'])
            end
    
        end
        if useMean
        nameStr = ['Comparative mean transient side sep'];
        else
        nameStr = ['Comparative median transient side sep'];
        end
        if plotMode == 1
            nameStr = [nameStr,'(Sep.)'];
        end
        set(gcf,'Name',nameStr);
    end



end


%--------------------------------------

%% End

if doAutomatedSavePlots
    figPath = baseFigPath; %Maybe in future make each dataset save to own separate folder?
    automatedSavePlots(figPath,'clearOldFigures',1)
    if alsoSaveOverStruct
        save([figPath,filesep,'overStruct.mat'],'overStruct')
        disp(['(overStruct also saved)'])
    end
end
