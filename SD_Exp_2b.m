% Sequential Dependancies for Moving Bar using Fixation Learning, then Open Loop and Closed Loop - SB
% Commenced 2/8/23 
% Runs a Seq Dependancies Experiment
% Uses UDP comms to get fly wing information (dWLE) and to send the experiment status and filename
% Sections  - User Input
%           - Parameters
%           - Set Up UDP, Screen, Stimulus
%           - Do Closed Loop Fixation - No SD
%           - End Fixation Start Seq Dep
%               - Qualification
%               - Start Seq Dep Open Loop
%               - Close the Loop and Measure SE time
%           - Save Files
%           - Plot Histograms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Real = 1;                   % 0 for testing - no user input for date fly exp
Plot = 0;                   % 1 for plot histograms when exp complete
screenNumber = 0;          % set to 1 for testing with 2nd monitor, set to 0 for rig %%%%%%%%%%%%%%%

clear uRx;                  % Close the UDP connection if open from before
clear uTx ;

if(Real)
    date = input('Date? ','s'); 
    fly = input('Fly? ','s');
    exp = input('Experiment? ','s');
    baseFilename = ['SD_', date, '_fly_', fly, '_exp_', exp];
else
    baseFilename = 'SD_Test';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Colour Param %%%%%%%%%%%
numberofParams = 24;             % Number of parameters to save
%global Param_Array; 
ParamArray = cell(numberofParams,2);
namesColour = {'colourBar', 'colourBack1', 'colourBack2'};
valuesColour = {'[0 0 255]','[128 128 128]','[128 128 128]'};
%%%%%%% General Param %%%%%%%%%%%
namesGen = {'intensityBar'  ,'intensityBack1'   ,'intensityBack2'   ,'widthBar' ,'widthBack1'   ,'gainProp' ,'angFriction'};
valuesGen = [0             ,1                  ,1                  ,20          ,16             ,0.002       ,.1];
%            0 to 1         0 to 1              0 to 1              pixels
namesFix = {'posStartFix'   ,'timeFix'  ,'trialsFix'};
valuesFix = [0             ,20          ,6];

namesSD = {'posStartSE'     ,'timeOutQual'  ,'qualWindow'   ,'qualTime' ,'qualTarget'   ,'timeOnOL' ,'timeOffOL','timeOnCL' ,'timeOffCL','flashDiv' ,'targetFix'    ,'trialsSE'};
valuesSD = [60              ,200            ,32             ,1000       ,5.33           ,0.05      ,0.1        ,4          ,0          ,0           ,16            ,120 ];
%           degrees         secs            pixels          ms          pixels          secs        secs        secs        not used    not used
Inc = 0;
for i = 1:length(namesColour)
    eval([namesColour{i}, ' = ', valuesColour{i}, ';']);
    ParamArray{i+Inc,1} = namesColour{i};
    ParamArray{i+Inc,2} = string(valuesColour(i));
end
Inc = length(namesColour);
for i = 1:length(namesGen)
    eval([namesGen{i}, ' = ', num2str(valuesGen(i)), ';']);
    ParamArray{i+Inc,1} = namesGen{i};
    ParamArray{i+Inc,2} = valuesGen(i);
end
Inc = Inc + length(namesGen);
for i = 1:length(namesFix)
    eval([namesFix{i}, ' = ', num2str(valuesFix(i)), ';']);
    ParamArray{i+Inc,1} = namesFix{i};
    ParamArray{i+Inc,2} = valuesFix(i);
end
Inc = Inc + length(namesFix);
for i = 1:length(namesSD)
    eval([namesSD{i}, ' = ', num2str(valuesSD(i)), ';']);
    ParamArray{i+Inc,1} = namesSD{i};
    ParamArray{i+Inc,2} = valuesSD(i);
end
Inc = Inc + length(namesSD);
Filename = mfilename('fullpath');
ParamArray{Inc,1} = 'Filename';
ParamArray{Inc,2} = Filename;  
%setParam(namesGen,valuesGen,2);
%setParam(namesFix,valuesFix,7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set Up UDP port(s), Screen and Stimulus %%%%%%%
UDPportNumber = 65438  ;          % UDP port number matching Python UDP port number
uRx = udpport("datagram",'LocalPort',UDPportNumber);
uTx = udpport("datagram");
write(uTx,baseFilename ,"127.0.0.1",UDPportNumber+1);

%%%%%% just set SupressAllWarnings to 0 to see the outputs
Screen('Preference', 'SuppressAllWarnings', 1);
Screen('Preference', 'SkipSyncTests', 1); % set it to '0' to see all of the warnings
Screen('Preference','VisualDebugLevel', 0); % can use values 0 to 5
%priorityLevel=MaxPriority(w); % Set maximum useable priority Level on this system

try
    AssertOpenGL;
    [w, screenRect]=Screen('OpenWindow',screenNumber, [0 0 0]); % Open a double buffered fullscreen window and set default background
    priorityLevel=MaxPriority(w); % Set maximum useable priority Level on this system
    
    %%%%%%%%%%%%%%%%% Stimulus Arrays &&&&&&&&&&&&&&&&&&&&&&&&&
    stimArrayMid = zeros(1,2571,3);        % create 1-D array
    stimArrayLeft = zeros(1,2571,3);        % create 1-D array
    stimArrayRight = zeros(1,2571,3);        % create 1-D array
    stimArrayBlank = zeros(1,2571,3);        % create 1-D array
    pixZero = 270*192/360;         % straight ahead
    pixOffset = posStartSE*192/360;  
    bar = zeros(1,1,3); bar(1,1,:) = intensityBar * colourBar;
    bkgd1 = zeros(1,1,3); bkgd1(1,1,:) = intensityBack1 * colourBack1;
    bkgd2 = zeros(1,1,3); bkgd2(1,1,:) = intensityBack2 * colourBack2;
    for i = 1:2571                  % Fill 1-D array
        if (mod(i,192) > (pixZero - widthBar/2) && mod(i,192) < (pixZero + widthBar/2)) % Modulo so that it wrapps around 0 to 192
            stimArrayMid(1,i,:) = bar;%intensity_bar;
            stimArrayLeft(1,i,:) = bkgd1;%intensity_background;
            stimArrayRight(1,i,:) = bkgd1;%intensity_background;
            stimArrayBlank(1,i,:) = bkgd1;%intensity_background;
        elseif(mod(i,192) > (pixZero - widthBar/2 + pixOffset) && mod(i,192) < (pixZero + widthBar/2 + pixOffset)) % Modulo so that it wrapps around 0 to 192
            stimArrayLeft(1,i,:) = bkgd1;%intensity_background;
            stimArrayMid(1,i,:) = bkgd1;%intensity_background;
            stimArrayRight(1,i,:) = bar;%intensity_bar;
            stimArrayBlank(1,i,:) = bkgd1;%intensity_background;
        elseif(mod(i,192) > (pixZero - widthBar/2 - pixOffset) && mod(i,192) < (pixZero + widthBar/2 - pixOffset)) % Modulo so that it wrapps around 0 to 192
            stimArrayLeft(1,i,:) = bar;%intensity_bar;
            stimArrayMid(1,i,:) = bkgd1;%intensity_background;
            stimArrayRight(1,i,:) = bkgd1;%intensity_background;
            stimArrayBlank(1,i,:) = bkgd1;%intensity_background;
        elseif(mod(i,192) < (pixZero - widthBack1/2) && mod(i,192) > (pixZero - 97))
            stimArrayMid(1,i,:) = bkgd2 ;%intensity_background;
            stimArrayLeft(1,i,:) = bkgd2;%intensity_background; 
            stimArrayRight(1,i,:) = bkgd2;%intensity_background;
            stimArrayBlank(1,i,:) = bkgd2;%intensity_background;
            % disp(i); 
        elseif(mod(i,192) > mod((pixZero + widthBack1/2),192)  && mod(i,192) < mod((pixZero + 97),192))
            stimArrayMid(1,i,:) = bkgd2 ;%intensity_background;          
            stimArrayLeft(1,i,:) = bkgd2;%intensity_background;
            stimArrayRight(1,i,:) = bkgd2;%intensity_background;
            stimArrayBlank(1,i,:) = bkgd2;%intensity_background;
            %disp(i); 
        else 
            stimArrayMid(1,i,:) = bkgd1;%intensity_background;
            stimArrayLeft(1,i,:) = bkgd1;%intensity_background;
            stimArrayRight(1,i,:) = bkgd1;%intensity_background;
            stimArrayBlank(1,i,:) = bkgd1;%intensity_background;
        end
    end

    stimTexMid=Screen('MakeTexture', w, stimArrayMid);   % Store 1-D single row grating in texture:
    stimTexLeft=Screen('MakeTexture', w, stimArrayLeft);   % Store 1-D single row grating in texture:
    stimTexRight=Screen('MakeTexture', w, stimArrayRight);   % Store 1-D single row grating in texture:
    stimTexBlank=Screen('MakeTexture', w, stimArrayBlank);   % Store 1-D single row grating in texture:
    width = 1920;               % Not sure why but it works;
    height = 1080;              % Not sure why but it works;
    dstRect=[0 0 width height];
    dstRect=CenterRect(dstRect, screenRect);
    %%%%%%%%%%%%%%%%%%%%%%%%%% Set waitframes %%%%%%%%%%%%%%%%%%%%%%%
    ifi=Screen('GetFlipInterval', w);       % Query duration of one monitor refresh interval:
    Display  = ['refresh rate', ifi];
    disp(Display);
    waitframes = 1;             % increase if GPU struggles
    waitduration = waitframes * ifi;    % Translate frames into seconds for screen update interval:
    %%%%%%%%%%%%%%%%%%%%%%%% Set up Parameters %%%%%%%%%%%%%%%%%%%%
    FixArray = zeros(720000, 6);      % 1 hour at 200Hz, time, lastUDPtime, lastWLE, stim degrees, on front, on front last 10s
    SEArray = zeros(720000, 8);      % 1 hour at 200Hz, time, lastUDPtime, lastWLE, stim degrees, Ang
    HistArray = zeros(720000, 1);      % For Histogram
    resultsArray = zeros(trialsSE,13);  % results - timeQual, posOL[1:5], timeSE
    %Save_Array = [];
    dWLE = 0;
    RXtime = 0;
    Maxdelay = 0; 
    loops = 0;
    noUDP = 0;
    DataIndexFix = 1;
    DataIndexSE = 1;
    AngAcc = 0;
    FixAv = 1000;
    Flying = 0;
    Flash = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Do Closed Loop Fixation - No SD %%%%%%%%%%%%%
    xoffset = posStartFix*192/360;
    for trial = 1:trialsFix
        flush(uRx);
        % *************************Start bar at initial pos and move it according to WLE for time_move secs unless abort via keypress. 
        vbl=Screen('Flip', w);
        vtimeMove = vbl + timeFix;
        while vbl < vtimeMove
            [Maxdelay, Flying, noUDP, pMarktime, dWLE, RXtime] = CheckUDP(Maxdelay, noUDP, Flying, dWLE, RXtime, uRx); %%%%%% check UDP port %%%
            [xoffset, AngAcc] = calcOffset(xoffset, AngAcc, dWLE, gainProp, angFriction);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Draw screen %%%%%%%%%%%%%
            srcRect=[xoffset 0 xoffset + width height];
            if(Flash == 1)
                Screen('DrawTexture', w, stimTexMid, srcRect, dstRect, 0,[],[]);
                %disp('on');
            else
                Screen('DrawTexture', w, stimTexBlank, srcRect, dstRect, 0,[],[]);
                %disp('off');
            end
            vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);
            if(flashDiv >0)
                if(mod(DataIndexFix,flashDiv) == 0)
                    if(Flash == 1)
                        Flash =0;
                    else
                        Flash =1;
                    end
                end
            end 
            if KbCheck
                break;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%% Sort out Parmaters and Arrays %%%
            if(loops == 0)  
                TimeLast = pMarktime;
                lastPrintTime = pMarktime;
                MaxLooptime = 0;
            end
            Looptime = pMarktime - TimeLast;
            TimeLast = pMarktime;
            if MaxLooptime < Looptime  
                MaxLooptime = Looptime;
            end
            loops = loops + 1;
            %if(Flying == 1)
            if(1)
                if (xoffset > -17 && xoffset < 17)
                    inFront = 1;
                    FixAv = FixAv + 1;
                else
                    inFront = 0;
                    FixAv = FixAv - 1;
                end
                if FixAv > 2000
                    FixAv = 2000;
                elseif FixAv <0
                    FixAv = 0;
                end
                FixArray(DataIndexFix,:) = [RXtime, pMarktime, dWLE, xoffset*360/192, inFront, FixAv];
                HistArray(DataIndexFix) = xoffset*360/192;
                DataIndexFix = DataIndexFix + 1;
            end
        end
        %%%%%%%%%%%%%%%%% End Fixation Trial %%%%%%%%%%%%%%%%%%
        fprintf('Block %d, Av Fix %d, Max UDP %.f, Max Loop %.1f, Av Loop %.2f no UDP %d \n', trial, FixAv, Maxdelay, MaxLooptime, timeFix*1000/loops, noUDP);    
        Maxdelay = 0; 
        loops = 0;
        noUDP = 0;
        FixAv = 1000;

        
        write(uTx,baseFilename ,"127.0.0.1",UDPportNumber+1);
        if KbCheck
            break;
        end
    end 
    disp( ' ' );
    
    %%%%%%%%%%%%%%%%% End Fixation Start Seq Dep %%%%%%%%%%%%%%%%%%%%%
    for trial = 1:trialsSE
        posOLbin = randi([0,1],1,5);    %%%%%%%%%%%%%%% Random Sequence for SD %%%%%%%%%
        %disp(posOLbin);
        posOL = posOLbin*2-1;           %%%%%%%%%%%% change from 0, to -1, 1 %%%%%%%%
        fprintf('Trial %d  ',trial);
        disp(posOL);
        posOLdec = posOLbin(1)+posOLbin(2)*2+posOLbin(3)*4+posOLbin(4)*8+posOLbin(5)*16;    %%% Decimal of the sequence %%%%
        posOLDiff = abs(diff(posOLbin)) ;
        posSD = posOLDiff(1)+posOLDiff(2)*2+posOLDiff(3)*4+posOLDiff(4)*8+1;        %%%%%%% SD position (1 to 16) %%%%%%
        %%%%%%%%%%%%%%%%%%%%%% Qualification %%%%%%%%%%%%%%%%%%%%%%%%%%
        flush(uRx);
        % *************************Start bar at initial pos and move it according to WLE for time_move secs unless abort via keypress. 
        vbl=Screen('Flip', w);
        vtimeMove = vbl + timeOutQual;
        timeStartQual = getBrisTime;
        timeUpdate = timeStartQual;
        timeQualWindow = timeStartQual;
        Attention = 0;
        while (vbl < vtimeMove) && (Attention == 0)
            [Maxdelay, Flying, noUDP, pMarktime, dWLE, RXtime] = CheckUDP(Maxdelay, noUDP, Flying, dWLE, RXtime, uRx); %%%%%% check UDP port %%%
            [xoffset, AngAcc] = calcOffset(xoffset, AngAcc, dWLE, gainProp, angFriction);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Draw screen %%%%%%%%%%%%%
            srcRect=[xoffset 0 xoffset + width height];
            if(Flash == 1)
                Screen('DrawTexture', w, stimTexMid, srcRect, dstRect, 0,[],[]);
                %disp('on');
            else
                Screen('DrawTexture', w, stimTexBlank, srcRect, dstRect, 0,[],[]);
                %disp('off');
            end
            vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);
            if(flashDiv >0)
                if(mod(DataIndexSE,flashDiv) == 0)
                    if(Flash == 1)
                        Flash =0;
                    else
                        Flash =1;
                    end
                end
            end 
            if KbCheck
                break;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%% Sort out Parmaters and Arrays %%%
            if(loops == 0)  
                TimeLast = pMarktime;
                lastPrintTime = pMarktime;
                MaxLooptime = 0;
            end
            Looptime = pMarktime - TimeLast;
            TimeLast = pMarktime;
            if MaxLooptime < Looptime  
                MaxLooptime = Looptime;
            end
            loops = loops + 1;
            %if(Flying == 1)
            if(1)
                if((abs(xoffset) > (qualWindow/2)) || (Flying == 0))
                    timeQualWindow = pMarktime;
                else
                    if((mod(loops,50) == 0)&& (Flying == 1))
                        fprintf('in qual window  %.f \n',pMarktime - timeQualWindow);
                    end
                end
                SEArray(DataIndexSE,:) = [RXtime, pMarktime, dWLE, xoffset*360/192, AngAcc, trial, 1, posOLdec];
                DataIndexSE = DataIndexSE + 1;
                if((abs(xoffset) < qualTarget) && (pMarktime > (timeQualWindow + qualTime)))
                    windowTime = pMarktime - timeQualWindow;
                    Attention = 1;  %%%%%%%%%%%%%%%%%%%%%%%%%% Fly paying attention %%%%%%%%%%%%%%%
                end
            end
            timeQual = pMarktime - timeStartQual;
            if(loops == 1000)
                fprintf('Qual %.1f, time %d, Max UDP %.f, Max Loop %.1f, Av Loop %.2f no UDP %d \n', trial, timeQual, Maxdelay, MaxLooptime, (pMarktime-timeUpdate)/loops, noUDP);    
                Maxdelay = 0; 
                loops = 0;
                noUDP = 0;
                timeUpdate = pMarktime;
                %disp( ' ' );
                write(uTx,baseFilename ,"127.0.0.1",UDPportNumber+1);
            end
        end
        posQual = SEArray(DataIndexSE-1,4);
        posPreQual = SEArray(DataIndexSE-11,4);
        fprintf('window time %.1f, Pos -50ms %.1f, Pos Qual %.1f\n',windowTime,posPreQual,posQual);
        %%%%%%%%%%%%%%%%% End Qualification %%%%%%%%%%%%%%%%%%
        timeQual = pMarktime - timeStartQual;
        fprintf('Qual %d, time %.1f, Max UDP %.f, Max Loop %.1f, Av Loop %.2f no UDP %d \n', trial, timeQual, Maxdelay, MaxLooptime, (pMarktime-timeUpdate)/loops, noUDP);    
        Maxdelay = 0; 
        loops = 0;
        noUDP = 0;
        %disp(pMarktime - timeQualWindow);
        disp( ' ' );
        write(uTx,baseFilename ,"127.0.0.1",UDPportNumber+1);
        %%%%%%%%%%%%%%%%%%%%% Start Seq Dep Open Loop %%%%%%%%%%%%%%%%%%%
        AngAcc = 0;
        for i = 1:4
            vbl=Screen('Flip', w);
            vtime_disp = vbl + timeOnOL;
            if(posOL(i) == 1)
                xoffset = posStartSE*192/360;
            else
                xoffset = -posStartSE*192/360;
            end
            while vbl < vtime_disp
                [Maxdelay, Flying, noUDP, pMarktime, dWLE, RXtime] = CheckUDP(Maxdelay, noUDP, Flying, dWLE, RXtime, uRx); %%%%%% check UDP port %%%
                %[xoffset, AngAcc] = calcOffset(xoffset, AngAcc, dWLE, gainProp, angFriction);
                srcRect=[xoffset 0 xoffset + width height];
                Screen('DrawTexture', w, stimTexMid, srcRect, dstRect, 0,[],[]);
                vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);
                pMarktime = getBrisTime;
                SEArray(DataIndexSE,:) = [RXtime, pMarktime, dWLE, xoffset*360/192, AngAcc, trial, 2, posOLdec];
                DataIndexSE = DataIndexSE +1;
                if KbCheck
                    break;
                end
            end
            vbl=Screen('Flip', w);
            vtime_disp = vbl + timeOffOL;
            while vbl < vtime_disp
                [Maxdelay, Flying, noUDP, pMarktime, dWLE, RXtime] = CheckUDP(Maxdelay, noUDP, Flying, dWLE, RXtime, uRx); %%%%%% check UDP port %%%
                %[xoffset, AngAcc] = calcOffset(xoffset, AngAcc, dWLE, gainProp, angFriction);
                srcRect=[0 0 0 + width height];
                Screen('DrawTexture', w, stimTexBlank, srcRect, dstRect, 0,[],[]);
                vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);
                pMarktime = getBrisTime;
                SEArray(DataIndexSE,:) = [RXtime, pMarktime, dWLE, 200, AngAcc, trial, 2, posOLdec];
                DataIndexSE = DataIndexSE + 1;
                if KbCheck
                    break;
                end
            end
        end
        %%%%%%%%%%%%%%% Close the Loop and Measure SE time %%%%%%%%%%%%%%%%%
        flush(uRx);
        if(posOL(5) == 1)
            xoffset = posStartSE*192/360;
        else
            xoffset = -posStartSE*192/360;
        end
        timeOffCL = 0;
        vbl=Screen('Flip', w);
        vtimeMove = vbl + timeOnCL;
        timeStartCL = getBrisTime;
        %reFix = 0;
        timeReFix = -1000;
        while (vbl < vtimeMove) %&& (reFix ==0)
            [Maxdelay, Flying, noUDP, pMarktime, dWLE, RXtime] = CheckUDP(Maxdelay, noUDP, Flying, dWLE, RXtime, uRx);
            [xoffset, AngAcc] = calcOffset(xoffset, AngAcc, dWLE, gainProp, angFriction);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Draw screen %%%%%%%%%%%%%
            srcRect=[xoffset 0 xoffset + width height];
            if(Flash == 1)
                Screen('DrawTexture', w, stimTexMid, srcRect, dstRect, 0,[],[]);
                %disp('on');
            else
                Screen('DrawTexture', w, stimTexBlank, srcRect, dstRect, 0,[],[]);
                %disp('off');
            end
            vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);
            if(flashDiv >0)
                if(mod(DataIndexSE,flashDiv) == 0)
                    if(Flash == 1)
                        Flash =0;
                    else
                        Flash =1;
                    end
                end
            end 
            if KbCheck
                timeReFix = -3000;
                break;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%% Sort out Parmaters and Arrays %%%
            if(loops == 0)  
                TimeLast = pMarktime;
                lastPrintTime = pMarktime;
                MaxLooptime = 0;
            end
            Looptime = pMarktime - TimeLast;
            TimeLast = pMarktime;
            if MaxLooptime < Looptime  
                MaxLooptime = Looptime;
            end
            loops = loops + 1;
            %pMarktime = getBrisTime;
            SEArray(DataIndexSE,:) = [RXtime, pMarktime, dWLE, xoffset*360/192, AngAcc, trial, 3, posOLdec];
            DataIndexSE = DataIndexSE + 1;
            if((abs(xoffset) > 80) && (timeReFix == -1000))
                timeReFix = -2000;        %%%%%%%%%%%%%%%%%%%% Round the Back %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                disp('round the back');
            end
            if((abs(xoffset) < targetFix) && (timeReFix == -1000))   %%%%%%%%%% Refixated on Bar %%%%%%%%%%
                %reFix = 1;
                timeReFix = pMarktime - timeStartCL;
                disp(posOL);
                disp('refixated');
                disp(timeReFix);
            end
        end
        %if (reFix == 0)
        %    timeReFix = -10;    %%%%%%%%%%%%%%%%%%%% Timeout %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %end
        %%%%%%%%%%%%%%%%% End Seq Dep Trial %%%%%%%%%%%%%%%%%%
        fprintf('SE Trial %d, Max UDP %.f, Max Loop %.1f, Av Loop %.2f no UDP %d \n', trial, Maxdelay, MaxLooptime, timeOnCL*1000/loops, noUDP);    
        %fprintf('time %f, loops %d', vbl-timeStartCL, loops);
        Maxdelay = 0; 
        loops = 0;
        noUDP = 0;
        resultsArray(trial, 1) = pMarktime;
        resultsArray(trial, 2) = timeQual;
        for i = 1:5
            resultsArray(trial, i+2) = posOL(i);
        end
        resultsArray(trial, 8) = posOLdec;
        resultsArray(trial, 9) = posSD;
        resultsArray(trial, 10) = timeReFix;
        resultsArray(trial, 11) = posPreQual;
        resultsArray(trial, 12) = posQual;
        resultsArray(trial, 13) = windowTime;
        disp(timeReFix);
        write(uTx,baseFilename ,"127.0.0.1",UDPportNumber+1);
        if KbCheck
            break;
        end
    end 
    %%%%%%%%%%%%%%%%% End Seq Dep %%%%%%%%%%%%%%%%%%%%%
    
    Priority(0);    % Restore normal priority
    sca;    %close onscreen and offscreen windows also close textures.
    %writematrix(Save_Array,FilenameCSV);
catch
    % Closes the onscreen window if its open on error.
    sca;
    Priority(0);
    psychrethrow(psychlasterror);
end %try..catch..
write(uTx,"Stop" ,"127.0.0.1",UDPportNumber+1);
write(uTx,"Stop" ,"127.0.0.1",UDPportNumber+1);
clear uRx;                    % Close the UDP connection
clear uTx ;

%%%%%%%%%%%%% Save Files %%%%%%%%%%%%%%%%%%%%
Save_Array_Fix = FixArray(1:DataIndexFix-1,1:6);
Save_Array_SeqDep = SEArray(1:DataIndexSE-1,1:8);
Save_Array_Results = resultsArray(1:trial,1:13);
disp('trial');
Folder = 'C:\Experiments\';
fixFilenameCSV = [Folder,'ML_Fix_',baseFilename,'.csv'];
seqdepFilenameCSV = [Folder,'ML_Data_',baseFilename,'.csv'];
resultsFilenameCSV = [Folder,'ML_Results_',baseFilename,'.csv'];
paramFilenameCSV = [Folder,'ML_Par_',baseFilename,'.csv'];
%stimulusFilenameMAT = [Folder,'ML_',baseFilename,'.mat'];
writematrix(Save_Array_Fix,fixFilenameCSV);
writematrix(Save_Array_SeqDep,seqdepFilenameCSV);
writematrix(Save_Array_Results,resultsFilenameCSV);
myTable = cell2table(ParamArray);
writetable(myTable, paramFilenameCSV, 'Delimiter', ',');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Histograms %%%%%%%%%%
if (Plot)
    Hist_Index = floor(DataIndexFix/4)*4;
    disp(DataIndexFix);
    disp(Hist_Index);
    if (Hist_Index > 20)
        Hist_Array1 = HistArray(1:Hist_Index/4);
        Hist_Array2 = HistArray(Hist_Index/4:Hist_Index/2);
        Hist_Array3 = HistArray(Hist_Index/2:Hist_Index*3/4);
        Hist_Array4 = HistArray(Hist_Index*3/4:Hist_Index);

        binEdges = [-180:3:180];
        subplot(2,2,1);
        histogram(Hist_Array1,'BinEdges', binEdges);
        title('Position Histogram1');
        xlabel('Degrees');
        ylabel('Frequency');
        subplot(2,2,2);
        histogram(Hist_Array2,'BinEdges', binEdges);
        title('Position Histogram2');
        xlabel('Degrees');
        ylabel('Frequency');
        subplot(2,2,3);
        histogram(Hist_Array3,'BinEdges', binEdges);
        title('Position Histogram3');
        xlabel('Degrees');
        ylabel('Frequency');
        subplot(2,2,4);
        histogram(Hist_Array4,'BinEdges', binEdges);
        title('Position Histogram4');
        xlabel('Degrees');
        ylabel('Frequency');
    end 
end
% We're done!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calc Xoffset AngAcc from dWLE using gain and friction parameters 
function [xoffset, AngAcc] = calcOffset(xoffsetIn, AngAccIn, dWLE, gainProp, angFriction)
    AngAcc = AngAccIn*(1-angFriction) + gainProp*dWLE;
    xoffset = xoffsetIn + AngAcc; 
    if(xoffset > 96)
        xoffset = xoffset - 192;        %Clip to +/-180deg
    end
    if(xoffset < -96)
        xoffset = xoffset + 192;
    end
end

% get the Posixtime for Brisbane
function timeBrisbane = getBrisTime()
    Marktime = datetime('now');
    Marktime.TimeZone = "Australia/Brisbane";
    timeBrisbane = posixtime(Marktime)*1000;
end

% Check the UDP from Python and update parametersaccordingly
function [Maxdelay, Flying, noUDP, pMarktime, dWLE, RXtime] = CheckUDP(Maxdelayin, noUDPin, Flyingin, dWLEin, RXtimein, u)
    noUDP = noUDPin;
    Maxdelay = Maxdelayin;
    Flying = Flyingin;
    dWLE = dWLEin;
    RXtime = RXtimein;
    if(u.NumDatagramsAvailable > 0)             % check for available UDP datagrams from Python
        datagram = read(u,1);                   % read datagram
        pMarktime = getBrisTime;
        RXnumber = str2double(char(datagram.Data));
        if(RXnumber > 500)
            RXtime = RXnumber;
            Delayms = (pMarktime - RXtime);
            %fprintf(' RX time %d, Mark time %d, delay %d', RXtime, pMarktime, Delayms);
            if(Maxdelayin < Delayms)
                Maxdelay = Delayms;
            end  
        elseif(RXnumber <360) 
            dWLE = RXnumber;
            Flying = 1;
            %disp(dWLE);
        else
            dWLE = 0;
            Flying = 0;
            %disp(dWLE);
        end
    else
        pMarktime = getBrisTime;
        noUDP = noUDP + 1;
    end
end