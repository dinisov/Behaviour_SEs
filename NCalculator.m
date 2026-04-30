%Script for calculating effective N from folder
clear
close all

%%

%folder = "C:\Users\uqmvan13\ScottCode\Data2"
folder = "C:\Users\uqmvan13\ScottCode\LaraRecapitulation\Split\Control"

%%

allFiles = dir([char(folder),filesep,'*.csv']);

totalFileCount = size(allFiles,1);


flyDeets = struct;
for subN = 1:totalFileCount
    thisFileName = allFiles(subN).name;
    %blirg
    disp(['File #',num2str(subN),' - ',thisFileName])
    %Initial name cleaning
    thisFileNameProc = thisFileName;
    thisFileNameProc = strrep(thisFileNameProc, ' ','');
    thisFileNameProc = strrep(thisFileNameProc, '.csv','');
    thisFileNameProc = split(thisFileNameProc,'_');
    %Identify elements
    prosDate = [];
    prosFly = [];
    prosExp = [];
    for fragN = 1:size(thisFileNameProc,1)
        temp = thisFileNameProc{fragN};
        %Date
            %Note: For Scott data, probably include a backup dateModified system here
        if size(temp,2) == 6 || size(temp,2) == 5 || size(temp,2) == 4
            try
                prosDate = str2num(temp); %Technically converting a date into a 6-digit non-date number, but eh
                %prosDate = temp; %Switch to string, in case of proximal 0 etc
            catch
                ['-# Caution: Error converting potential date to number #-']
            end
        end
        %Fly
        if isequal(temp,'fly') && fragN < size(thisFileNameProc,1)
            try
                prosFly = str2num(thisFileNameProc{fragN+1}); %Technically converting a date into a 6-digit non-date number, but eh
            catch
                ['-# Caution: Error obtaining fly number #-']
            end
        end
        %Exp
        if isequal(temp,'exp') && fragN < size(thisFileNameProc,1)
            try
                prosExp = str2num(thisFileNameProc{fragN+1}); %Technically converting a date into a 6-digit non-date number, but eh
            catch
                ['-# Caution: Error obtaining exp number #-']
            end
        end

        %Backup for data that lacks a date in the name
        if isempty(prosDate)
            %antenora
            try
                prosDate = str2num( datestr(datetime( allFiles(subN).date ),'ddmmyy') );
                disp(['(Date of exp had to be derived from dateModified)'])
            catch
                ['-# Caution: Error deriving date from dateModified #-']
            end
        end

    end
    disp(['Derived date: ',num2str(prosDate),', fly: ',num2str(prosFly),', exp: ',num2str(prosExp)])
    %disp(['Derived date: ',prosDate,', fly: ',num2str(prosFly),', exp: ',num2str(prosExp)]) 
    flyDeets(subN).date = prosDate;
    flyDeets(subN).fly = prosFly;
    flyDeets(subN).exp = prosExp;
end

%--------------------

disp(['Apparent unique dates: ',num2str(unique([flyDeets.date]))])

dateIndex = unique([flyDeets.date]);
flyIndex = unique([flyDeets.fly]);
expIndex = unique([flyDeets.exp]);

%flyColl = cell(2,numel(dateIndex));
%flyColl = cell(1,99);
flyColl = zeros(nanmax(dateIndex),99,99); %Highly inefficient, but eh
%QA
%if numel(dateIndex) >= 99
%    ['## Alert: Large number of prospective dates detected; Analysis methods unprepared for situation ##']
%    crash = yes %If this happens, increase arbitrary values I guess?
%end
if nanmax(flyIndex) >= 99 || nanmax(expIndex) >= 99
    ['## Alert: Large number of prospective flies or exp # detected; Analysis methods unprepared for situation ##']
    crash = yes %If this happens, increase arbitrary values I guess?
end
for subN = 1:size(flyDeets,2)
    thisDateInd = dateIndex( find(dateIndex == flyDeets(subN).date) );
    if isempty(thisDateInd)
        disp(['-# Date was not apparently derivable for dataset #',num2str(subN),' #-'])
        thisDateInd = 99;
    end
    if numel(thisDateInd) > 1
        ['-# Alert: Date overfind #-']
        crash = yes
    end

    if ~isempty(flyDeets(subN).fly)
        thisFlyInd = flyDeets(subN).fly;
    else
        thisFlyInd = 99; 
        disp(['-# Caution: Fly # was not apparently derivable for dataset #',num2str(subN),'; Arbitrary number assigned #-'])
    end

    if ~isempty(flyDeets(subN).exp)
        thisExpInd = flyDeets(subN).exp;
    else
        thisExpInd = 99; 
        disp(['-# Caution: Fly # was not apparently derivable for dataset #',num2str(subN),'; Arbitrary number assigned #-'])
    end

    %Assign
    flyColl(thisDateInd,thisFlyInd,thisExpInd) = flyColl(thisDateInd,thisFlyInd,thisExpInd) + 1;


end

%------

disp(char(10))

prosUniqueFlyNum = size(flyDeets,2);

flatColl = nansum(flyColl,[3]); %Flatten exp i.e. bring all exp to surface in columns

dupls = nansum( flatColl == 2, 'all' );
subtractionValue = 0;
if dupls > 0
    disp([num2str(dupls),' apparent fly extra representation/s (i.e. More than one exp) detected']);
    temp = find( flatColl == 2 );
    for du = 1:numel(temp)
        [i,j,k] = ind2sub( size(flatColl), temp(du)); %Find coords of duplications
        temp2 = find( flyColl(i,j,:) == 1 );
        disp([num2str(i),' - fly ',num2str(j),' - exps ',num2str(temp2')])
        subtractionValue = subtractionValue + (numel(temp2)-1);
    end
    disp(['Original unique count (',num2str(prosUniqueFlyNum),') will be subtracted by number of duplicates (',num2str(subtractionValue),')'])
        %Note here that 'duplicate' just means one fly had multiple experiments included
end

actualN = prosUniqueFlyNum - subtractionValue;

%---

disp([char(10),'Calculated unique N =',num2str(actualN)])




