%Script for matching two arbitrary sets of CSVs
%i.e. "Data" and "Fix" etc
clear
close all

%%

sourceFolder = 'C:\Users\uqmvan13\ANALYSIS\Scott\Data\Scott\Data2';
matchFolder = 'C:\Users\uqmvan13\ANALYSIS\Scott\Data\Scott\FIX\ALL';

sinkFolder = 'C:\Users\uqmvan13\ANALYSIS\Scott\Data\Scott\FIX\Data2';

%%

fileType = '.csv';

sourceDescriptorDeletion = 'Data'; %This will be deleted in the acquired names during search
matchDescriptorDeletion = 'Fix'; %Ditto, match folder

actuallyCopy = 1; %Whether to enact a copy (or not)

%%

sourceFiles = dir([sourceFolder,filesep,'*',fileType]);
disp([num2str(size(sourceFiles,1)),' files found in source folder'])

matchFiles = dir([matchFolder,filesep,'*',fileType]);
disp([num2str(size(matchFiles,1)),' files found in match folder'])

if isempty(sourceFiles) || isempty(matchFiles)
    ['## Alert: One or both folders did not contain valid files ##']
    crash = yes
end

%%

%Match
matchYes = 0;
matchNo = 0;
for file = 1:size(sourceFiles,1)
    thisFile = sourceFiles(file).name;
    disp(['Finding match for ',thisFile])

    thisFile = strrep(thisFile,sourceDescriptorDeletion,''); %It doesn't matter this isn't a valid name; It will be matched with another invalid name

    matches = [];
    for subFile = 1:size(matchFiles,1)
        thisMatchFile = matchFiles(subFile).name;
        thisMatchFile = strrep(thisMatchFile,matchDescriptorDeletion,'');
        if isequal(thisFile,thisMatchFile)
            disp(['Match found (',matchFiles(subFile).name,')'])
            matches = [matches,subFile];
        end
    end
    if isempty(matches)
        disp(['No match found'])
        matchNo = matchNo + 1;
        continue
    else
        matchYes = matchYes + 1;
    end
    if size(matches,2) > 1
        ['-# Caution: More than one apparent match found #-']
    end

    if actuallyCopy && ~isempty(matches)
        [status,msg] = copyfile([matchFiles(matches(1)).folder,filesep,matchFiles(matches(1)).name],...
            [sinkFolder,filesep,matchFiles(matches(1)).name]);
    end


end

disp([char(10),'Of ',num2str(size(sourceFiles,1)),' source files, ', num2str(matchYes),' matches were found (and ',num2str(matchNo),' non-finds)'])
if actuallyCopy
    disp(['(Files were [apparently] copied)'])
end