%behaviour SEs

close all; clear;

addpath('./Functions');

%% import data

files = dir('./Data');

files = files(~[files.isdir],:);

datas = cell(length(files),1);

for i = 1:length(files)
    datas{i} = table2array(readtable(fullfile('Data',files(i).name)));
end

data = cell2mat(datas);

%there are some zero rows for some reason
data = data(logical(sum(data,2)),:);

%% analyse data

% groups (1,32),(2,31),(3,30), etc, as representing the same pattern
% avoids very costly flip() operations later
auxSeq = [1:16 16:-1:1];

% get binary value of sequence (up to pattern)
for trial = 1:size(data,1) 
    % put binary value in 11th column
    data(trial,11) = auxSeq(bin2dec(num2str(data(trial,3:7) > 0)) + 1);
end

SEProfile = zeros(16,1);

SEMSEs = zeros(16,1);

for s = 1:16
    SEProfile(s) = mean(data(data(:,11)==s,10));
    SEMSEs(s) = std(data(data(:,11)==s,10))/sqrt(size(data(data(:,11)==s,10),1));
end

create_seq_eff_plot(SEProfile(seq_eff_order(5)),[],'errors',SEMSEs);