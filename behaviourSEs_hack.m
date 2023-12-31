%behaviour SEs

close all; clear;

addpath('./Functions');

%% import data

files = dir('./Data');

files = files(~[files.isdir],:);

dataFlies = cell(length(files),1);

for i = 1:length(files)
    thisFlyData = table2array(readtable(fullfile('Data',files(i).name)));
    dataFlies{i} = zeros(size(thisFlyData,1),7);
    dataFlies{i}(:,1:5) = thisFlyData(:,3:7);
    dataFlies{i}(:,7) = thisFlyData(:,10);
end

data = cell2mat(dataFlies);

%there are some zero rows for some reason
data = data(logical(sum(data,2)),:);

data = data(data(:,7) > 0,:);

%% analyse data

% groups (1,32),(2,31),(3,30), etc, as representing the same pattern
% avoids very costly flip() operations later
auxSeq = [1:16 16:-1:1];

% get binary value of sequence (up to pattern)
for trial = 1:size(data,1) 
    data(trial,6) = auxSeq(bin2dec(num2str(data(trial,1:5) > 0)) + 1);
end

dataSeq = cell(1,16); %cell array containing all data points for each sequence
nSeq = zeros(1,16); % n data points for each sequence

SEProfile = zeros(16,1);
SEMSEs = zeros(16,1);

for seq = 1:16
    dataSeq{seq} = data(data(:,6)==seq,7);
    nSeq(seq) = size(dataSeq{seq},1);
    SEProfile(seq) = mean(dataSeq{seq});
    SEMSEs(seq) = std(dataSeq{seq})/sqrt(nSeq(seq));
end

dataSeq = dataSeq(seq_eff_order(5));

create_seq_eff_plot(SEProfile(seq_eff_order(5)),[],'errors',SEMSEs(seq_eff_order(5)));

%% perform ANOVA, show global fixTime distribution, and distributions per sequence

anovaTable = zeros(max(nSeq),16);

for s = 1:16
   anovaTable(1:nSeq(s),s) = dataSeq{s}; 
end

anovaTable(anovaTable == 0) = NaN;

anova1(anovaTable);

% global fixTime distribution
histogram(data(:,7)); title ('Global Fix. Time Dist.');

% distributions for all sequences separately
figure;

load('binomial_x_labels_latex_alt_rep.mat','binomial_x_labels_latex');

%this is just to help turn horizontal sequences into vertical ones
ind_horiz = sub2ind(size(binomial_x_labels_latex{1}),1:4,[1 1 1 5]);

nBins = 5;

for seq = 1:16
   subplot(4,4,seq);
   histogram(dataSeq{seq},nBins);
   title(binomial_x_labels_latex{seq}(ind_horiz));
end