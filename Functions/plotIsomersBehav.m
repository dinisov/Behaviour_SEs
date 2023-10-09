function plotIsomersBehav(fixTimes)
%plotIsomers Plot isomers
%   If isomer input is empty it calculates "natural" isomers (all ending with 0/1)
%   "isomer" defines a set of 16 sequences, the other isomer being 1-isomer

        
index1 = zeros(1,2^5); index2 = zeros(1,2^5);

for i = 1:2:31
    index1(i) = i; 
    index2(i+1) = i+1;
end

index1 = index1+fliplr(index1); index1(17:end) = []; index1 = index1(seq_eff_order(5));
index2 = index2+fliplr(index2); index2(17:end) = []; index2 = index2(seq_eff_order(5));

%         disp(sort(union(index1,index2)));

%grab all the ERPs for each isomer, already in SE standard order
fixTimes1 = fixTimes(:,index1); 
fixTimes2 = fixTimes(:,index2);

R1 = calculateSEsBehav(fixTimes1);
R2 = calculateSEsBehav(fixTimes2);

figure; create_seq_eff_plot([R1.meanFixTimes.' R2.meanFixTimes.'],[],'errors',[R1.semFixTimes.' R2.semFixTimes.']);
h = legend({'First Isomer','Second Isomer'});
set(h,'FontSize',6);

h = findobj(gca,'Type','ErrorBar');
set(h(1),'color','r');

end

