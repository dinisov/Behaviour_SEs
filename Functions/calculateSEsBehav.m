function R = calculateSEsBehav(fixTimes)

    R = struct;
    
    fixTimes(fixTimes == 0) = NaN;
    
    meanFixTimes = nanmean(fixTimes);
    medianFixTimes = nanmedian(fixTimes);
    
    stdFixTimes = std(fixTimes,[],1,'omitnan');
    nFixTimes = sum(~isnan(fixTimes), 1);
    semFixTimes = stdFixTimes ./ sqrt(nFixTimes);
    
    R.nFixTimes = nFixTimes;
    R.meanFixTimes = meanFixTimes;
    R.semFixTimes = semFixTimes;
    R.medianFixTimes = medianFixTimes;

end