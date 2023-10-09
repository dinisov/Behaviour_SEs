function SEAnovaBehav(fixTimes)
%SEAnova Performs an ANOVA on SE data

    if any(isnan(fixTimes(:)))
        anova1(fixTimes);
    else
        fixTimes(~fixTimes) = NaN;
        anova1(fixTimes);
    end
        
end