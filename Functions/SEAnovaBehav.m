function [P,ANOVATAB,STATS] = SEAnovaBehav(fixTimes)
%SEAnova Performs an ANOVA on SE data

    if any(isnan(fixTimes(:)))
        [P,ANOVATAB,STATS] = anova1(fixTimes);
    else
        fixTimes(~fixTimes) = NaN;
        [P,ANOVATAB,STATS] = anova1(fixTimes);
    end
        
end