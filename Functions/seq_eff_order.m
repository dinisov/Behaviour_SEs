function [order] = seq_eff_order(n_back)

[~,ind] = sortrows(fliplr(abs(diff(double(dec2bin(0:2^n_back-1))-48,1,2))));

order = ind(ind <= 2^n_back/2);