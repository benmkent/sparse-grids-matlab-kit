function [inds,env] = coeff_envelope(i,x)
% [inds,env]=COEFF_ENVELOPE(i,x) 
%  
% Computes an envelope env above the data set i \mapsto x for integer i and
% real numbers x.
%
%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 B. Kent, L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------
inds = min(i):max(i);
env = zeros(size(inds));
for jj = 1:length(inds)
    env(jj) = max(abs(x(i>=inds(jj))));
end
end
