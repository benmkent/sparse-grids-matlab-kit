function [inds,env] = coeff_envelope(i,x)
inds = min(i):max(i);
env = zeros(size(inds));
for jj = 1:length(inds)
    env(jj) = max(abs(x(i>=inds(jj))));
end
end
