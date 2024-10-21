function pt = findchangeptssimple(x)
% FINDCHANGEPTSSIMPLE(x) A simple change point detection implementation
%
% pt = FINDCHANGEPTSSIMPLE(x)
% For a input data series z segments the data using two piecewise linear models if appropriate.
% This is a simple alternative for the signal processing toolkit function findchangepts.
% Input:
%   x:    data series to be modelled by two linear models
% Output:
%   pt:   change pt at which the linear model changes
%
% This is not thoroughly tested.
%
% See also: findchangepts
%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

    penalty = 0.0;
    if isempty(x)
        pt = 1;
        return
    end
    n = length(x);
    % Consider no change pt
    ii=n;
    ind1 = 1:ii;
    partition1 = x(ind1);
    beta1 = [ind1(:) , ones(ii,1)] \ partition1(:);
    res = partition1(:) - [ind1(:) , ones(ii,1)] * beta1;
    totalres(1) = sqrt(sum(res.^2));

    % Consider one change pt
    for ii = 1:length(x)
        %Chooses a point and divides the signal into two sections.
        ind1 = 1:ii;
        ind2 = (ii+1):n;
        partition1 = x(ind1);
        partition2 = x(ind2);

        %Computes an empirical estimate of the desired statistical property for each section.
        % Linear fit
        beta1 = [ind1(:) , ones(ii,1)] \ partition1(:);
        beta2 = [ind2(:) , ones(n-ii,1)] \ partition2(:);

        %At each point within a section, measures how much the property deviates from the empirical estimate. Adds the deviations for all points.
        res1 = partition1(:) - [ind1(:) , ones(ii,1)] * beta1;
        res2 = partition2(:) - [ind2(:) , ones(n-ii,1)] * beta2;

        %Adds the deviations section-to-section to find the total residual error.
        totalres(ii+1) = sqrt(sum(res1.^2) + sum(res2.^2)) + penalty;
    end
    %Varies the location of the division point until the total residual error attains a minimum.
    [minresidual, minind] = min(totalres);
    pt = minind;
end