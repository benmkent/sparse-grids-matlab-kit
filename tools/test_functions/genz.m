function f = genz(n, C, W, T, N)
% GENZ Generates the Genz test functions based on specified parameters.
%
% INPUTS:
%   n - Number of dimensions (positive integer).
%   C - Scaling coefficient (scalar).
%   W - Shift parameter (scalar).
%   T - Decay type ('nodecay', 'quadraticdecay', 'quarticdecay',
%       'exponentialdecay', 'squaredexponentialdecay').
%   N - Function type ('oscillatory', 'productpeak', 'productpeaknormalised',
%       'cornerpeak', 'gaussianpeak', 'c0continuous', 'discontinuous').
%
% OUTPUT:
%   f - Handle to the generated function.

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

% Minimum coefficient
cmin = 5e-6;

% Generate coefficients based on the decay type
switch T
    case "nodecay"
        c = (1:n + 0.5) / n; % Linear decay
    case "quadraticdecay"
        c = 1 ./ ((1:n) + 1).^2; % Quadratic decay
    case "quarticdecay"
        c = 1 ./ ((1:n) + 1).^4; % Quartic decay
    case "exponentialdecay"
        c = exp(log(cmin) * ((1:n) + 1) ./ n); % Exponential decay
    case "squaredexponentialdecay"
        c = 10.^(log10(cmin) * ((1:n) + 1) / n); % Squared exponential decay
    otherwise
        error("unknown coefficients"); % Error for unknown decay types
end

% Normalize coefficients and apply scaling
c = C * c(:) ./ sum(c);

% Initialize weights
w = ones(size(c)) * W;

% Define the function based on the specified type
switch N
    case "oscillatory"
        f = @(theta) cos(2*pi*w(1) + sum(c .* theta(:))); % Oscillatory function
    case "productpeak"
        f = @(theta) prod(c.^(-2) + (theta - w).^2)^(-1); % Product peak function
    case "productpeaknormalised"
        fpp = genz(n, C, W, T, "productpeak"); % Generate product peak
        f = @(theta) fpp(theta) * prod(c.^(-2)); % Normalized product peak
    case "cornerpeak"
        f = @(theta) (1 + sum(c .* theta)).^(-(n + 1)); % Corner peak function
    case "gaussianpeak"
        f = @(theta) exp(-sum(c.^2 .* (theta - w).^2)); % Gaussian peak function
    case "c0continuous"
        f = @(theta) exp(-sum(c .* abs(theta - w))); % Continuous function
    case "discontinuous"
        f = @(theta) (theta(1) <= w(1)) .* (theta(2) <= w(2)) .* exp(sum(c .* theta)); % Discontinuous function
    otherwise
        error("unknown integrand family"); % Error for unknown function types
end
end