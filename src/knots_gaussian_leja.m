function [x,w]=knots_gaussian_leja(n)

% [x,w] = knots_gaussian_leja(n)
%
% DEPRECATED!! This function has been replaced by KNOTS_NORMAL_LEJA in release 22.2 and will disappear in future releases
%
%
% returns the collocation points (x) and the weights (w) 
% for the weighted Leja sequence for integration 
% w.r.t to the weight function 
%
% rho(x)=1/sqrt(2*pi)*exp(-x^2/2) 
%
% i.e. the density of a gaussian random variable 
% with mean 0 and standard deviation 1.
%
% Knots and weights have been precomputed (up to 150) following the work
% 
% A. Narayan, J. Jakeman, "Adaptive Leja sparse grid constructions for stochastic collocation and high-dimensional approximation"
% SIAM Journal on Scientific Computing,  Vol. 36, No. 6, pp. A2952–A2983, 2014
%
% an error is raised if more than 150 points are requested.
%
% knots are sorted increasingly before returning (weights are returned in the corresponding order)

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile, B. Sprungk
% See LICENSE.txt for license
%----------------------------------------------------

error('SparseGKit:RenamedFunction',['knots_gaussian_leja has been replaced by knots_normal_leja in release 22.2, '...
    'which takes as input also the mean and the standard deviation of the normal (Gaussian) pdf '  ...
    'to allow for integration w.r.t. general normal densities (and not only w.r.t. the standard Gaussian). '...
    'This message will disappear in future releases of the sparse-grid-matlab-kit.'])

if n>150
   error('SparseGKit:OutOfTable',strcat('this number of points is not available:',num2str(n)))

else
    
    % load points and weights
    % SS = load('GaussianLejaPrecomputedKnotsAndWeights90.mat');
    % x = SS.X(1:n);
    % w = SS.W(1:n,n);

    %[X,W] = get_GaussianLejaPrecomputedKnotsAndWeights90;
    [X,W] = get_GaussianLejaPrecomputedKnotsAndWeights150;
    x = X(1:n);
    w = W(1:n,n);
    
    
    % sort knots increasingly and weights accordingly. Weights need to be row vectors
    [x,sorter]=sort(x);
    w=w(sorter)';

end

end