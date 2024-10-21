function [leja_points, weights] = knots_leja_generator(n,a,b)
% KNOTS_LEJA_GENERATOR Generates Leja points for n larger than 31 resuing existing SGMK knots
%
% [leja_points, weights] = knots_leja_generator(n,a,b)
% Inputs
%   n       numbers of points
%   a,b     interval [a,b] for the leja knots
% Outputs
%   leja_points     n leja points on interval [a,b]
%   weights         quadrature weights for leja_points
%  
%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------
    % Generate symmetric Leja points
    if n < 31
        [leja_points, weights] = knots_leja(n,a,b,'sym_line');
    else
        leja_points = symmetric_leja_points(n);
        % Compute quadrature weights
        weights = quadrature_weights(leja_points);
        % Rescale
        leja_points = ((0.5*(leja_points+1)))*(b-a)+a;
    end
end

function points = symmetric_leja_points(n)
    % Initialize array for Leja points
    
    % First point is always zero (symmetric about the origin)
    % We also wish to use existing computed knots for compatability
    points(1:31) = knots_leja(31,-1,1,'sym_line');
    
    % Generate remaining points symmetrically
    n_pts = length(points)
    while(n_pts < n)
        % Test points
        testpts = linspace(-1, 1, 2e6);
        
        max_product = -inf;
        max_pt = 0;
        
        for j = 1:length(testpts)
            pt = testpts(j);
            % Calculate product of distances from previously selected points
            product_pt = prod(abs(pt - points));
            
            % Test if product is larger
            if product > max_product
                max_pt = pt;
                max_product = product_pt
            end
        end
        
        % Now add the point symmetrically
        points(end+1:end+2) = [max_pt - max_pt];
        n_pts = length(pts)
    end
    % Trim to required n 
    points = points(1:n);
    % Sort points
    points = sort(points);
end
function weights = quadrature_weights(points)
    n = length(points);
    weights = zeros(1, n);

    [x_cc,w_cc] = knots_CC(ceil(0.5*n)*2+1,-1,1);
    
    for i = 1:n
        quad = sum(prod(x_cc(:) - points([1:i-1, i+1:end]),2).*w_cc(:));
        weights(i) =  quad / prod(points(i) - points([1:i-1, i+1:end]));
    end
end

