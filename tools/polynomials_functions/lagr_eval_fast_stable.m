function L = lagr_eval_fast(current_knot,other_knots,ok_len,non_grid_points,ng_size)

% L = LAGR_EVAL_FAST_STABLE(current_knot,other_knots,ok_len,non_grid_points,ng_size)
%
% builds the lagrange function L(x) s.t.
% L(current_knot)=1;
% L(other_knots)=0;
%
% and returns L=L(non_grid_points)
%
% where ok_len = length(other_knots), ng_size = size(non_grid_points);
% 
% this is essentially the same function as LAGR_EVAL_FAST, but improved stability.
% Based upon the barycentric formulation, see e.g.
 % Berrut,  Jean-Paul and Trefethen,  Lloyd N., Barycentric Lagrange Interpolation, 10.1137/s0036144502417715
% Higham,  N. J., The numerical stability of barycentric Lagrange interpolation, 10.1093/imanum/24.4.547
%
%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


% each monodim lagrange function is a product of K terms like (x-x_k)/(current_knot-x_k),
% so i compute it with a for loop

% this is the number of iterations
% ok_len=length(other_knots);

% L is the result. It is a column vector, same size as non_grid_points. It is initialized to 1, and then iteratively multiplied by
% each factor
% L=ones(size(non_grid_points));
L=ones(ng_size);

% for k=1:ok_len
%     % these are the non current-knot, one at a time
%     knots_k=other_knots(k);
%     % here it comes the computation of the lagrangian factor (x-x_k)/(current_knot-x_k)
%     L=L.*(non_grid_points-knots_k)/(current_knot-knots_k);
% end


    ng_minus_grid = non_grid_points(:) - [current_knot, other_knots];
    lx = prod(ng_minus_grid,2);
    wjm1 = prod(current_knot - other_knots);
    L = lx*wjm1^(-1)./(non_grid_points - current_knot);
    L(current_knot == non_grid_points) = 1;
end