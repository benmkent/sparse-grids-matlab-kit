function pp_adaptf = post_process_adapt(adaptf, knots, lev2knots)
% POST_PROCESS_ADAPT Postprocesses an adapted sparse grid to remove the contributions from reduced margin.
%
% pp_adaptf = post_process_adapt(adaptf, knots, lev2knots)
% Input: adaptf     adapted sparse grid
%        knots      knot_fn used in sparse grid construction
%        lev2knots  lev2knots function  used in sparse grid construction
% Output: pp_adaptf postprocess adapt_sparse_grid structure
%
% See also ADAPT_SPARSE_GRID
%
%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license

% Use the multi-index set private.I to define the approximation.
actual_MI = adaptf.private.I;
actual_S = create_sparse_grid_multiidx_set(actual_MI, knots, lev2knots);
actual_Sr = reduce_sparse_grid(actual_S);
actual_nb_pts = actual_Sr.size;
f_on_actual_Sr = interpolate_on_sparse_grid(adaptf.S,adaptf.Sr,adaptf.f_on_Sr, actual_Sr.knots);

% Redefine the structure data with reduced margin.
pp_adaptf = adaptf;
pp_adaptf.S = actual_S;
pp_adaptf.Sr = actual_Sr;
pp_adaptf.nb_pts = actual_nb_pts;
pp_adaptf.f_on_Sr = f_on_actual_Sr;
pp_adaptf.inf = quadrature_on_sparse_grid(f_on_actual_Sr, actual_Sr);

pp_adaptf.private.T = adaptf.S;
pp_adaptf.private.Tr = adaptf.Sr;
pp_adaptf.private.f_on_Tr = adaptf.f_on_Sr
end