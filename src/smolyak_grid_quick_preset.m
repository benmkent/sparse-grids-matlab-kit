function [S,Sr,C] = smolyak_grid_quick_preset(N,w)

% [S,Sr,C] = smolyak_grid_quick_preset(N,w)
%
% DEPRECATED!! This function has been replaced by CREATE_SPARSE_GRID_QUICK_PRESET in release 23.5 and will disappear in future releases
%
% SMOLYAK_GRID_QUICK_PRESET is a function meant to be used for quick creation of a preset, "vanilla" sparse grid:
%  Clenshaw--Curtis points in [-1,1] with lev2knots_doubling and multi-index set: sum(ii-1) \leq w 
%  (cf define_functions_for_rule('SM') )
%
%
% [S,Sr,C] = SMOLYAK_GRID_QUICK_PRESET(N,W) returns the Smolyak-type grid in N dimensions at level w. 
%       In practice, it is a short-hand for
%
%       knots = @(n) knots_CC(n,-1,1);
%       [S,C] = create_sparse_grid(N,w,knots,@lev2knots_doubling);
%       Sr = reduce_sparse_grid(S);


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2022 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

error('SparseGKit:RenamedFunction', ['smolyak_grid_quick_preset has been replaced by create_sparse_grid_quick_preset in release 23.5.'...
    ' This message will disappear in future releases of the sparse-grid-matlab-kit.'])

knots = @(n) knots_CC(n,-1,1);

[S,C] = smolyak_grid(N,w,knots,@lev2knots_doubling);

Sr = reduce_sparse_grid(S);