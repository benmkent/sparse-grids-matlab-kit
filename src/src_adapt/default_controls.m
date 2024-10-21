function controls = default_controls(controls,N_full)
% DEFAULT_CONTROLS  Set default controls structure.
%
% controls = DEFAULT_CONTROLS(controls,N_full)
% Set default controls structure. Any controls input within controls are kept otherwise default values are assigned.
% Inputs:
%   controls    input controls stucture with values to be preserved
%   N_full      dimension of function domain
% Outputs:
%   controls    full control structure completed with default values
%
%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

if ~isfield(controls,'pts_tol')
    controls.pts_tol=1e-14;
end
if ~isfield(controls,'max_pts')
    controls.max_pts=1000;
end
if ~isfield(controls,'prof_tol')
    controls.prof_tol=1e-14;
end
if ~isfield(controls,'tol_init')
    controls.tol_init=controls.prof_tol;
end
if ~isfield(controls,'paral')
    controls.paral=NaN;
end
if ~isfield(controls,'prof')
    controls.prof='Linf/new_points';
end
if ~isfield(controls,'plot')
    controls.plot=false;
end
if ~isfield(controls,'nested')
    error('controls must specify the value of ''nested'' field')
end
if ~isfield(controls,'op_vect')
    controls.op_vect = @(A,B) sqrt(sum((A - B).^2,1));
end
if strcmp(controls.prof,'weighted Linf/new_points') || strcmp(controls.prof,'weighted Linf')
    if ~isfield(controls,'pdf')
        error('you need to set the field ''pdf'' to use ''weighted Linf'' and ''weighted Linf/new_points'' profits')
    end
end
if ~isfield(controls,'var_buffer_size')
    controls.var_buffer_size = N_full;
elseif isfield(controls,'var_buffer_size') && controls.var_buffer_size > N_full
    controls.var_buffer_size = N_full;
    warning('SparseGKit:BuffGTNfull','controls.var_buffer_size cannot be greater than N_full. The code will proceed with controls.var_buffer_size = N_full;') 
    pause
end

if ~isfield(controls,'recycling')
    controls.recycling = 'priority_to_evaluation';
end
switch controls.recycling
    case {'priority_to_evaluation','priority_to_recycling'}
    otherwise
        error('unknown value of field controls.recycling')
end

if ~isfield(controls,'polyseries')
    controls.polyseries = false;
end
if ~isfield(controls,'polytype')
   controls.polytype = 'chebyshev'; 
end
if ~isfield(controls,'plateau_gradient')
    controls.plateau_gradient = 0.1;
end
if ~isfield(controls,'safety_factor')
    controls.safety_factor = 2.0;
end
if ~isfield(controls,'burn_in')
    controls.burn_in=0;
end
if ~isfield(controls,'burn_out')
    controls.burn_out=0;
end
end