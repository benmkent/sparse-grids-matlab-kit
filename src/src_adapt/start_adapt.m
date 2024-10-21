function [N,N_log,var_with_pts,S,Sr,f_on_Sr,I,I_log,poly_log,idx,maxprof,idx_bin,profits,G,G_log,coeff_G,Hr,f_on_Hr,...
    nb_pts,nb_pts_log,num_evals,intf, prof_tol] = start_adapt(f,N,knots,lev2knots,prev_adapt,controls)
% START_ADAPT Initialises or resumes the ADAPT_SPARSE_GRID algorithm
%
% [N,N_log,var_with_pts,S,Sr,f_on_Sr,I,I_log,poly_log,idx,maxprof,idx_bin,profits,G,G_log,coeff_G,Hr,f_on_Hr,nb_pts,nb_pts_log,num_evals,intf, prof_tol] = START_ADAPT(f,N,knots,lev2knots,prev_adapt,controls)
%
% --> I         : is the set of explored indices (the grid is actually larger, it includes as well their neighbours)   
% --> I_log     : I must be sorted lexicographically for software reasons. I_log is the same set of indices, but
%                   sorted in the order in which they are chosen by the algorithm
% --> idx       : is the idx with the highest profit, whose neighbour we will next explore
% --> maxprof   : the corresponding profit
% --> idx_bin   : is the set of idx whose profit has been computed. They have been added to the grid but their neighbour is yet to be explored
% --> profits   : is the corresponding set of profits
% --> G         : is the set of the grid.
% --> G_log     : same as I_log,  but for G
% --> coeff_G   : coefficients of the combination technique applied to G
% --> nb_pts    : the number of points in the grid
% --> num_evals : the number of function evaluations
% --> nb_pts_log: for each iteration, the current nb_pts
% --> S         : create_sparse_grid_multiidx_set(G,knots,lev2knots);
% --> Sr        : reduce_sparse_grid(S);
% --> f_on_Sr   : evaluate_on_sparse_grid(f,Sr);
% --> intf      : approx of integral of f using Sr
% --> Hr        : all the points visited by the algo, stored as a reduced grid to be able to use ; only useful for non-nested grids, where it differs from Sr.
% --> var_with_pts: vector of variables in which we have actually put points. 
%                   length(var_with_pts) + controls.var_buffer_size = N_curr
% --> N_log     : for every iteration, the value of N_curr
%
% See also ADAPT_SPARSE_GRID
%
%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

% declare a global variable controlling verbosity
global MATLAB_SPARSE_KIT_VERBOSE

if isempty(MATLAB_SPARSE_KIT_VERBOSE)
    MATLAB_SPARSE_KIT_VERBOSE=1;
end


if isempty(prev_adapt)

    % it's a fresh start
    %--------------------------------------------

    var_with_pts = []; % we have put no points in no variables for now
    N_log = N;
    I=ones(1,N);
    I_log = ones(1,N);
    poly_log = ones(1,N);
    idx = ones(1,N);
    maxprof = Inf;

    idx_bin=[];
    profits=[];
    
    G = I;    
    G_log = G;
    coeff_G = 1;
    S  = create_sparse_grid_multiidx_set(G,knots,lev2knots);
    Sr = reduce_sparse_grid(S,controls.pts_tol);
    f_on_Sr = evaluate_on_sparse_grid(f,Sr); % here we don't need controls.pts_tol, there is no check on new/old points

    Hr = Sr;
    f_on_Hr = f_on_Sr; % it is a matrix of size VxM where M is the number of points and f:R^N->R^V

    intf=f_on_Sr*Sr.weights';
    nb_pts=size(f_on_Sr,2); 
    nb_pts_log = nb_pts;
    num_evals = nb_pts;

    prof_tol = controls.prof_tol;
else

    % we are resuming from a previous run
    %--------------------------------------------

    if MATLAB_SPARSE_KIT_VERBOSE
        disp('adapt--recycling')
    end

    N = prev_adapt.N;
    N_log = prev_adapt.private.N_log;
    var_with_pts = prev_adapt.private.var_with_pts;
    I=prev_adapt.private.I;
    I_log = prev_adapt.private.I_log;
    poly_log = prev_adapt.private.poly_log;
    idx = prev_adapt.private.idx;
    maxprof = prev_adapt.private.maxprof;

    idx_bin=prev_adapt.private.idx_bin;
    profits=prev_adapt.private.profits;
    
    G = prev_adapt.private.G;    
    G_log = prev_adapt.private.G_log;    
    coeff_G = prev_adapt.private.coeff_G;    
    S = prev_adapt.S;
    Sr= prev_adapt.Sr;
    f_on_Sr = prev_adapt.f_on_Sr;
    
    Hr=prev_adapt.private.Hr;
    f_on_Hr = prev_adapt.private.f_on_Hr;
    
    intf = prev_adapt.intf;
    nb_pts=prev_adapt.nb_pts;  
    nb_pts_log=prev_adapt.private.nb_pts_log;
    num_evals = prev_adapt.num_evals;

    prof_tol = prev_adapt.private.prof_tol;
    
end

end