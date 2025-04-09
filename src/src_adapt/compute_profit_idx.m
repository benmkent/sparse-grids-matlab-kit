function [nb_pts,num_evals,nb_pts_log,Prof_temp,f_on_Tr,Hr,f_on_Hr,intnew] = ...
      compute_profit_idx(ng_idx,f,S,T,Tr,Sr,Hr,f_on_Sr,f_on_Hr,intf,nb_pts,num_evals,nb_pts_log,knots,lev2knots,controls)
% COMPUTE_PROFIT_IDX Computes the profit for the index ng_idx using precomputed grids
%
% [nb_pts,num_evals,nb_pts_log,Prof_temp,f_on_Tr,Hr,f_on_Hr,intnew] = COMPUTE_PROFIT_IDX(ng_idx,f,S,T,Tr,Sr,Hr,f_on_Sr,f_on_Hr,intf,nb_pts,num_evals,nb_pts_log,knots,lev2knots,controls)
%
% This was part of the adapt_sparse_grid loop in old releases.
%
% For a description of inputs and ouputs refer to ADAPT_SPARSE_GRID
%
% See also ADAPT_SPARSE_GRID
%
%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

N = size(ng_idx,2);

if controls.nested

    % here we we evaluate on new points only. Note that finding which points have been evaluated already
    % relies on multiindex info almost exclusively (because the points are nested) so this is quite efficient
    
    [f_on_Tr,new_points,idx_newp] = evaluate_on_sparse_grid(f,T,Tr,f_on_Sr,S,Sr,controls.paral,controls.pts_tol);
    intnew = f_on_Tr*Tr.weights';
    
    newp = length(idx_newp);
    nb_pts=nb_pts + newp;
    nb_pts_log(end+1)=nb_pts; 
    num_evals = nb_pts;
    Hr=[];
    f_on_Hr=[];
    
else
    
    % in this case, we need to keep track of all the points explored, even those that have been discarded in
    % previous iterations
    
    if strcmp(controls.recycling,'priority_to_evaluation') 
        %here we allow for multiple evaluations of the same point because we recycle from the previous grid only.
        % if the function evaluation is "cheap" this is much faster, because the search for common points relies
        % on multiindices and not on comparison of coordinates
        [f_on_Tr,~,idx_newp] = evaluate_on_sparse_grid(f,T,Tr,f_on_Sr,S,Sr,controls.paral,controls.pts_tol);
        intnew = f_on_Tr*Tr.weights';
        
        Hr.knots = [Hr.knots, Tr.knots(:,idx_newp)];
        f_on_Hr = [f_on_Hr, f_on_Tr(:,idx_newp)];

        newp = prod(lev2knots(ng_idx));
        nb_pts = size(Tr.knots,2); 
        num_evals = num_evals + length(idx_newp);
        nb_pts_log(end+1)=nb_pts;
        
    else
        % here we want to make sure no multiple evaluations of the same point occur. Thus we look in 
        % all points ever visited, but this is expensive because we rely on point coordinates only!

        [f_on_Tr,~,idx_newp] = evaluate_on_sparse_grid(f,T,Tr,f_on_Hr,[],Hr.knots,controls.paral,controls.pts_tol); 
        intnew = f_on_Tr*Tr.weights';
        
        Hr.knots = [Hr.knots, Tr.knots(:,idx_newp)];
        f_on_Hr = [f_on_Hr, f_on_Tr(:,idx_newp)];
        
        newp = prod(lev2knots(ng_idx));
        nb_pts = size(f_on_Tr,2); 
        num_evals = size(Hr.knots,2);
        nb_pts_log(end+1)=nb_pts;
    end
    % moreover, if profit is of type Linf, we need to evaluate the new grid on the ``nominally new points'',
    
    switch controls.prof
        
        case {'Linf/new_points','Linf','weighted Linf/new_points','weighted Linf'}
            
            Tx = tensor_grid(N,lev2knots(ng_idx),knots);
            new_points= Tx.knots;
            Tr_on_new_pts = interpolate_on_sparse_grid(T,Tr,f_on_Tr,new_points);
            
            
        case {'deltaint/new_points','deltaint'}
            
            % no need of new points
            
        otherwise
            
            error('do we need new points in this case? fix code here')
    end
    
end % closes if controls.nested


switch controls.prof
    
    case 'Linf/new_points'
        Sr_on_new_pts = interpolate_on_sparse_grid(S,Sr,f_on_Sr,new_points);
        if controls.nested
            Prof_temp = max(controls.op_vect(f_on_Tr(:,idx_newp),Sr_on_new_pts))/newp;
        else
            Prof_temp = max(controls.op_vect(Tr_on_new_pts,Sr_on_new_pts))/newp;
        end
        
    case 'Linf'
        Sr_on_new_pts = interpolate_on_sparse_grid(S,Sr,f_on_Sr,new_points);
        if controls.nested
            Prof_temp = max(controls.op_vect(f_on_Tr(:,idx_newp),Sr_on_new_pts));
        else
            Prof_temp = max(controls.op_vect(Tr_on_new_pts,Sr_on_new_pts));
        end
        
    case 'deltaint/new_points'
        delta_int  = controls.op_vect(intnew,intf);
        Prof_temp = delta_int / newp;
        
    case 'deltaint'
        delta_int  = controls.op_vect(intnew,intf);
        Prof_temp = delta_int;
        
    case 'weighted Linf/new_points'
        Sr_on_new_pts = interpolate_on_sparse_grid(S,Sr,f_on_Sr,new_points);
        if controls.nested
            Prof_temp = max( controls.op_vect(f_on_Tr(:,idx_newp),Sr_on_new_pts).*controls.pdf(new_points) )/newp;
        else
            Prof_temp = max( controls.op_vect(Tr_on_new_pts,Sr_on_new_pts).*controls.pdf(new_points) )/newp;
        end
        
    case 'weighted Linf'
        Sr_on_new_pts = interpolate_on_sparse_grid(S,Sr,f_on_Sr,new_points);
        if controls.nested
            Prof_temp = max( controls.op_vect(f_on_Tr(:,idx_newp),Sr_on_new_pts).*controls.pdf(new_points) );
        else
            Prof_temp = max( controls.op_vect(Tr_on_new_pts,Sr_on_new_pts).*controls.pdf(new_points) );
        end
    otherwise
        error('unknown profit indicator. Check spelling')
end

end