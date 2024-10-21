function [profits, idx_bin, Ng, Prof_temp, profit_poly_td, polyseries_I, polyseries_coeffs, U_new,poly_already_exists] = compute_polyseries_profits(T,Tr, f_on_Tr,controls, I, N_full, poly_log)
% [profits, idx_bin, Ng, Prof_temp, profit_poly_td, polyseries_I, polyseries_coeffs, U_new,poly_already_exists] = COMPUTE_POLYSERIES_PROFITS(T,Tr, f_on_Tr,controls, I, N_full, poly_log)
%
% [profits, idx_bin, Ng, Prof_temp, profit_poly_td, polyseries_I, polyseries_coeffs, U_new,poly_already_exists] = COMPUTE_POLYSERIES_PROFITS(T,Tr, f_on_Tr,controls, I, N_full, poly_log)
% Outputs:
%   profits:        maximum spectral polynomial coefficient norm introduced by each idx in idx_bin
%   idx_bin:        index bin
%   Ng:             new idxs in idx_bin
%   Prof_temp:      empty
%   profit_poly_td: total polynomial degree for spectral coefficient used to define profit for each idx
%   polyseries_I:   spectral polynomial degrees
%   polyseries_coeffs: spectal polynomial coefficients
%   U_new:          Spectral polynomial substructures for new multi-indices
%   poly_already_exists: index vector to identify previous existing polynomial terms in spectal represntation.
%
% Inputs:
%   T:              full sparse grid (including reduced margin)
%   Tr:             reduced sparse grid T
%   f_on_Tr:        fn evaluations on Tr
%   controls:       controls structure
%   I:              current sparse grid approximation
%   N_full:         full domain dimension
%   poly_log:       log of polynomials added to approximation I


% Computes profit indicators for each multi-index based on the maximum norm
% of the spectral polynomial coefficients that it introduces to the
% approximation.
% This redefines the profit on each iteration for all multindices
%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------
inds = min(i):max(i);

% Compute polynomial coefficients
polytype = controls.polytype; % Can be cell array of different polynomial types.
[polyseries_coeffs,polyseries_I, U] = convert_to_modal(T,Tr,f_on_Tr,controls.domain,polytype);
polyseries_norms = controls.op_vect(polyseries_coeffs.',0*polyseries_coeffs.');

% Map polyomial series multi-indices to sparse grid multi-indices.
map=zeros(size(polyseries_I,1),length(U));
for ii = 1:size(polyseries_I,1)
    poly_mi = polyseries_I(ii,:);
    for jj=1:length(U)
        if ismember(poly_mi, U(jj).multi_indices,'rows')
            map(ii,jj) = 1; % Polynomial ii is in grid jj
        end
    end
end

% Identify new multi-indices in T.idx
mi_set_T = reshape([T.idx].',[N_full,length(T)]).';
[new_mi , new_mi_idx] = setdiff(mi_set_T,I,'rows');
idx_mask = zeros(length(U),1);
idx_mask(new_mi_idx)=1;

% filter the map of the polynomial indices to only consider new multi-indices (Ng?) by deleting columns
map_new_mi_only = map(:,idx_mask==1);
U_new = U(idx_mask==1);

% zero the rows corresponding to polynomials that already exist in the approximation.
poly_already_exists = ismember(polyseries_I,poly_log,'rows');
%poly_already_exists = sum(map(:,idx_mask==0),2) > 0;
map_new_mi_only(poly_already_exists,:) = 0;

% compute ALL profits
% for each SG multi-index not in I
profit_poly_td=[];
Prof_temp = [];
for ii = 1:size(map_new_mi_only,2)
    % identify the polynomial series coefficients that would be added
    idx_valid = map_new_mi_only(:,ii)==1;
    % add define the profit by the maximum series coefficient norm
    [Prof_temp(ii), max_idx] = max(polyseries_norms(idx_valid));
    polyseries_I_valid = polyseries_I(idx_valid,:);
    profit_poly_td(ii) = sum(polyseries_I_valid(max_idx,:),2);
end
profits = Prof_temp; % delete old profit
Prof_temp =[];
idx_bin = new_mi;
Ng = []; % all MI are considered new
%profit_poly_td contains the corresponding polynomial total degrees. Used for plotting.
end