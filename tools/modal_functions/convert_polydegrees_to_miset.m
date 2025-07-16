function [I_mi] = convert_polydegrees_to_miset(I_poly, lev2knots)
%CONVERT_POLYDEGREES_TO_MISET Converts polynomial degrees to multi-indices.
%
%   [I_mi] = CONVERT_POLYDEGREES_TO_MISET(I_poly, lev2knots) takes a matrix
%   of polynomial degrees (I_poly) and a mapping lev2knots that.
%   The function identifies the multi-indices required for each polynomial
%   degree to be included in the sparse grid defined with I_mi.
%
%   Inputs:
%       I_poly    - A matrix where each row represents a set of polynomial
%                   degrees for a specific dimension.
%       lev2knots - A mapping function or array that specifies the
%                   relationship between levels and knots.
%
%   Outputs:
%       I_mi      - Required multi-index set.
%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2025 B. Kent, L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------
I_mi = nan(size(I_poly));

[n_poly,n_d] = size(I_poly);

% Each dimension may have a different lev2knots
if length(lev2knots) == 1
    for ii=1:n_d
        lev2knots_cell{ii} = lev2knots;
    end
elseif length(lev2knots) == n_d
    lev2knots_cell = lev2knots;
else
    error();
end

% For each dimension map levels to poly degree
max_degrees = max(I_poly,[],1);
max_degree = max(max_degrees);

mapmatrix = zeros(max_degree+1,n_d);
for ii=1:n_d
    jj=1;
    degreeprev=1;
    lev2knots_fun = lev2knots_cell{ii};
    while lev2knots_fun(jj-1) <= max_degree+1
        lev2knots_degrees{ii}(jj) = lev2knots_fun(jj);
        degree = min([lev2knots_fun(jj),max_degree+1]);
        mapmatrix(degreeprev:degree,ii) = jj;
        degreeprev = degree+1;
        jj=jj+1;
    end
end

% For each polynomial map to a level
for ii=1:n_poly
    for jj=1:n_d
        I_mi(ii,jj) = mapmatrix(I_poly(ii,jj)+1,jj);
    end
end
I_mi = unique(I_mi,'rows');
[~,I_mi] = check_set_admissibility(I_mi);
end

