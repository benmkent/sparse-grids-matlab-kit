function [sum_modal_sparse_grid] = add_modal_sparse_grids(modal_sparse_grid_1, modal_sparse_grid_2, subtract_flag)
%ADD_MODAL_SPARSE_GRIDS Adds or subtracts two sparse grid approximations
%   This function combines two sparse grid approximations in their spectral
%   formulation by adding or subtracting their polynomial coefficients for
%   matching terms. The result is a combined sparse grid function.
%
%   Inputs:
%       modal_sparse_grid_1 - The first sparse grid approximation in modal form.
%       modal_sparse_grid_2 - The second sparse grid approximation in modal form.
%       subtract_flag       - A boolean flag indicating whether to subtract
%                             the second sparse grid from the first. If true,
%                             subtraction is performed; otherwise, addition.
%
%   Outputs:
%       sum_modal_sparse_grid - The resulting sparse grid approximation in
%                               modal form after addition or subtraction.
%
%   See also: convert_to_modal

if nargin == 2
    subtract_flag = false;
end

% Make coefficients negative if we wish to subtract.
if subtract_flag == true
    modal_sparse_grid_2.coeffs = -modal_sparse_grid_2.coeffs;
end

%%
% Check poly_types are the same
if strcmp(modal_sparse_grid_1.poly_types, modal_sparse_grid_2.poly_types) == false
    error("Modal sparse grids must use same polynomial expansion")
end

% Check domain and poly_types are the same
if norm(modal_sparse_grid_1.domain - modal_sparse_grid_2.domain) > eps()
    error("Modal sparse grids must use same domain")
end

%%
% Now combine coefficients
sum_modal_sparse_grid = modal_sparse_grid_1;

for ii = 1:size(modal_sparse_grid_2.poly_degrees, 1)
    idx = modal_sparse_grid_2.poly_degrees(ii, :);
    % First try to match term
    [isMatch, rowIndex] = ismember(idx, sum_modal_sparse_grid.poly_degrees, 'rows');

    if isMatch == true
        sum_modal_sparse_grid.coeffs(rowIndex, :) = sum_modal_sparse_grid.coeffs(rowIndex, :) + modal_sparse_grid_2.coeffs(ii, :);
    else
        sum_modal_sparse_grid.poly_degrees(end + 1, :) = idx;
        sum_modal_sparse_grid.coeffs(end + 1, :) = modal_sparse_grid_2.coeffs(ii, :);
        sum_modal_sparse_grid.n = sum_modal_sparse_grid.n + 1;
    end

end

end
