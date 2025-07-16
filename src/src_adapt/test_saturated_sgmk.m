function [saturated] = test_saturated_sgmk(idx, saturated_flag, saturated_ipt, lev2knots)
% TEST_SATURATED Tests introduced polynomials for idx are of TD greater than saturation point.
%
%
% Syntax:
%   [saturated] = TEST_SATURATED_SGMK(idx, saturated, saturated_ipt, lev2knots)
%
% Inputs:
%   idx               - Multi-index row matrix.
%   saturated_flag         - True if saturated.
%   saturated_ipt     - Vector of saturation TD thresholds for model.
%   lev2knots         - Level to knots function cell array.
%
% Outputs:
%   saturated         - Vector indicating if idx rows are `saturated'.
%
% Description:
%   For each multi-index in idx checks if the polynomials introduced
%   are saturated by comparing their total degree against the corresponding
%   TD threshold in saturated_ipt. If all polynomials have a total degree
%   greater than the threshold, the index is marked as "saturated".
%
%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2023-2025 L. Tamellini, B.M. Kent
% See LICENSE.txt for license
%----------------------------------------------------

% Initialize the saturated logical array
saturated = false(size(idx, 1), 1);
n = size(idx,2)- d;
Id = eye(n);

if saturated_flag == false
    % All indices are permitted
    return
else
    if ~iscell(lev2knots)
        lev2knots_cell = cell(n,1);
        lev2knots_cell(:) = {lev2knots};
        lev2knots = lev2knots_cell;
    end

    for ii = 1:size(idx, 1)
        % Calculate the lowest total degree of the newly introduced polynomials
        for jj = 1:n
            degrees(jj) =lev2knots{jj}(idx(ii,jj)-1);
        end
        ranges = arrayfun(@(n) 0:n, degrees, 'UniformOutput', false);
        all_poly_degrees = table2array(combinations(ranges{:}));

        % We must now setdiff away existing polynomial terms.
        % I think we must only remove those from reduced margin
        for jj = 1:n
            for kk = 1:n
                degrees(kk) = lev2knots{kk}(idx(ii, kk) - Id(jj,kk))-1;
            end
            ranges = arrayfun(@(n) 0:n, degrees, 'UniformOutput', false);
            all_poly_degrees_jj = table2array(combinations(ranges{:}));
            all_poly_degrees = setdiff(all_poly_degrees,all_poly_degrees_jj,'rows');
        end

        % Check if the minimum new total degree meets the saturation threshold
        if min(sum(all_poly_degrees,2)) >= saturated_ipt(location)
            saturated(ii) = true; % Mark as saturated
        else
            saturated(ii) = false; % Mark as not saturated
        end
    end
end
end
