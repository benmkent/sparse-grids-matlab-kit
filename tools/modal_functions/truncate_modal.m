function Mr = truncate_modal(Mr, tol, mass)
% TRUNCATE_MODAL removes modal coefficients with coefficient norm less than
% or equal to tolerance tol.
%
% Mr = truncate_modal(Mr, tol, mass)
%   Mr is the modal struture representing the polynomial approximation
%   tol is the coefficient tolerance
%   mass is the mass matrix used to evaluate coefficient norms.
% See also evaluate_modal, modal_sparse_grid
%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2025 B. Kent, L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

if nargin < 2
    tol = 1e-10;
end
if nargin < 3
    mass = 1;
end

% Check coefficient norms
cnorm = zeros(Mr.n,1);
for ii=1:Mr.n
    c = Mr.coeffs(ii,:);
    cnorm(ii) = sqrt(c*mass*c.');
end

% Find rows to remove
poly_to_delete = cnorm < tol;

% Trim structure
Mr.coeffs(poly_to_delete,:) = [];
Mr.poly_degrees(poly_to_delete,:) = [];
Mr.n = size(Mr.coeffs,1);

