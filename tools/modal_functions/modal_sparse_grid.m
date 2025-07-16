function Mr=modal_sparse_grid(S,Sr,nodal_values,domain,flags,~)
% MODAL_SPARSE_GRID wraps the convert_to_modal function to return a modal
% structure.
% This structure is enough to reconstruct the polynomial, and also to
% exploit its spectral representation.
%
% Mr = MODAL_SPARSE_GRID(S,Sr,nodal_values,domain,flags,~)
%       Mr is a structure of modal polynomial degrees and corresponsing spectral coefficients
%       S is the sparse grid
%       Sr is the reduced sparse grid
%       nodal_values is the function evaluations on the reduced sparse grid
%       domain is the function domain
%       flags is the desired spectral polynomials. Either a single string
%       to be used for all parameter dimensions, or a cell array containing
%       a polynomial type for each dimension.
%
% See also convert_to_modal
%
%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2025 B. Kent, L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------
if nargin==4
    error('SparseGKit:WrongInput',strcat('not enough input arguments.',errmsg))
end

if any(~ismember(flags,{'legendre','chebyshev','hermite','laguerre','generalized laguerre','jacobi'}))
    error('SparseGKit:WrongInput',strcat('One or more strings in FLAGS unrecognized. ',errmsg));
end

if iscell(flags) && ~iscell(domain)
    errmsg = ['Input argument DOMAIN must be a cell array. ' ...
        'Please note that CONVERT_TO_MODAL has been changed after release 18.10. ' ...
        'The domain for the case of polynomials of "mixed" type is now a cell array, each cell containing the domain for the corresponding polynomial. ' ...
        'Type help convert_to_modal for help. '...
        'This message will disappear from future relesases of SPARSE-GRID-MATLAB-KIT.'];
    error('SparseGKit:WrongInput',strcat(errmsg));
end

% Convert to modal
[modal_coeffs,K] = convert_to_modal(S,Sr,nodal_values,domain,flags);

% Construct a structure containing spectral data
Mr.coeffs = modal_coeffs;
Mr.poly_degrees = K;
Mr.poly_types = flags;
Mr.domain = domain;
Mr.n = size(modal_coeffs,1);