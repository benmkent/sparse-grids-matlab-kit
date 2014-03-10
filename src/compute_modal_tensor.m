function U = compute_modal_tensor(S,S_values,domain,flag)

% COMPUTE_MODAL_TENSOR given a tensor grid and the values on it, re-express the interpolant 
% as a modal expansion.
% 
% U=COMPUTE_MODAL_TENSOR(S,S_VALUES,DOMAIN,'legendre') considers the tensor grid S on the
%       hyper-rectangle DOMAIN with associated point evaluations S_VALUES and converts the
%       resulting lagrangian multivariate interpolant to a sum of Legendre polynomials.
%       DOMAIN is a 2xN matrix = [a1, a2, a3, ...; b1, b2, b3, ...] 
%       defining the hyper-rectangluar domain of the sparse grid: (a1,b1) x (a2,b2) x ...  
%       U is a struct with fields U.size (the number of Legendre polynomials needed), 
%       U.multi_indices (one multi-index per Legendre polynomial), U.coeffs
%
% U=COMPUTE_MODAL_TENSOR(S,S_values,domain,'chebyshev') works as the previous call, using 
%       Chebyshev polynomials
%
% U=COMPUTE_MODAL_TENSOR(S,S_values,domain,'hermite') works as the previous call, using 
%       Hermite polynomials. Here DOMAIN is a 2XN matrix = [mu1, mu2, mu3, ...; sigma1, sigma2, sigma3,...]
%       such that the first variable has normal distribution with mean mu1 and std sigma1
%       and so on.


try
    ismember(flag,{'legendre','chebyshev','hermite'})
catch
    error(['Input argument FLAG unrecognized. '...
        ' Please note that COMPUTE_MODAL_TENSOR does not accept INTERVAL_MAP '...
        'input argument any longer. '...
        'Type help convert_to_modal for help. '...
        'This error message will not be shown in future releases of SPARSE-GRID-MATLAB-KIT'])
end


% I will need the knots in each dimension separately. 
% As the number of knots is different in each direction, I use a cell array

nb_dim=size(S.knots,1);
knots_per_dim=cell(1,nb_dim);

for dim=1:nb_dim
    knots_per_dim{dim}=unique(S.knots(dim,:));
end

% The modal expansion in i-th direction uses up to degree k, with k as follows

degrees=zeros(1,nb_dim);
for dim=1:nb_dim
    degrees(dim) = length(knots_per_dim{dim})-1;
end

% the modal polynomials to be used are s.t. the corresponding multi-indices have
% all components less than or equal to the maximum degree

I = multiidx_box_set(degrees,0);

% return multiindex_set as 
U.multi_indices=I;

nb_multiindices=size(I,1);
U.size=nb_multiindices;

% safety check. I have solve a system, so I need the vandermonde matrix to be squared!
rows = S.size; % one equation per point
cols = nb_multiindices; % one unknown per polynomial

if rows~=cols, error('vandermonde matrix will not be square!'), end


% now create the vandermonde matrix with evaluation of each multi-index in every point
% of the grid

V = zeros(rows,cols);

for c=1:cols
    k = I(c,:);
    %vc = lege_eval_multidim(interval_map(S.knots),k,domain(1,:),domain(2,:));
    switch flag
        case 'legendre'
            vc = lege_eval_multidim(S.knots,k,domain(1,:),domain(2,:));
        case 'hermite'
            vc = herm_eval_multidim(S.knots,k,domain(1,:),domain(2,:));
        case 'chebyshev'
            vc = cheb_eval_multidim(S.knots,k,domain(1,:),domain(2,:));            
        otherwise
            error('unknown family of polynomials')
    end    
    V(:,c)=vc';
end

% now solve the system

U.modal_coeffs = V \ S_values;



