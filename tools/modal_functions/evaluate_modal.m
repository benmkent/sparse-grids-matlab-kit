function f_values = evaluate_modal(Mr,non_grid_points)
% EVALUATE_MODAL evaluates the modal polynomial approximation (surrogate
% model) on a generic point of the parameter space.
%
% F_VALUES = EVALUATE_MODAL(Mr,NON_GRID_POINTS) evaluates the
%       modal approximation a vector-valued function F: R^N -> R^V based on
%       the modal represntation Mr
%       NON_GRID_POINTS is the set of points where one wants to evaluate the polynomial
%       approximation. It is a matrix, each column is a different point).
%
% This code is based upon the convert to model code for constructing
% Vandermode matrices.
%
%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 B. Kent, L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

% declare a global variable controlling verbosity
global MATLAB_SPARSE_KIT_VERBOSE
if isempty(MATLAB_SPARSE_KIT_VERBOSE)
    MATLAB_SPARSE_KIT_VERBOSE=1;
end

% we have changed the dimensions of non_grid_points to follow
% the convention that points are stored as columns. Again, we throw an
% error if the other convention is detected (this also prevents accidental input errors)
[npoly,udim] = size(Mr.coeffs);
ydim = size(Mr.poly_degree,2);
if size(non_grid_points,1)~=ydim
    error('SparseGKit:WrongInput',[ 'Incompatible sizes. Starting from release 14.4 non_grid_points needs to store points as columns, '...
        'i.e. its dimensions must be N x number_of_queried_evaluations (that is, following the same convention as points'...
        'stored in the ''knots'' sparse grids field). '...
        'Type help INTERPOLATE_ON_SPARSE_GRID for details.'])
end


% the other sizes. Observe that, although we have changed the orientation of the inputs for consistency
% matlab processes faster matrices when working column-wise, so we transpose everything to gain efficiency
nb_pts   = size(non_grid_points,2);
% non_grid_points = non_grid_points';

%V = [p1(x1) p2(x1) p3(x1)... ; p1(x2) p2(x2) ...]
V = zeros(nb_pts,npoly);
I = Mr.poly_degree;
flags = Mr.poly_types;
domain = Mr.domain;
if length(flags)==1 || ischar(flags) % the second condition for when the function is called on one sigle family of polynomials
    for c=1:npoly
        k = I(c,:);
        %vc = lege_eval_multidim(interval_map(S.knots),k,domain(1,:),domain(2,:));
        switch flags
            case 'legendre'
                vc = lege_eval_multidim(non_grid_points,k,domain(1,:),domain(2,:));
            case 'hermite'
                vc = herm_eval_multidim(non_grid_points,k,domain(1,:),domain(2,:));
            case 'chebyshev'
                vc = cheb_eval_multidim(non_grid_points,k,domain(1,:),domain(2,:));
            case 'laguerre'
                vc = lagu_eval_multidim(non_grid_points,k,domain);
            case 'generalized laguerre'
                vc = generalized_lagu_eval_multidim(non_grid_points,k,domain(1,:),domain(2,:));
            case 'jacobi'
                vc = jacobi_prob_eval_multidim(non_grid_points,k,domain(1,:),domain(2,:),domain(3,:),domain(4,:));
            otherwise
                error('SparseGKit:WrongInput','unknown family of polynomials')
        end
        V(:,c)=vc';
    end
else
    for c=1:npoly
        k = I(c,:);
        vc=ones(1,rows);
        for n=1:nb_dim
            switch flags{n}
                case 'legendre'
                    % vc = vc.*lege_eval(S.knots(n,:),k(n),domain(1,n),domain(2,n));
                    vc = vc.*lege_eval(non_grid_points(n,:),k(n),domain{n}(1),domain{n}(2));
                case 'hermite'
                    % vc = vc.*herm_eval(S.knots(n,:),k(n),domain(1,n),domain(2,n));
                    vc = vc.*herm_eval(non_grid_points(n,:),k(n),domain{n}(1),domain{n}(2));
                case 'chebyshev'
                    % vc = vc.*cheb_eval(S.knots(n,:),k(n),domain(1,n),domain(2,n));
                    vc = vc.*cheb_eval(non_grid_points(n,:),k(n),domain{n}(1),domain{n}(2));
                case 'laguerre'
                    vc = vc.*lagu_eval_multidim(non_grid_points(n,:),k(n),domain{n});
                case 'generalized laguerre'
                    vc = vc.*generalized_lagu_eval_multidim(non_grid_points(n,:),k(n),domain{n}(1),domain{n}(2));
                case 'jacobi'
                    vc = jacobi_prob_eval_multidim(non_grid_points(n,:),k(n),domain{n}(1),domain{n}(2),domain{n}(3),domain{n}(4));
                otherwise
                    error('SparseGKit:WrongInput','unknown family of polynomials')
            end
        end
        V(:,c)=vc';

    end
end % of for loop on tensor grid

% Now assemble coefficients with polynomials evaluations
f_values = zeros(udim,nb_pts);
for ii=1:nb_pts
    for jj=1:npoly
        f_values(:,ii) = f_values(:,ii)+Mr.coeffs(jj,:).'*V(ii,jj);
    end
end

