function [x,w] = knots_leja(n,x_a,x_b,type,whichrho)

% [X,W] = KNOTS_LEJA(N,A,B,<type>,'prob') returns N Leja points of type <type> between A and B, 
%   for the approximation of
%
%   \int_a^b f(x) * 1/(b-a) dx
%
% [X,W] = KNOTS_LEJA(N,A,B,<type>,'nonprob') returns N Leja points of type <type> between A and B,
%   for the approximation of
%
%   \int_a^b f(x) dx
%
% [X,W] = KNOTS_LEJA(N,A,B,<type>)  is the same as [X,W] = KNOTS_LEJA(N,A,B,<type>,'prob')
%
% Follows the description of the choices of <type>.
%
% -----------------------------------------------------------
%
% [X,W] = KNOTS_LEJA(N,A,B,'line') given X(1)=B, X(2)=A, X(3)=(A+B)/2 recursively
%   defines the n-th point by
%
%   X_n= argmax_[A B] prod_{k=1}^{n-1} abs(X-X_k)
%
%
% [X,W] = KNOTS_LEJA(N,A,B,'sym_line') given X(1)=(A+B)/2, X(2)=B, X(3)=A recursively
%   defines the n-th and (n+1)-th point by
%
%   X_n= argmax_[A B] prod_{k=1}^{n-1} abs(X-X_k)
%   X_(n+1) = symmetric point of X_n with respect to (A+B)/2
%
%
%  [X,W] = KNOTS_LEJA(N,A,B,'p_disk') given  i=sqrt(-1) (imaginary unit), y(1)=1, y(2)=-1, y(3)=i, y(4)=-i,
%   recursively defines the m-th and (m+1)-th point by
%
%   X_m= argmax_[complex unit ball] prod_{k=1}^{m-1} abs(y-y_k)
%
%   and X is obtained by projecting the y sequence on the real interval [-1,1], collecting the first n unique
%   points and translating them to the interval [A,B].


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


if nargin==4
    whichrho='prob';
end

switch type

    %--------------------------------------------------------
    case 'line'
        [x,w] = line_leja_tab(n);


        %--------------------------------------------------------
    case 'sym_line'
        [x,w] = sym_line_leja_tab(n);


        %--------------------------------------------------------
    case 'p_disk'
        [x,w] = p_disk_leja_tab(n);
        %--------------------------------------------------------
    case 'bk'
        [x,w] = leja_bk(n);

        %--------------------------------------------------------
    otherwise
        error('SparseGKit:WrongInput','unknown leja type')

end


% Leja points have been precomputed on the interval (-1,1) and assuming a probabilistic weight, i.e.
% the resulting quadrature rule is a discretization of \int_{-1}^1 f(x) 1/2 dx.
% So now we need to rescale points and weights, if needed


% translation to a generic interval
if x_b~=1 || x_a~=-1
    scale_fact=(x_b-x_a)/2;
    tras=(x_b+x_a)/2;
    x=scale_fact*x+tras;
end

% fix weights
switch whichrho
    case 'prob'
        % as mentioned above, weights precomputed are already probabilistic, so in this case there's nothing to do
    case 'nonprob'
        % in this case, we need to rescale. The 1/2 is already included in the precomputed weights, so we only need to
        % to multiply by (x_b-x_a)
        w=w*(x_b-x_a);
    otherwise
        error('SparseGKit:WrongInput','4th argument must be either prob or nonprob')
end

x=x';
w=w';