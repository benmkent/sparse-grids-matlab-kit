%%
clc
addpath(genpath('..'))
%% Setup a sparse grid interpolant
% Construct a SG approximation for a test function.
n=4;
k=6;
domain = [0*ones(1,n);1*ones(1,n)];
S = create_sparse_grid(n,n+k,@(n) knots_CC(n,0,1,'prob'),@(l)lev2knots_doubling(l), @(ii) sum(ii));
Sr = reduce_sparse_grid(S);

f =  genz(n, 10.0, 0.5, "quadraticdecay", "gaussianpeak");

f_on_Sr = evaluate_on_sparse_grid(f,Sr);

figure(1);
plot_sparse_grids_interpolant(S,Sr,domain,f_on_Sr);

%% Convert to modal structure
Mr=modal_sparse_grid(S,Sr,f_on_Sr,domain,'chebyshev');
% This gives a modal structure for the approximation.
display(Mr)

% This can be evaluated and is equal to the classic structure.
nmc = 1e4;
samples = [rand(nmc,n).*diff(domain,1) + domain(1,:)].';
sg_eval = interpolate_on_sparse_grid(S,Sr,f_on_Sr,samples);

mr_eval = evaluate_modal(Mr, samples);

e_sg = norm(f(samples) - sg_eval)/nmc;
e_mr = norm(f(samples) - mr_eval)/nmc;
fprintf('==================================================================\n')
fprintf('Approximation error for function f\n')
fprintf('------------------------------------------------------------------\n')
fprintf('Mean Square Error in sparse grid approximations:   %E\n', e_sg )
fprintf('Mean Square Error in spectral approximations:      %E\n\n', e_mr )

%% Two spectral approximations can be summed together
k2 = 5; % not equal to k
S2 = create_sparse_grid(n,n+k2,@(n) knots_CC(n,0,1,'prob'),@(l)lev2knots_doubling(l), @(ii) sum(ii));
Sr2 = reduce_sparse_grid(S2);

g = genz(n, 1.0, 0.5, "quadraticdecay", "cornerpeak");
g_on_Sr2 = evaluate_on_sparse_grid(g,Sr2);
Mr2=modal_sparse_grid(S2,Sr2,g_on_Sr2,domain,'chebyshev');
sg2_eval = interpolate_on_sparse_grid(S2,Sr2,g_on_Sr2,samples);
mr2_eval = evaluate_modal(Mr2, samples);

% Check approximation error
norm(g(samples) - sg2_eval)/nmc;
e_mr2 = norm(g(samples) - mr2_eval)/nmc;

% Add two approximations
Mr_sum = add_modal_sparse_grids(Mr,Mr2);
mr_sum_eval = evaluate_modal(Mr_sum,samples);
e_sum = norm(f(samples)+g(samples)- mr_sum_eval)/nmc;

fprintf('\n==================================================================\n')
fprintf('Approximation error for sum of approximations for f and g \n')
fprintf('------------------------------------------------------------------\n')
fprintf('                           MSE \t\tn dof\n')
fprintf('f spectral:                %1.2E\t%d\n', e_mr, Mr.n)
fprintf('g spectral:                %1.2E\t%d\n', e_mr2, Mr2.n)
fprintf('f+g spectral:              %1.2E\t%d\n', e_sum, Mr_sum.n)

%% Truncation
tol = 1e-4;
Mr_truncated = truncate_modal(Mr, tol);
Mr2_truncated = truncate_modal(Mr2, tol);
Mr_sum_truncated = truncate_modal(Mr_sum, tol);

% Check approximation error
mr_t_eval = evaluate_modal(Mr_truncated, samples);
mr2_t_eval = evaluate_modal(Mr2_truncated, samples);
mr_t_sum_eval = evaluate_modal(Mr_sum_truncated, samples);

e_mr_t = norm(f(samples) - mr_t_eval)/nmc;
e_mr2_t = norm(g(samples) - mr2_t_eval)/nmc;
e_sum_t = norm(f(samples)+g(samples)- mr_t_sum_eval)/nmc;

fprintf('\n==================================================================\n')
fprintf('Approximation error for truncated approximations\n')
fprintf('------------------------------------------------------------------\n')
fprintf('Spectral term tolerance:           %E\n',tol)
fprintf('------------------------------------------------------------------\n')
fprintf('                           MSE \t\tn dof\n')
fprintf('f spectral:                %1.2E\t%d\n', e_mr, Mr.n)
fprintf('f spectral (truncated):    %1.2E\t%d\n', e_mr_t, Mr_truncated.n)
fprintf('g spectral:                %1.2E\t%d\n', e_mr2, Mr2.n)
fprintf('g spectral (truncated):    %1.2E\t%d\n', e_mr2_t, Mr2_truncated.n)
fprintf('f+g:                       %1.2E\t%d\n', e_sum, Mr_sum.n)
fprintf('f+g (truncated):           %1.2E\t%d\n', e_sum_t, Mr_sum_truncated.n)

%% Reconstruction
% To turn this back into a sparse grid approximation we identify the
% necessary MI to include the polynomials.
MI_truncated = convert_polydegrees_to_miset(Mr_sum_truncated.poly_degrees, @(l)lev2knots_doubling(l));
St = create_sparse_grid_multiidx_set(MI_truncated, @(n) knots_CC(n,0,1,'prob'),@(l)lev2knots_doubling(l));
St_r = reduce_sparse_grid(St);
fg_on_St_r = evaluate_modal(Mr_sum_truncated,St_r.knots);

St_eval = interpolate_on_sparse_grid(St,St_r,fg_on_St_r,samples);
e_sg_t = norm(f(samples)+g(samples)-St_eval)/nmc;

fprintf('\n==================================================================\n')
fprintf('Approximation error for function f+g (reconstructed) \n')
fprintf('------------------------------------------------------------------\n')
fprintf('Spectral term tolerance:           %E\n',tol)
fprintf('------------------------------------------------------------------\n')
fprintf('                                   MSE\t\tn dof\n')
fprintf('f+g (reconstructed sparse grid):   %1.2E\t%d\n',e_sg_t,St_r.size)
fprintf('f+g (truncated spectral):          %1.2E\t%d\n',e_sum_t,Mr_sum_truncated.n)