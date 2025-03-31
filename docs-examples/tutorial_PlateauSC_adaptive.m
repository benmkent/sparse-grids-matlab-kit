
%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


%% a convergence study in L-inf norm (max in interpolation error)
clear

% function to be interpolated
f=@(x) 1./(x(1,:).^2+x(2,:).^2 + 0.3) + 1e-4*(2*rand(1)-1);

% domain is [a,b]^N
N=2;
a=-1;
b=1;

% settings for sparse grids
knots=@(n) knots_CC(n,a,b);
lev2knots=@lev2knots_doubling;
controls.paral=NaN; 
controls.prof_toll = 1e-10;

controls.nested=true;

% evaluate error as max error over 100 random points in [a,b]^2. Note that here we have hard-coded that a=-1,  b=1
nb_rand_pts = 1e3;
Rand_pts = 2*rand(2,nb_rand_pts)-1;
truef_evals = f(Rand_pts);

%%
% generate a sequence of sparse grids with these many points (approximately), for each save values of interest 
% (exact nb pts, error,  approximation of integral of f)
max_pts = [5 7 13 21 29 50 80 200 400 600 1000];
PP = length(max_pts);
quadf_vals = zeros(1,PP);
sg_pts = zeros(1,PP);
sg_err =zeros(1,PP);

prev_adapt = [];

% the loop over the sparse grids
k=1;
for p = max_pts
    controls.max_pts=p;
    adapt1 = adapt_sparse_grid(f,N,knots,lev2knots,prev_adapt,controls);
    sg_eval = interpolate_on_sparse_grid(adapt1.S,adapt1.Sr,adapt1.f_on_Sr,Rand_pts);
    sg_err(k) = max(abs(sg_eval - truef_evals));
    quadf_vals(k)=adapt1.intf;
    sg_pts(k) = adapt1.nb_pts;
    if p<max_pts(end)
        prev_adapt=adapt1;
        k=k+1;
    end
end

% take the last computed integral as reference integral
quadf_ref = quadf_vals(end);
quadf_vals(end)=[];
max_pts(end)=[];

% error plots
figure(1)
loglog(sg_pts(1:end-1),abs(quadf_vals-quadf_ref),'-ob','LineWidth',2,'MarkerFaceColor','b','DisplayName','adaptive sg')
title('quadrature error')
grid on

figure(2)
loglog(sg_pts,sg_err,'-ob','LineWidth',2,'MarkerFaceColor','b','DisplayName','adaptive sg')
title('interp error')
grid on

%% use plateau detection
% generate a sequence of sparse grids with these many points (approximately), for each save values of interest 
% (exact nb pts, error,  approximation of integral of f)
max_pts = [5 7 13 21 29 50 80 200 400 600 1000];
PP = length(max_pts);
quadf_vals = zeros(1,PP);
sg_pts = zeros(1,PP);
sg_err =zeros(1,PP);

controls.detect_plateau = 1;
controls.domain = [-1,-1;1,1];
controls.min_plateau_length = 3;
controls.burn_in = 3;
controls.burn_out = 3;

prev_adapt = [];

% the loop over the sparse grids
k=1;
for p = max_pts
    controls.max_pts=p;
    adapt1 = adapt_sparse_grid(f,N,knots,lev2knots,prev_adapt,controls);
    sg_eval = interpolate_on_sparse_grid(adapt1.S,adapt1.Sr,adapt1.f_on_Sr,Rand_pts);
    sg_err(k) = max(abs(sg_eval - truef_evals));
    quadf_vals(k)=adapt1.intf;
    sg_pts(k) = adapt1.nb_pts;
    if p<max_pts(end)
        prev_adapt=adapt1;
        k=k+1;
    end
end

% take the last computed integral as reference integral
% quadf_ref = quadf_vals(end);
quadf_vals(end)=[];
max_pts(end)=[];

% error plots
figure(1); hold on;
loglog(sg_pts(1:end-1),abs(quadf_vals-quadf_ref),'--xr','LineWidth',2,'MarkerFaceColor','r','DisplayName','adaptive sg')
hold off;

figure(2); hold on;
loglog(sg_pts,sg_err,'--xr','LineWidth',2,'MarkerFaceColor','r','DisplayName','adaptive sg')
hold off;
