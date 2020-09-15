%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile, B. Sprungk
% See LICENSE.txt for license
%----------------------------------------------------


%% Some tests for Gamma-Leja Nodes
% Standard Gamma distribution (beta=1): rho(x)=x^alpha*exp(-x) 

clear
close all

%% Plot the weights of each Gamma-Leja quadrature rule 

alpha = -0.3; 
 
[X,W]=compute_GammaLejaKnotsAndWeights50(alpha);

% plot(X) when X is matrix plots the columns of X, therefore

figure
semilogy(abs(W),'-') % plots weights of each quadrature rule

figure
plot(W','-') % plots trend of weights of each quadrature rule


%% Testing each Gamma-Leja quadrature on computation of moments of a Gamma random variable

clear

alpha = -0.3; 

% Moments of a Gamma distribution rho(x)=x^alpha*exp(-x): 
% Gamma(alpha1+n)/(Gamma(alpha+1)*beta^n), here beta=1
p_max = 50;
mom = zeros(1,1+p_max);
for p = 0:1:p_max
    mom(p+1) = gamma(alpha+1+p)/gamma(alpha+1);
end

% Quadrature error for polynomials
imax = 50;
err = zeros(1+p_max,imax);
errGLagu = zeros(1+p_max,imax);

% Leja knots 
[X,W]=compute_GammaLejaKnotsAndWeights50(alpha);

% for each formula, we test its approximation of increasing moments
for n=1:50
    % Gamma-Leja quadrature rule using n nodes
    % select the first n knots and the corresponding weights from X,W compute above
    
    [x_Lj,w_Lj]=select_weighted_leja(n,X,W); 
        
    % Gauss-Laguerre quadrature of same accuracy
    [x_GLagu,w_GLagu] = knots_gamma(ceil(n/2),alpha,1); % beta = 1
    
    for p=0:p_max
        if p<n+5 % if the degree is "not too much" compute error
            err(1+p,n) = abs(mom(1+p) - dot(x_Lj.^p,w_Lj) );
            errGLagu(1+p,n) = abs(mom(1+p) - dot(x_GLagu.^p,w_GLagu) );
        else % otherwise,  error is just too much,  we  set it to NaN
            err(1+p,n) = NaN;
            errGLagu(1+p,n) = NaN;    
        end
    end
    
    % Plotting errors
    figure
    semilogy(0:p_max, err(:,n),'o',0:p_max, errGLagu(:,n),'+','LineWidth',2)
    grid on
    set(gca,'FontSize',14)
    legend(sprintf('Gamma-Leja (n=%i)',n),sprintf('Gauss-Laguerre (n=%i)',ceil(n/2)))
    set(legend,'Location','NorthWest')
    xlabel('Polynomial degree p')
    ylabel('Absolute quadrature error')
    
    pause
    
end

%% 1d quadrature - convergence test: increase number of points in the quadrature rule

clear

alpha = 1; 

imax=50;

% function to be integrated

% the function resembles a discontinuous (step function) the higher the factor in the argument of the exp is
f = @(x) 1./(1+exp(0.5*x)); 

% oscillatory function: change the period to see different speed of convergence 
% the higher the factor the more oscillatory the function 
% f = @(x) cos(0.5*x); 

% peak function: the factor in front of the norm acts on the steepness of the peak 
% the higher the number the steeper the peak
% f = @(x) 1./(1+0.1*x.^2); 

quad_Lj = zeros(1,imax);
quad_GLagu = zeros(1,imax);
nb_pts =zeros(1,imax);

% Leja knots 
[X,W]=compute_GammaLejaKnotsAndWeights50(alpha);

% refining quad rule
for i=1:imax
    
    n = lev2knots_lin(i);
    nb_pts(i) = n;

    [x_Lj,w_Lj] = select_weighted_leja(n,X,W);   
    [x_GLagu,w_GLagu] = knots_gamma(n,alpha,1);
    
    quad_Lj(i) = dot(f(x_Lj),w_Lj);
    quad_GLagu(i) = dot(f(x_GLagu),w_GLagu);
end

% exact integral
[x_GLagu,w_GLagu] = knots_gamma(100,alpha,1);
exact = dot(f(x_GLagu),w_GLagu);
err_Lj = abs(quad_Lj - exact);
err_GLagu = abs(quad_GLagu - exact);


% Plotting errors
figure
semilogy(nb_pts, err_Lj,'-xr','LineWidth',2,'DisplayName','Exponential-Leja pts')
grid on
hold on
semilogy(nb_pts, err_GLagu,'-ob','LineWidth',2,'DisplayName','Gauss-Laguerre pts')

legend show
set(legend,'Location','SouthWest')


ylim([1e-16 10])


%% 2d quadrature - convergence test

clear

alpha = 3; 

% dimension of space
N=2;

% we use a simple TD rule, up to this level
w_max=15;

% function to be integrated
f = @(x) 1/(1+exp(0.1*sum(x))); 
% f = @(x) cos(0.2*sum(x)); 
% f = @(x) 1/(1+0.1*norm(x)^2); 

% Leja knots 
[X,W]=compute_GammaLejaKnotsAndWeights50(alpha);

knots_Lj = @(n) select_weighted_leja(n,X,W);   
knots_GLagu = @(n) knots_gamma(n,alpha,1);

quad_Lj=zeros(1,w_max);
quad_GLagu=zeros(1,w_max);

nb_pts_Lj =zeros(1,w_max);
nb_pts_GLagu =zeros(1,w_max);

% we introduce auxiliary containers to recycle previous evaluations and speed up the computation
S_Lj_old=[];
Sr_Lj_old=[];
evals_Lj_old=[];

S_GLagu_old=[];
Sr_GLagu_old=[];
evals_GLagu_old=[];

% the convergence loop for Leja and Gauss Laguerre
for w=1:w_max       
    
    disp('Gamma-Leja');
    S_Lj = smolyak_grid(N,w,knots_Lj,@lev2knots_2step, @(i) sum(i-1), S_Lj_old); % using 2step rule to ramp faster to rules with high number of points
    Sr_Lj = reduce_sparse_grid(S_Lj);
    [res,evals]= quadrature_on_sparse_grid(f, S_Lj, Sr_Lj, evals_Lj_old, S_Lj_old, Sr_Lj_old);
    quad_Lj(w) = res;
    evals_Lj_old = evals;
    S_Lj_old=S_Lj;
    Sr_Lj_old = Sr_Lj;
    nb_pts_Lj(w) = Sr_Lj.size;

    disp('Gauss-Laguerre');
    S_GLagu = smolyak_grid(N,w,knots_GLagu,@lev2knots_lin, @(i) sum(i-1), S_GLagu_old);
    Sr_GLagu = reduce_sparse_grid(S_GLagu);    
    [res, evals]  = quadrature_on_sparse_grid(f, S_GLagu, Sr_GLagu, evals_GLagu_old, S_GLagu_old, Sr_GLagu_old);
    quad_GLagu(w) = res;
    evals_GLagu_old = evals;
    S_GLagu_old = S_GLagu;
    Sr_GLagu_old = Sr_GLagu;    
    nb_pts_GLagu(w) = Sr_GLagu.size;
end


% exact integral
disp('computing reference solution');
S_GLagu = smolyak_grid(N,w_max+4,knots_GLagu,@lev2knots_lin, @(i) sum(i-1), S_GLagu_old);
Sr_GLagu = reduce_sparse_grid(S_GLagu);
exact = quadrature_on_sparse_grid(f,S_GLagu,Sr_GLagu,evals_GLagu_old,S_GLagu_old,Sr_GLagu_old);

% errors
err_Lj=abs(quad_Lj - exact);
err_GLagu=abs(quad_GLagu - exact);

figure
loglog(nb_pts_Lj, err_Lj,'-xr','LineWidth',2,'DisplayName','Gamma-Leja pts')
grid on
hold on
loglog(nb_pts_GLagu, err_GLagu,'-ob','LineWidth',2,'DisplayName','Gauss Laguerre pts')

legend show
set(legend,'Location','SouthWest')

%% interpolation error in 1D

clear

alpha = -0.2; 

% Leja knots 
[X,W]=compute_GammaLejaKnotsAndWeights50(alpha);

% function to be interpoled
f1 = @(x) 1./(1+0.1*x.^2); 
f2 = @(x) 1./(2+exp(x)); 
f3 = @(x) cos(3*x+1); 

ff = {f1;f2;f3}; 

figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9])

for j=1:3
    
    f = ff{j}; 

    % sample the function 
    samplesize = 50;
    sampleset = (linspace(0,2,samplesize)); 
    F_samples = f(sampleset);

    % the convergence analysis
    imax=50;

    err_Lj=zeros(1,imax);
    err_GH=zeros(1,imax);

    nb_pts =zeros(1,imax);

    for i=1:imax

        n = lev2knots_lin(i);
        nb_pts(i) = n;
        nnn = 1:n;

        % here we build the lagrange interpolant for Leja and evaluate the error
        [x_Lj,w_Lj] = select_weighted_leja(n,X,W);       
        interp_Lj = zeros(1,samplesize);
        for k=nnn
            interp_Lj =  interp_Lj + f(x_Lj(k))*lagr_eval(x_Lj(k), x_Lj(nnn~=k),sampleset);
        end    
        err_Lj(i) = max(abs(F_samples - interp_Lj)); 


        % repeat for Gauss Hermite
        [x_GLagu,w_GLagu]=knots_gamma(n,alpha,1);
        interp_GLagu = zeros(1,samplesize);
        for k=nnn
            interp_GLagu =  interp_GLagu + f(x_GLagu(k))*lagr_eval(x_GLagu(k), x_GLagu(nnn~=k),sampleset);
        end    
        err_GLagu(i) = max(abs(F_samples - interp_GLagu)); 


    end
    
    subplot(2,3,j)
    plot(sampleset, F_samples,'-k','LineWidth',1.5,'DisplayName','f')
    hold on 
    plot(sampleset, interp_Lj,'*r','DisplayName','interp - Leja pts ')
    plot(sampleset, interp_GLagu,'+b','DisplayName','interp - Gauss pts ')
    legend show
    set(legend,'Location','SouthWest')
    title(strrep(char(ff{j}),'@(x)','f='))

    % Plotting errors
    subplot(2,3,3+j)
    semilogy(nb_pts, err_Lj,'-xr','LineWidth',1.5,'DisplayName','Gamma-Leja pts')
    hold on
    semilogy(nb_pts, err_GLagu,'-ob','LineWidth',1.5,'DisplayName','Gauss-Laguerre pts')
    grid on
    legend show
    set(legend,'Location','SouthWest')

end

