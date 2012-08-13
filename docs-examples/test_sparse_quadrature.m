clear

% function to be interpolated, in (-1,1)^N. input column points, out row vector

f = @(x,b) prod(1./sqrt(x+b));
b=3;
N = 4;
I_1d=(2*sqrt(1+b)-2*sqrt(-1+b));
I_ex = I_1d^N;


% define sparse grid
[lev2knots,idxset]=define_functions_for_rule('TD',N);
knots=@(n) knots_uniform(n,-1,1,'nonprob');


% for loop

w_max=6;
q_error=zeros(1,w_max);
work=zeros(1,w_max);


for w=0:w_max

    disp(w)

    % create grid
    [S,C]=smolyak_grid(N,w,knots,lev2knots,idxset);
    
    Sr=reduce_sparse_grid(S);

    I=quadrature_on_sparse_grid(@(x)f(x,b),S);
    
    %compute convergence error
    q_error(w+1)=abs(I_ex-I);

    % compute work
    work(w+1)=size(Sr.knots,2);
        
end


% error w.r.t. level w
figure
semilogy(0:w_max,q_error,'-or','DisplayName','Numerical Error, w.r.t. grid level');
hold on
semilogy(0:w_max,1./exp(0:w_max),'--o','DisplayName','exp(-level)')
legend show

% error w.r.t. nb. points
figure
loglog(work,q_error,'-or','DisplayName','Numerical Error, w.r.t. #points');
legend show

