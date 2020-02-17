function main_l1reg(p,n,s,maxit,noise_sigma)
%% Recovering sparse signal by solving min 0.5 ||y-Ax||_2^2 + lbd ||x||_1 using proximal methods w.r.t to l1 vs l2 norm.

%Generate sparse signal
xtrue = zeros(p,1); 
xS= randsample(p,s);
xtrue(xS)= randn(s,1);

xtrue = xtrue(:); %/norm(xtrue,1);

correlated = 0;
if correlated
        c = 0.5 ; % ~correlation b/w cols
        LL = triu(c + .01*randn(p),1);
        Q = LL + LL' + eye(p);
        min_eig  = min(eig(Q));
        if min_eig <0 %make sure Q is psd
            Q = (Q - (min_eig-1e-7)*eye(p))/(1 - (min_eig-1e-7));
        end
        R = chol(Q);
        A = 1/sqrt(n) * randn(n, p)*R;
else
        A = randn(n,p)/sqrt(n); 
end

w = noise_sigma*randn(n,1); 
b = A*xtrue + w;

lbd = noise_sigma * sqrt(log(p)); 
tol = 1e-9;

gradf = @(x)A'*(A*x-b); 
f = @(x) 0.5*norm(A*x-b,2)^2;
g = @(x) lbd*norm(x,1);

x0 = zeros(p,1);
L_2= normest(A)^2;
L_1 = full(max(max(abs(A'*A))));

cvx_flag = 1;
parameter.verbose = 0;

%parameters for acceleration

eps = 0.03; %log(p) - 1 - sqrt( (log(p) - 1)^2 -1);
d = @(x) 0.5*norm(x- x0,1+eps)^2;
d_sigma = eps/p^(2*eps/(1+eps)); 
%% cvx
xcvx = zeros(p,1);
if cvx_flag
    cvx_begin
    cvx_precision best
       variable xcvx(p)
       minimize 1/2*square_pos(norm(A*xcvx-b,2))+lbd*norm(xcvx,1)
    cvx_end
end

nnz(abs(xcvx)>1e-9)

%%
format long
fprintf('p = %d, n = %d, s = %d, L_2||xcvx||^2_2 = %2.5f, L_1||xcvx||^2_1 = %2.5f, L_1 d(xcvx)/sigma =  %2.5f, ||xcvx - xtrue||_2 = %2.5f \n',p,n,s,L_2*norm(xcvx,2)^2,L_1*norm(xcvx,1)^2, 2*L_1*d(xcvx)/d_sigma, norm(xcvx - xtrue,2))

%% GPM with l1-norm

proxg = @(grad,x_prev,L,iter) ProxL1L1(grad,x_prev,lbd,L);
p_norm = @(x) norm(x,1);
d_norm = @(x) norm(x,inf);

FCmin = @(x,tol) 0;

parameter.x_init = x0;
parameter.x_true = xtrue;
parameter.x_cvx = xcvx;
parameter.maxit = maxit;
parameter.tol = tol;
parameter.Lips = L_1;
parameter.LS = 0;
parameter.FC = 0;
parameter.save = 1;

[x_gpm1,info_gpm1] = GPM(f,gradf,g,proxg,FCmin,p_norm,d_norm,parameter);

parameter.LS = 1;
[x_gpm1_ls,info_gpm1_ls] = GPM(f,gradf,g,proxg,FCmin,p_norm,d_norm,parameter);

%% Accelerated GPM with l1-norm and d(x) = 0.5 ||x - x_init||_{1+eps}^2

proxg = @(grad,x_prev,L,iter) ProxL1L1(grad,x_prev,lbd,L);
p_norm = @(x) norm(x,1);
d_norm = @(x) norm(x,inf);

parameter.Lips = L_1;
parameter.mu = 0;
parameter.sigma = d_sigma;
parameter.heuristic = 0;
parameter.restart = 0;
parameter.tau = 0;

%emin for tau = 0, d(x) = 0.5||x - x0||^2_{1+eps}
emin = @(p,alpha,tau) x0 -sign(p).*abs(p).^(1/eps)./(alpha*norm(p, 1+1/eps)^(-1+1/eps)); 
%emin = @(p,alpha,tau) x0 - p/alpha; %for tau = 0, eps = 1.
parameter.LS = 0;

[x_agpm1_b,info_agpm1_b] = accGPM(f,gradf,g,proxg,d,emin,p_norm,d_norm,parameter);


%% GPM with l2-norm (i.e., ISTA)
proxg = @(grad,x_prev,L,iter) ProxL2L1(grad,x_prev,lbd,L);
p_norm = @(x) norm(x,2);
d_norm = @(x) norm(x,2);
FCmin = @(x,tol) 0;

parameter.x_init = x0;
parameter.x_true = xtrue;
parameter.x_cvx = xcvx;
parameter.maxit = maxit;
parameter.tol = tol;
parameter.Lips = L_2;
parameter.LS = 0;
parameter.FC = 0;

[x_ista,info_ista] = GPM(f,gradf,g,proxg,FCmin,p_norm,d_norm,parameter);

parameter.LS = 1;
[x_ista_ls,info_ista_ls] = GPM(f,gradf,g,proxg,FCmin,p_norm,d_norm,parameter);

%% FISTA 

parameter.LS = 0;
parameter.restart = 0;
[x_fista, info_fista] = FISTA(f,gradf,g,proxg,parameter);

%% Save results
fsave = sprintf('Results/main_l1reg_p%d_n%d_s%d_maxit%d_noisesig%1.0e_synthetic%d',p,n,s,maxit,noise_sigma, 1);
fprintf('saving to %s \n', fsave);
save(fsave);