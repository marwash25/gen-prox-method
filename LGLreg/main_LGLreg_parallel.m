function [info_fista, info_gpm_inf,info_agpm_inf,info_agpm_inf_ls, info_gpm_1,info_agpm_1,info_agpm_1_ls] = main_LGLreg_parallel(p,n,s,maxit,sigma,algos,sizeGrp,randomGrps,tol_GPM, tol_FISTA,adapt_tol)
%% Recovering group sparse signal by solving min 0.5 ||y-Ax||_2^2 + lbd ||x||_LGL using proximal methods with different ways to compute the prox.
% where ||x||_LGL = min_{w \in R_+} 1^Tw s.t. Bw >= |x|
%% Generate group sparse signal

runGPM = algos.runGPM;
runaccGPM = algos.runaccGPM;
runaccGPM_ls = algos.runaccGPM_ls;
runGPM_1 = algos.runGPM_1;
runaccGPM_1 = algos.runaccGPM_1;
runaccGPM_ls_1 = algos.runaccGPM_ls_1;
runFISTA = algos.runFISTA;


if randomGrps 
    %% Generate random groups

    alpha = 2.5; % ratio between p and expanded dimension alpha p 
    M = ceil(alpha*p/sizeGrp); % number of groups
    M = p;
    B = zeros(p,M);
    for i = 1:M
       B(randsample(p,sizeGrp),i) = 1;
    end

    %make sure all variables are covered by at least one group
    sumB = sum(B,2);
    ind = sumB == 0;
    if sum(ind) ~=0
        B = [B,ind];
        M = M + 1;
    end

else
    %% Generate interval groups with overlap
   overlap = floor(sizeGrp/3); % size of overlap between grps
   last_ind = 1;
    M=0;
    while last_ind<p
        M=M+1;
        if last_ind+sizeGrp-1 < p
        G{M} = last_ind:last_ind+sizeGrp-1;
        else 
        G{M} = last_ind:p; break
        end
        last_ind = last_ind+sizeGrp-overlap;

    end
    display('Number of groups:')
    M

    B = zeros(p,M); % Construct Group matrix
    for i=1:M
    B(G{i},i)=1;
    end

end

B = sparse(B);

Gsupp = randsample(M,s);
xtrue = min(sum(B(:,Gsupp),2),1); 

bd = norm(xtrue, inf); %upper bound on the latent variables
%bd = Inf;

correlated = 0;
 if correlated
            c = 0.5; % ~correlation b/w cols
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
 
w = sigma*randn(n,1);
b = A*xtrue + w;

tol = 1e-6;
proxmaxit = 1e7;

gradf = @(x)A'*(A*x-b); 
f = @(x) 0.5*norm(A*x-b,2)^2;


nlbds = 5; 
lambdas = logspace(-6,-1,5); %sigma * sqrt(2*log(p));
errors = zeros(nlbds,1);
xcvxs = zeros(p,nlbds);
wcvxs = zeros(M,nlbds);

for i = 1:nlbds
lbd = lambdas(i);
cvx_begin quiet
    cvx_precision best
       variables xcvx(p) wcvx(M)
       minimize (1/2*square_pos(norm(A*xcvx-b,2)) +lbd*sum(wcvx))
       subject to
            0<= wcvx %<=bd
            %-bd <= xcvx <= bd
            B*wcvx>=abs(xcvx)
cvx_end

xcvxs(:,i) = xcvx;
wcvxs(:,i) = wcvx;
errors(i) = norm(xcvx - xtrue,2);
end
[errormin, indmin] = min(errors);
xcvx = xcvxs(:,indmin);
wcvx = wcvxs(:,indmin);
lbd = lambdas(indmin);

format long
fprintf('p = %d, n = %d, s = %d, nnz(xtrue) = %d, lbd = %2.7f, ||xcvx - xtrue|| = %2.5f \n',p,n,s,nnz(xtrue), lbd, norm(xcvx - xtrue,2))

%% build model for g
clear model
model.obj = ones(M,1);

model.A = sparse(B);

model.ub = bd*ones(M,1); 
model.lb = zeros(M,1); 

model.sense = '>';
clear params;
params.Presolve = 2;
params.TimeLimit = 100;
params.outputflag = 0;

g = @(x) lbd*normLGLinf(x,model,params); 
%%
x0 = zeros(p,1);
L_2 = normest(A)^2
L_1 = max(max(abs(A'*A)))
L_inf = min([p*L_2,p^2*L_1,sum(sum(abs(A'*A)))]) %this is just an upper bound, the actual operator norm is NP hard to compute

errFcn = @(x) Inf;
%% Solve using accGPM with l1 prox solved using BS + LP
%if runGPM_1 || runaccGPM_1 || runaccGPM_ls_1
eps = 1;
d_l1 = @(x) 0.5*norm(x- x_init,1+eps)^2;
 
%% build model for proxg

clear model_l1
model_l1.obj = [lbd*ones(M,1);zeros(2*p,1)];

model_l1.A = sparse([[eye(p),-B,zeros(p,2*p)];[-eye(p),-B,zeros(p,2*p)];[eye(p),zeros(p,M),-eye(p),eye(p)];[-eye(p),zeros(p,M),eye(p),-eye(p)];[zeros(1,p+M),ones(1,2*p)];[zeros(1,p+M),-ones(1,2*p)]]);
model_l1.rhs = ([zeros(p,1);zeros(p,1)]);

model_l1.ub = [Inf*ones(p,1);bd*ones(M,1);Inf*ones(2*p,1)]; 
model_l1.lb = [-Inf*ones(p,1);zeros(M+2*p,1)]; 

model_l1.sense = ['<'];
clear params;
params.Presolve = 2;
params.TimeLimit = 100;
params.outputflag = 0;
    

prox_tol = tol_GPM;
prox_maxit = proxmaxit;
method = 1;
acc = prox_tol;
accflag = 1;

clear modelpk
modelpk.sense = ['<'];
save = 0;

proxg_L1LGLacc = @(grad,x_prev,L,iter) ProxL1LGL(grad,x_prev,L,model_l1,modelpk,params,B,bd,lbd,prox_tol/iter^adapt_tol,prox_maxit,accflag,acc,errFcn,save,method); 


FCmin = @(x,tol) 0;
accflag = 0;

proxg_L1LGL = @(grad,x_prev,L,iter) ProxL1LGL(grad,x_prev,L,model_l1,modelpk,params,B,bd,lbd,prox_tol/iter^adapt_tol,prox_maxit,accflag,acc,errFcn,save,method); 


%%
p_norm_l1 = @(x) norm(x,1);
d_norm_l1 = @(x) norm(x,inf);

parameter_L1.x_init = x0;
parameter_L1.x_true = xtrue;
parameter_L1.x_cvx = xcvx;
parameter_L1.maxit = maxit;
parameter_L1.tol = tol;
parameter_L1.Lips = L_1;
parameter_L1.LS = 0;
parameter_L1.save = 1;
parameter_L1.mu = 0;
parameter_L1.sigma = eps/p^(2*eps/(1+eps));
parameter_L1.heuristic = 0;
parameter_L1.restart = 0;
parameter_L1.tau = 0;
parameter_L1.verbose = 0;
parameter_L1.FC = 0;
parameter_L1_ls = parameter_L1;
parameter_L1_ls.LS = 1;

emin_l1 = @(p,alpha,tau) x0 -sign(p).*abs(p).^(1/eps)./(alpha*norm(p, 1+1/eps)^(-1+1/eps)); 

%end
%% Solve using accGPM with linf prox solved using BS + LP
%if runGPM || runaccGPM || runaccGPM_ls
    eps = 1; %use 2-norm
    d = @(x) 0.5*norm(x- x0,1+eps)^2;

    %% build model for proxLinfLGL

    clear model
    model.obj = lbd*ones(M,1);

    model.A = sparse([[eye(p),-B];[-eye(p),-B];[eye(p),zeros(p,M)];[-eye(p),zeros(p,M)]]);
    model.rhs = ([zeros(p,1);zeros(p,1)]);

    model.ub = [Inf*ones(p,1);bd*ones(M,1)]; 
    model.lb = [-Inf*ones(p,1);zeros(M,1)]; 

    model.sense = ['<'];
    clear params;
    params.Presolve = 2;
    params.TimeLimit = 100;
    params.outputflag = 0;

    prox_tol = tol_GPM;
    prox_maxit = proxmaxit;
    method = 1;
    acc = prox_tol;
    accflag = 1;

    clear modelpk
    modelpk.sense = ['<'];
    save = 0;

    proxg_LinfLGLacc = @(grad,x_prev,L,iter) ProxLinfLGL(grad,x_prev,L,model,modelpk,params,B,bd,lbd,prox_tol/iter^adapt_tol,prox_maxit,accflag,acc,errFcn,save,method); 
    
    FCmin = @(x,tol) 0;
    accflag = 0;

    proxLinfLGL = @(grad,x_prev,L,iter) ProxLinfLGL(grad,x_prev,L,model,modelpk,params,B,bd,lbd,prox_tol/iter^adapt_tol,prox_maxit,accflag,acc,errFcn,save,method); 
    
    %% parameters for LGL inf
    p_norm = @(x) norm(x,inf);
    d_norm = @(x) norm(x,1);

    parameter_Linf.x_init = x0;
    parameter_Linf.x_true = xtrue;
    parameter_Linf.x_cvx = xcvx;
    parameter_Linf.maxit = maxit;
    parameter_Linf.tol = tol;
    parameter_Linf.Lips = L_inf;
    parameter_Linf.LS = 0;
    parameter_Linf.save = 1;
    parameter_Linf.mu = 0;
    parameter_Linf.sigma = 1;
    parameter_Linf.heuristic = 0;
    parameter_Linf.restart = 0;
    parameter_Linf.tau = 0;
    parameter_Linf.verbose = 0;
    parameter_Linf.FC = 0;

    parameter_Linf_ls = parameter_Linf;
    parameter_Linf_ls.LS = 1;
                
    emin = @(p,alpha,tau) x0 -sign(p).*abs(p).^(1/eps)./(alpha*norm(p, 1+1/eps)^(-1+1/eps)); 
%end


%% Solve using FISTA with l2 prox solved using cyclic projections
%if runFISTA
prox_tol = tol_FISTA;
prox_maxit = proxmaxit;

proxg_L2LGL = @(grad,x_prev,L,iter) ProxL2LGL(grad,x_prev,lbd,L,B,prox_tol/iter^adapt_tol,prox_maxit,errFcn,save);


parameter_L2.x_init = x0;
parameter_L2.x_true = xtrue;
parameter_L2.x_cvx = xcvx;
parameter_L2.maxit = maxit;
parameter_L2.tol = tol;
parameter_L2.Lips = L_2;
parameter_L2.LS = 0;
parameter_L2.restart = 0;
parameter_L2.save = 1;
parameter_L2.verbose = 1;
%end

%% Run algos in parallel

info_algs = cell(9,1);

parfor algoInd = 1:9
    
    switch algoInd
        case 1
            if runaccGPM 
                [x_agpm_inf,info_algs{algoInd}] = accGPM(f,gradf,g,proxg_LinfLGLacc,d,emin,p_norm,d_norm,parameter_Linf);
            end
        case 2    
            if runaccGPM_ls
                [x_agpm_inf_ls,info_algs{algoInd}] = accGPM(f,gradf,g,proxg_LinfLGLacc,d,emin,p_norm,d_norm,parameter_Linf_ls);
            end
        case 3
            if runGPM
                [x_gpm_inf,info_algs{algoInd}] = GPM(f,gradf,g,proxLinfLGL,FCmin,p_norm,d_norm,parameter_Linf);
            end
        case 4
            if runFISTA
                [x_fista, info_algs{algoInd}] = FISTA(f,gradf,g,proxg_L2LGL,parameter_L2);
            end
        case 5
%             if runCGD
%                 [x_cg, info_algs{algoInd}] = condgradient(f,gradf,b,A,g,LO,parameter_cg);
%             end
        case 6
%             if runCGD_ls
%                 [x_cg_ls, info_algs{algoInd}] = condgradient(f,gradf,b,A,g,LO,parameter_cg_ls);
%             end
        case 7
            if runaccGPM_1
                [x_agpm_1,info_algs{algoInd}] = accGPM(f,gradf,g,proxg_L1LGLacc,d_l1,emin_l1,p_norm_l1,d_norm_l1,parameter_L1);
            end
        case 8  
            if runaccGPM_ls_1
                [x_agpm_1_ls,info_algs{algoInd}] = accGPM(f,gradf,g,proxg_L1LGLacc,d_l1,emin_l1,p_norm_l1,d_norm_l1,parameter_L1_ls);
            end
        case 9
            if runGPM_1
                [x_gpm_1,info_algs{algoInd}] = GPM(f,gradf,g,proxg_L1LGL,FCmin,p_norm_l1,d_norm_l1,parameter_L1);
            end
            
    end
end


info_agpm_inf = info_algs{1};
info_agpm_inf_ls = info_algs{2};
info_gpm_inf = info_algs{3};
info_fista = info_algs{4};
info_agpm_1 = info_algs{7};
info_agpm_1_ls = info_algs{8};
info_gpm_1 = info_algs{9};

end
