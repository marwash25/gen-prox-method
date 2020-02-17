function [info_cyclic,info_bslp,info_bslp_l1] = testproxLGL(p,sizeGrp,tol,N,randomGrps)
%% Test prox of Latent Group Lasso
% x \in argmin g^Tx + L/2 ||x - u||_inf^2 + lbd ||x||_LGL
% ||x||_LGL = min_{w \in [0,bd]} 1^Tw s.t. Bw >= |x| (if we use my TU
% relaxation for ||x||_inf <=1)
% Otherwise ||x||_LGL = min_{w \in R_+} 1^Tw s.t. Bw >= |x|

%N = 20; %number of monte Carlo runs

error_cyclic = zeros(N,1);
time_cyclic = zeros(N,1);
dist_cyclic = zeros(N,1);

error_bslp = zeros(N,1);
time_bslp = zeros(N,1);
time_pk_bslp = zeros(N,1);
dist_bslp = zeros(N,1);

error_bslp_l1 = zeros(N,1);
time_bslp_l1 = zeros(N,1);
time_pk_bslp_l1 = zeros(N,1);
dist_bslp_l1  = zeros(N,1);

maxit = 1e7;
save = 1;
    
info_cyclic = [];
info_bslp = [];

for iter = 1:N 
     

    if randomGrps 
        %% Generate random groups

        alpha = 2.5; % ratio between p and expanded dimension alpha p 
        M = ceil(alpha*p/sizeGrp); % number of groups

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
        last_ind=1;
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
    %% generate a signal with small group support
    bd = Inf; 

    %choice of setup from ProxLGL paper
    u = randn(p,1);
    L = 1;
    % pick lbd such that all groups are active, 0.8 min_G ||u_G||_1
    lbd = 0.8*min(abs(u)'*B); 
    g = randn(p,1);

   %% build model for LGL-norm
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

    LGLnorm = @(x) lbd*normLGLinf(x,model,params); 

    if true
    %% Solve the l2-prox using cvx 
    cvx_begin quiet
    cvx_precision best
       variables xcvx2(p) wcvx2(M)
       minimize (g'*xcvx2 + L/2*square_pos(norm(xcvx2-u,2))+lbd*sum(wcvx2))
       subject to
            0<= wcvx2 %<=bd
            %-bd <= xcvx2 <= bd
            B*wcvx2>=abs(xcvx2)
    cvx_end
    
    errFcn_cyclic = @(x) norm(x - xcvx2)/norm(xcvx2);
    %% Solve the l2-prox using cyclic projections (this only considers unbounded LGL)

    [x_cyclic,~,info] = ProxL2LGL(g,u,lbd,L,B,tol,maxit,errFcn_cyclic,save);
    time_cyclic(iter) = info.time;
    error_cyclic(iter) = info.error;
    dist_cyclic(iter) = info.dist;      
    end
    

    %% Solve the linf-prox using cvx 
    cvx_begin quiet
    cvx_precision best
       variables xcvx(p) wcvx(M)
       minimize (g'*xcvx + L/2*square_pos(norm(xcvx-u,inf))+lbd*sum(wcvx))
       subject to
            0<= wcvx %<=bd
            %-bd <= xcvx <= bd
            B*wcvx>=abs(xcvx)
    cvx_end
    
    errFcn_bslp = @(x) norm(x - xcvx)/norm(xcvx);
    %%
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

    clear modelpk
    modelpk.sense = ['<'];
    
    accflag = 1; %compute pk 
    acc = tol;
    method = 1;

    [x_bslp,~,info] = ProxLinfLGL(g,u,L,model,modelpk,params,B,bd,lbd,tol,maxit,accflag,acc,errFcn_bslp,save,method);
    
    time_bslp(iter) = info.time;
    time_pk_bslp(iter) = info.time_pk;
    error_bslp(iter) = info.error;
    dist_bslp(iter) = info.dist;
    
     %% Solve the l1-prox using cvx  
    cvx_begin quiet
    cvx_precision best
       variables xcvx1(p) wcvx1(M)
       minimize (g'*xcvx1 + L/2*square_pos(norm(xcvx1-u,1))+lbd*sum(wcvx1))
       subject to
            0<= wcvx1 %<=bd
            %-bd <= xcvx <= bd
            B*wcvx1>=abs(xcvx1)
    cvx_end

    errFcn_bslp_l1 = @(x) norm(x - xcvx1)/norm(xcvx1);
   
    %% Express the constraint ||x - u||_1 <= t as x - u = [eye(p), - eye(p)] alpha, sum(alpha) = t. our variables are: [x,w,alpha]
    
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
    
    accflag = 1; %compute pk 
    acc = tol;
    method = 1;

    [x_bslp_l1,~,info] = ProxL1LGL(g,u,L,model_l1,modelpk,params,B,bd,lbd,tol,maxit,accflag,acc,errFcn_bslp_l1,save,method);

    
    time_bslp_l1(iter) = info.time;
    time_pk_bslp_l1(iter) = info.time_pk;
    error_bslp_l1(iter) = info.error;
    dist_bslp_l1(iter) = info.dist;

end

info_bslp.time = mean(time_bslp);
info_bslp.time_pk = mean(time_pk_bslp);
info_bslp.error = mean(error_bslp);
info_bslp.dist = mean(dist_bslp);

info_bslp_l1.time = mean(time_bslp_l1);
info_bslp_l1.time_pk = mean(time_pk_bslp_l1);
info_bslp_l1.error = mean(error_bslp_l1);
info_bslp_l1.dist = mean(dist_bslp_l1);

info_cyclic.time = mean(time_cyclic);
info_cyclic.error = mean(error_cyclic);
info_cyclic.dist = mean(dist_cyclic);
end
