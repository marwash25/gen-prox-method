function [x,pk,info] = ProxL2LGL(g,u,lbd,L,B,tol,maxit,errFcn,save)
%% solve prox of LGLinf w.r.t L2-norm using cyclic projection
% x \in argmin g^Tx + L/2 ||x - u||_2^2 + lbd ||x||_LGL
%       argmin L/2 ||x - (u - g/L)||_2^2 + lbd ||x||_LGL
% ||x||_LGL = min_{w >= 0} 1^Tw s.t. Bw >= |x|
% errFcn: distance to true solution if available, otherwise set to Inf.

% we use moreau's decompostion
% proxLGL(z) = z - Proj{|| z_Gi ||_1 <= lbd, \forall i} 

%B = logical(B);

if save
    % Set the clock.
    time1       = tic;
    % Start clock from here to account for the initialization time
    timestart   = toc(time1); 
else
    info = [];
end

z = (u - g/L);
[p,M] = size(B);
ind = 1:M;
w = z;
w_proj = zeros(p,1);
B = logical(B);

% Identify active groups ||z_Gi||_1 > lbd/L
ActiveGrps = ind(abs(z)'*B > lbd/L);
nActGrps = length(ActiveGrps);


if nActGrps == 0
    x = zeros(p,1);
    if save
     	info.count = 0;
        info.time  = toc(time1) - timestart;
        info.error = errFcn(x);
        info.dist  = 0;
    end
    return
end

for iter = 1:maxit 
    w_prev = w;
    % select group to project on its ell1-ball
    grp = ActiveGrps(mod(iter-1,nActGrps)+1);
    indG = B(:,grp);
    
    % project on ell1-ball restricted to this group
    w_proj(~indG) = w_prev(~indG);
    w_proj(indG) = ProjL2L1(zeros(nnz(indG),1),w_prev(indG),lbd/L,1);%project on L1-ball of radius lbd/L
    
    w = z/(iter+1) + w_proj * iter/(iter+1);
    x = z - w;

    if save
    info.IterTime(iter) = toc(time1) - timestart;
    end
    
    %use the 2nd stopping criterion if true solution not available
    if (errFcn(x) <= tol ) || (errFcn(x) == Inf && norm(w - w_prev,2)<= tol)
        break
    end

    if save
        % Start the clock.
        timestart   = toc(time1);
    end
end
%fprintf('Consecutive Dist = %1.6e\n', norm(w - w_prev,2));

if save
        info.time  = sum(info.IterTime);
        info.count = iter;
        info.error = errFcn(x);
        info.dist = norm(w - w_prev)/norm(w_prev); 
end

pk = zeros(size(x));
end