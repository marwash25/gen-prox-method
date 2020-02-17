function [x,pk,info]=ProxL1LGL(g,u,L,model,modelpk,params,B,bd,lbd,tol,maxit,accflag,acc,errFcn,save,method)
%% solve prox of LGLinf w.r.t Linf-norm using BS + LP
% x \in argmin g^T(x-u) + L/2 ||x - u||_1^2 + lbd ||x||_LGL
% ||x||_LGL = min_{0 <=w <= bd} 1^Tw s.t. Bw >= |x|
% errFcn: distance to true solution if available, otherwise set to Inf.

if save
    % Set the clock.
    time1       = tic;
    % Start clock from here to account for the initialization time
    timestart   = toc(time1); 
else
    info = [];
end

[p,M] = size(B);
model.obj = [g;model.obj];
rhs = [model.rhs;u;-u];
eps = 1e-8;
% the max slope coming from decreasing the LGL norm is at most the slope 
% corresponding to the ell_1 norm, i.e., we can use this bound: norm(g+lbd,1)
% At t = 0, i.e., x =u, the objective is 0+||u||_LGL
% At t = ||u||_inf, i.e., x= 0, the objective is g^T(-u) + ||0||_LGL
% So the max slope is ...?
maxSlope = min(norm(g+lbd,1), bd*max(abs(g)'*B)); 
%maxSlope = min(norm(g,1), bd*max(abs(g)'*B)); 

t_max = maxSlope/L;

%t_max = 2*norm(xcvx-u,inf);
t_min = 0; %smallest feasible t, 0 if u is feasible
t_feasible = t_min;

fx = zeros(1,3);
xvec = zeros(p,3);
wvec = zeros(M,3);
x = xvec(:,2);
x_prev = ones(p,1);

%We use Binary search to find zero of the fct f(t) = t - slope/L;
for iter = 1:maxit 
   
    t = .5*(t_min + t_max);

    if t-eps<t_feasible
        tvec = [t+eps,t];
    else
        tvec = [t+eps,t,t-eps];
    end   

    for i = 1:3
        %models(i) = model;
        model.rhs = [rhs;tvec(i);-tvec(i)];
        %time = tic;
        result = gurobi(model, params);
        %time = toc(time)

        if strcmp(result.status, 'INFEASIBLE')
            fx(i) = Inf;
            xvec(:,i) = zeros(p,1);
            wvec(:,i) = zeros(M,1);
        else
            fx(i) = result.objval;
            xvec(:,i) = result.x(1:p);
            wvec(:,i) = result.x(p+1:p+M);
            %alpha = result.x(p+M+1:end);
        end
    end

    if t-eps<t_feasible
        slope2 = Inf;
    else
        slope2 = (fx(3) - fx(2))/eps; %absolute value of the slope
    end
    
    slope1 = (fx(2) - fx(1))/eps;
    
    if slope1 <= L*t && L*t <= slope2
        break %reached optimal t
    end
    slope = slope1; %which should also be equal to slope2

    if t - slope/L>0
        t_max = t;
    else
        t_min = t;
    end
    x_prev = x;
    x = xvec(:,2);
    
    if save
    info.IterTime(iter) = toc(time1) - timestart;
    end
    
    %use the 2nd stopping criterion if true solution not available
    if (errFcn(x) <= tol ) || (errFcn(x) == Inf && norm(x - x_prev,2) <= tol) || (t_max == t_min)
        break
    end
    
    if save
        % Start the clock.
        timestart   = toc(time1);
    end
end
% fprintf('Consecutive Dist = %1.6e\n', norm(x - x_prev,2) );
if save
        info.time  = sum(info.IterTime);
        info.count = iter;
        info.error = errFcn(x);
        info.dist = norm(x - x_prev)/norm(x_prev); 
end

pk = zeros(p,1);
time_pk = tic;
if accflag
w = wvec(:,2);
%% compute pk

t = norm(x - u,1); %this should already be true but just in case
if t > acc
infeasible = true;
signxu = sign(x-u);
signx = sign(x);

if method ==1
   %%% solve only for pk on indices where x == u, the rest = -L t sign(x-u)
   pk = -L*t*signxu;
   indxu = abs(x-u) <= acc; %find indices where x==u
   indxu_c = ~indxu;
   nindxu = nnz(indxu);
   %xu_nonzero = (x(indxu)-u(indxu))/(-L*t);
   xpk_c = x(indxu_c)'*pk(indxu_c);
   xg = x'*g;
   normLGLx = lbd*sum(w);
   s_pkcg = signx.*(pk.*indxu_c - g);
   BT_pkcg = B'*s_pkcg;
   BTindxu = B(indxu,:)'.*repmat(signx(indxu)',[M,1]);
      
   modelpk.lb = -Inf*ones(nindxu,1);
   modelpk.ub = Inf*ones(nindxu,1);
   modelpk.obj = zeros(nindxu,1);
   modelpk.A = sparse([ diag(signxu(indxu)/(-L*t));x(indxu)';-x(indxu)';BTindxu; diag(signxu(indxu)); -diag(signx(indxu))]);
   %xu_nonzero';-xu_nonzero';
else  
   %%% solve for full pk
   xu =  (x-u)/(-L*t);
   xg = x'*g;
   normLGLx = lbd*sum(w);
   sg = (signx.*g);
   BTg = B'*sg;
   BTsign = B'.*repmat(signx',[M,1]);
   
   modelpk.lb = -Inf*ones(p,1);
   modelpk.ub = Inf*ones(p,1);
   modelpk.obj = zeros(p,1);
   modelpk.A = sparse([xu';-xu'; diag(signxu/(-L*t));x';-x';BTsign; diag(signxu); -diag(signx)]);
end

while infeasible && acc<=1
    
    if method == 1        
        modelpk.rhs = [ones(nindxu,1) + acc;normLGLx+xg-xpk_c+acc;-normLGLx-xg+xpk_c+acc;-BT_pkcg + lbd + acc; acc* ones(nindxu,1); s_pkcg(indxu) + acc];
    else %full pk
        modelpk.rhs = [t+acc;-t+acc; ones(p,1) + acc;normLGLx+xg+acc;-normLGLx-xg+acc;BTg + lbd + acc; acc* ones(p,1); -sg + acc];
    end
    
    resultpk = gurobi(modelpk, params);
    if strcmp(resultpk.status, 'INFEASIBLE')
        acc = acc*10;
%         display('Infeasible, increase accuracy to')
%         acc
    else 
        infeasible = false;
        if method == 0 
            s = resultpk.x;
            pk(indxu) = z(indxu).*s;
        elseif method ==1
            pk(indxu) = resultpk.x;
        else
            pk = resultpk.x;
        end

    end

end

end
end
if save
        time_pk = toc(time_pk);
        info.time_pk  = time_pk;
end

   