function [x,pk,info]=ProxLinfLGL(g,u,L,model,modelpk,params,B,bd,lbd,tol,maxit,accflag,acc,errFcn,save,method)
%% solve prox of LGLinf w.r.t Linf-norm using BS + LP
% x \in argmin g^T(x-u) + L/2 ||x - u||_inf^2 + lbd ||x||_LGL
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
rhs = model.rhs;
eps = 1e-8;
% the max slope coming from decreasing the LGL norm is at most the slope 
% corresponding to the ell_1 norm, i.e., we can use this bound: norm(g+lbd,1)
% At t = 0, i.e., x =u, the objective is 0+||u||_LGL
% At t = ||u||_inf, i.e., x= 0, the objective is g^T(-u) + ||0||_LGL
maxSlope = min(norm(g+lbd,1), bd*max(abs(g)'*B)); 

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
        model.rhs = [rhs;tvec(i)*ones(p,1)+u;tvec(i)*ones(p,1)-u];
        %time = tic;
        result = gurobi(model, params);
        %time = toc(time)

        if strcmp(result.status, 'OPTIMAL') %strcmp(result.status, 'INFEASIBLE')
            fx(i) = result.objval;
            xvec(:,i) = result.x(1:p);
            wvec(:,i) = result.x(p+1:end);
        else
            fx(i) = Inf;
            xvec(:,i) = zeros(p,1);
            wvec(:,i) = zeros(M,1);
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

t = norm(x - u,inf); %this should already be true but just in case
if t > acc
infeasible = true;
signxu = sign(x-u);
signx = sign(x);

if method==0
    nactiveGrps = nnz(w);
    activeGrps = find(w~= 0);
    V = zeros(p,nactiveGrps); % decomposition of x into v's
    Vnorm = zeros(p,nactiveGrps); % V/norm(V,inf)
    res = abs(x);

    for k= 1:nactiveGrps
       i = activeGrps(k);
       V(:,k) = min(abs(res.*B(:,i)),w(i)); 
       normInfV = norm(V(:,k),inf);
       if normInfV~= 0
        Vnorm(:,k) = V(:,k)/normInfV;
       end
       res = res - V(:,k);
    end

    z = -L*t*signxu;
    BTg = B'*(signx.*g);
    VnormTg = Vnorm'*(signx.*g);
elseif method ==1
   %%% solve only for pk on indmax, the rest is zero
   indmax = abs(abs(x-u) - t) <= acc; 
   nindmax = nnz(indmax);
   xumax = (x(indmax)-u(indmax))/(-L*t);
   xg = x'*g;
   normLGLx = lbd*sum(w);
   sg = (signx.*g);
   BTg = B'*sg;
   BTmax = B(indmax,:)'.*repmat(signx(indmax)',[M,1]);
      
   modelpk.lb = -Inf*ones(nindmax,1);
   modelpk.ub = Inf*ones(nindmax,1);
   modelpk.obj = zeros(nindmax,1);
   modelpk.A = sparse([xumax';-xumax'; signxu(indmax)'/(-L*t);x(indmax)';-x(indmax)';BTmax; diag(signxu(indmax)); -diag(signx(indmax))]);
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
   modelpk.A = sparse([xu';-xu'; signxu'/(-L*t);x';-x';BTsign; diag(signxu); -diag(signx)]);
end

while infeasible && acc<=1
    
    if method==0
        indmax = abs(abs(x-u) - t) <= acc; 
        %pk(~indmax) = 0; % this is their final correct value
        nindmax = nnz(indmax);

        BTmax = B(indmax,:)'.*repmat((signx(indmax).*z(indmax))',[M,1]);

        VnormTmax = Vnorm(indmax,:)'.*repmat((signx(indmax).*z(indmax))',[nactiveGrps,1]);

        modelpk.obj = zeros(nindmax,1);
        modelpk.lb = zeros(nindmax,1);
        modelpk.A = sparse([ones(1,nindmax);-ones(1,nindmax);-diag(signx(indmax).*z(indmax));BTmax(setdiff(1:M,activeGrps),:)]);
        modelpk.rhs = [1+acc;-1+acc;-signx(indmax).*g(indmax);acc+lbd+BTg(setdiff(1:M,activeGrps))];
        if nactiveGrps >0
            modelpk.A = sparse([modelpk.A;BTmax(activeGrps,:);-BTmax(activeGrps,:);VnormTmax;-VnormTmax]);
            modelpk.rhs = [modelpk.rhs;acc+lbd+BTg(activeGrps); acc-lbd-BTg(activeGrps); acc+lbd+VnormTg; acc-lbd-VnormTg];
        end
    elseif method == 1     
        modelpk.rhs = [t+acc;-t+acc; 1 + acc;normLGLx+xg+acc;-normLGLx-xg+acc;BTg + lbd + acc; acc* ones(nindmax,1); -sg(indmax) + acc];
    else
        modelpk.rhs = [t+acc;-t+acc; 1 + acc;normLGLx+xg+acc;-normLGLx-xg+acc;BTg + lbd + acc; acc* ones(p,1); -sg + acc];
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
            pk(indmax) = z(indmax).*s;
        elseif method ==1
            pk(indmax) = resultpk.x;
        else
            pk = resultpk.x;
        end

    end

end

end
if save
        time_pk = toc(time_pk);
        info.time_pk  = time_pk;
end
end

   