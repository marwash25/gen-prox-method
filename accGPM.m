function [x,info]=accGPM(f,gradf,g,proxg,d,emin,p_norm,d_norm,parameter)
%% General Proximal Descent Method
% solves min f(x) + g(x) by accelerated proximal method w.r.t to the norm p_norm
% Input:
% f: smooth, L-Lip gradient fct
% g: non-smooth fct with efficient prox w.r.t p_norm
% gradf: gradient of f
% proxg: prox of g w.r.t p_norm, should also return pk, if tau = 0, pk is
% not used, so can return anything
% d: sigma-strongly cvx fct w.r.t p_norm
% emin: estimate sequence minimizer 
% emin(p,alpha,tau) = argmin_w alpha d(w) + w^T p + tau g(w) 
%
% parameter:
%   x_init: initial solution
%   x_true: true solution 
%   maxit: max number of iteration
%   tol: accuracy for stopping criterion
%   Lips: Lipschitz constant
%   LS: do line-search
%   x_cvx: cvx solution 
%   tau: estimate sequence choice, 0 for Bang's sequence, 1 for Nesterov's
%   sigma: strong cvx constant of d w.r.t to p-norm
%   mu: strong cvx constant of g w.r.t to p-norm
%   heuristic: don't use strong cvx constant sigma, instead guess from
%   previous iteration what rho >= ||pk||^2_{1+1/eps}/d_norm(pk)^2 will be, this might
%   improve performance.
%   save: 1 to do bookeeping, 0 othewise 
%   restart: 1 to do restart, 0 otherwise
%   verbose: 1 to print errors, 0 otherwise
% Output:
% x: Solution
% info:
% ObjErr: f(xk) - f*
% EstErr_p: p_norm(xk - x*)
% EstErr_d: d_norm(xk - x*)
% EstErr_2: ||xk - x*||_2
% IterDist_p = p_norm(xk+1 - xk)
%%
if parameter.save
    fprintf('%s\n', repmat('*', 1, 68));
    fprintf('Running Accelerated GPM w.r.t \n')
    p_norm
    if parameter.LS 
        fprintf('with line search \n')
    end
    if parameter.tau
        fprintf('with Nesterov Estimate Sequence  \n')
    else
        fprintf('with Bang Estimate Sequence  \n')
    end
    if parameter.heuristic
        fprintf('using heuristic bound on strong convexity \n')
    else
        fprintf('using worst case bound on strong convexity \n')
    end
    fprintf('%s\n', repmat('*', 1, 68));
    
    info.time = zeros(1,parameter.maxit);
    info.itertime = zeros(1,parameter.maxit);
    info.ObjErr = zeros(1,parameter.maxit);
    info.EstErr2 = zeros(1,parameter.maxit);
    info.EstErr_p = zeros(1,parameter.maxit);
    info.EstErr_d = zeros(1,parameter.maxit);
    info.IterDist_p = zeros(1,parameter.maxit);
    info.OptErr2 = zeros(1,parameter.maxit);
    %info.count = zeros(1,parameter.maxit);
    
    % Set the clock.
    time1       = tic;
else
    info = [];
end

    
    x = parameter.x_init;
    w = x;
    L = parameter.Lips;
    tau = parameter.tau;
    mu = parameter.mu;
    
    beta0 = 1;
    
    if parameter.heuristic
        p = length(x);
        eps = log(p) - 1 - sqrt( (log(p) - 1)^2 -1);
        rho = p^(2*eps/(1+eps));
        parameter.sigma = eps;
    else
        rho = 1;
    end
    Alphas = beta0/parameter.sigma;
    beta = beta0;
    pkacc = 0;
    gradacc = 0;
    Alphas_acc = 0;
    
    
    for iter = 1:parameter.maxit
        
%         if iter ==30
%             keyboard
%         end
        if parameter.save
            % Start the clock.
            timestart   = toc(time1);
        end
        
        x_prev = x;
        bemma = beta/(rho*L);
        alpha = 0.5*(sqrt(bemma^2 + 4*bemma)-bemma);
        beta = (1-alpha)*beta + alpha*tau*mu;
        
        y = (1-alpha)*x + alpha*w;
        
        gradfy = gradf(y);
        %x_prev = x;

        if parameter.LS
            fy = f(y);
            L = L/2;
            while L<= parameter.Lips
                [x,pk] = proxg(gradfy,y,L,iter);
                
                if f(x) <= fy + gradfy'*(x-y) + L/2*p_norm(x-y)^2
                    break;
                end
                L = 2*L;
            end
        else 
            [x,pk] = proxg(gradfy,y,L,iter);
        end
        
            
        if(norm(x-y)<parameter.tol)
            break
        end
        
        if parameter.heuristic
            rho = norm(pk,1+ 1/eps)^2/d_norm(pk)^2;
        end
        
        Alphas = (1- alpha)*Alphas;
        pkacc = alpha * pk + (1-alpha)*pkacc;
        gradacc = alpha * gradfy + (1-alpha)*gradacc;
        Alphas_acc = alpha + (1-alpha)*Alphas_acc;
        
        %w_k+1 \in argmin e_k+1(w) = argmin Alphas d(w) + w^T ((1-tau) pkacc + tau gradacc) + tau*Alphas_acc g(w)
        w = emin((1-tau)*pkacc + tau*gradacc,Alphas,tau*Alphas_acc);
        
        if parameter.save
           info.itertime(iter)  = toc(time1) - timestart;
           info.IterDist_p(iter) = p_norm(x - y);
           info.ObjErr(iter) = f(x) + g(x) - f(parameter.x_cvx) - g(parameter.x_cvx);
           if parameter.verbose % Print information.
             fprintf('Iter = %4d, f(x) = %5.3e\n', iter, info.ObjErr(iter));
           end
           info.EstErr2(iter) = norm(x - parameter.x_true,2);
           info.EstErr_p(iter) = p_norm(x - parameter.x_true);
           info.EstErr_d(iter) = d_norm(x - parameter.x_true);
           info.OptErr2(iter) = norm(x - parameter.x_cvx,2);
           %info.count(iter) = count;
           
           if parameter.restart && iter>1 && info.ObjErr(iter) > info.ObjErr(iter-1)
                Alphas = beta0/parameter.sigma;
                beta = beta0;
                pkacc = 0;
                gradacc = 0;
                Alphas_acc = 0;
                x = x_prev;
                w = x;
               
               info.IterDist_p(iter) = p_norm(x - y);
               info.ObjErr(iter) = f(x) + g(x) - f(parameter.x_cvx) - g(parameter.x_cvx);
               if parameter.verbose % Print information.
                    fprintf('Restarted, Iter = %4d, f(x) = %5.3e\n', iter, info.ObjErr(iter));
               end
               info.EstErr2(iter) = norm(x - parameter.x_true,2);
               info.EstErr_p(iter) = p_norm(x - parameter.x_true);
               info.EstErr_d(iter) = d_norm(x - parameter.x_true);
               info.OptErr2(iter) = norm(x - parameter.x_cvx,2);
           end
        end

    end
    
    if parameter.save
        info.time    = cumsum(info.itertime);
    end
end
