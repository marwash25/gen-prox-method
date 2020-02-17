function [x,info]=GPM(f,gradf,g,proxg,FCmin,p_norm,d_norm,parameter)
%% General Proximal Descent Method
% solves min f(x) + g(x) by proximal method w.r.t to the norm p_norm
% Input:
% f: smooth, L-Lip gradient fct
% g: non-smooth fct with efficient prox w.r.t p_norm
% gradf: gradient of f
% proxg: prox of g w.r.t p_norm
% FCmin: minimizer for the fully corrective, only used if FC =1.
% parameter:
%   x_init: initial solution
%   x_true: true solution 
%   maxit: max number of iteration
%   tol: accuracy for stopping criterion
%   Lips: Lipschitz constant
%   LS: do line-search
%   x_cvx: cvx solution 
%   FC: do a fully corrective step
%   save: 1 to do bookeeping, 0 othewise 
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
    fprintf('Running GPM w.r.t \n')
    p_norm
    if parameter.LS 
        fprintf('with line search \n')
    end
    if parameter.FC
        fprintf('with fully corrective step \n')
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
    L = parameter.Lips;

    for iter = 1:parameter.maxit

        if parameter.save
            % Start the clock.
            timestart   = toc(time1);
        end
        gradfx = gradf(x);
        x_prev = x;
        
        if parameter.LS
            
            fx_prev = f(x);
            L = L/2;
            %niter_ls = 0;
            while L<= parameter.Lips
                %niter_ls = niter_ls + 1;
                [x,~] = proxg(gradfx,x_prev,L,iter);
                
                if f(x) <= fx_prev + gradfx'*(x-x_prev) + L/2*p_norm(x-x_prev)^2;
                    break;
                end
                L = 2*L;
            end
            %niter_ls
        else
            [x,~] = proxg(gradfx,x_prev,parameter.Lips,iter);
        end
        
       if(norm(x-x_prev)<parameter.tol)
            break
       end
        
        if parameter.FC
            x = FCmin(x,parameter.tol/iter);
        end

       if parameter.save
           info.itertime(iter)  = toc(time1) - timestart;
           
           info.IterDist_p(iter) = p_norm(x - x_prev);
           info.ObjErr(iter) = f(x) + g(x) - f(parameter.x_cvx) - g(parameter.x_cvx);
           info.EstErr2(iter) = norm(x - parameter.x_true,2);
           
           if parameter.verbose % Print information.
             fprintf('Iter = %4d, f(x) = %5.3e, ||x - x_prev||=%5.3e, ||x - x_true|| = %5.3e \n', iter, info.ObjErr(iter),norm(x-x_prev),  info.EstErr2(iter));
           end
           
           info.EstErr_p(iter) = p_norm(x - parameter.x_true);
           info.EstErr_d(iter) = d_norm(x - parameter.x_true);
           info.OptErr2(iter) = norm(x - parameter.x_cvx,2);
           %info.count(iter) = count;
       end
    
    end
    
    if parameter.save
        info.time    = cumsum(info.itertime);
    end
    
    
end
