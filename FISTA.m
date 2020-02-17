function [x, info] = FISTA(f,gradf,g,proxg,parameter)
%% solves min f(x) + g(x) by proximal method w.r.t to the l2-norm
% Input:
% f: smooth, L-Lip gradient fct
% g: non-smooth fct with efficient prox w.r.t l2_norm
% gradf: gradient of f
% proxg: prox of g w.r.t l2_norm
% parameter:
%   x_init: initial solution
%   x_true: true solution 
%   maxit: max number of iteration
%   tol: accuracy for stopping criterion
%   Lips: Lipschitz constant
%   LS: do line-search
%   restart: do restart
% Output:
% x: Solution
% info:
% ObjErr: f(xk) - f*
% EstErr_2: ||xk - x*||_2
%
%%
if parameter.save
    fprintf('%s\n', repmat('*', 1, 68));
    fprintf('Running FISTA \n')
    fprintf('%s\n', repmat('*', 1, 68));
    
    info.ObjErr = zeros(1,parameter.maxit);
    info.EstErr2 = zeros(1,parameter.maxit);
    info.OptErr2 = zeros(1,parameter.maxit);
    info.time = zeros(1,parameter.maxit);
    info.itertime = zeros(1,parameter.maxit);
    % Set the clock.
    time1       = tic;
else
    info = [];
end
    % Initialize x, y and t.
    x       = parameter.x_init;
    y       = x;
    t       = 1;
   
    L = parameter.Lips;
    
    % Main loop.
    for iter = 1:parameter.maxit
        
        if parameter.save
            % Start the clock.
            timestart   = toc(time1);
        end
        x_prev      = x;      
        t_prev      = t;
        
        gradfy = gradf(y);
         
        % Update next iteration 
         if parameter.LS
            fy = f(y);
            L = L/2;
            while L<= parameter.Lips
                x = proxg(gradfy,y,L,iter);
                
                if f(x) <= fy + gradfy'*(x-y) + L/2*p_norm(x-y)^2;
                    break;
                end
                L = 2*L;
            end
         else 
            
            x = proxg(gradfy,y,L,iter);

         end
        
         if(norm(x-y)<parameter.tol)
            break
         end
        
        t = .5*(1+ sqrt(4*t^2+1));
        gamma  = (t_prev - 1)/t;
        y = x + gamma *(x - x_prev);

       if parameter.save
           info.itertime(iter)  = toc(time1) - timestart;
           info.ObjErr(iter) = f(x) + g(x) - f(parameter.x_cvx) - g(parameter.x_cvx);
           info.EstErr2(iter) = norm(x - parameter.x_true,2);
           info.OptErr2(iter) = norm(x - parameter.x_cvx,2);
       end
       
       if parameter.restart && iter>1 && info.ObjErr(iter) > info.ObjErr(iter-1)

            x = x_prev;
            y = x;
            t = 1;
            if parameter.save
               info.ObjErr(iter) = f(x) + g(x) - f(parameter.x_cvx) - g(parameter.x_cvx);
               info.EstErr2(iter) = norm(x - parameter.x_true,2);
               info.OptErr2(iter) = norm(x - parameter.x_cvx,2);
            end
       end
    end

    if parameter.save
        info.time    = cumsum(info.itertime);
    end
end
