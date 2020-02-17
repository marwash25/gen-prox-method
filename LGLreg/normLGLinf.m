function [fx,w] = normLGLinf(x,model,params)

%% Compute LGLinf norm of x, i.e., min_{0<=w <=ub} 1^Tw : B w>= |x|
model.rhs = abs(full(x));

result = gurobi(model, params);

if strcmp(result.status, 'INFEASIBLE')
    fx = Inf;
else
fx = result.objval;
end
w = result.x;
end