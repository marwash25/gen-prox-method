function [x,pk] = ProxL2L1(g,u,lbd,L)
%% Solves the proximal operator of g(x) = lbd ||x||_1 w.r.t l1-norm
% x \in argmin g^Tx + lbd ||x||_1 + L/2 ||x - u||_2^2
%     = argmin 0.5 ||x - (u - g/L)||_2^2 + lbd/L ||x||_1
% pk \in subdiff( g + lbd ||x||_1) && subdiff(-L/2 ||x - u||_2^2)
%%
z = u - g/L;
x = max(abs(z)-lbd/L,0).*sign(z);
pk = -L*(x - u);

end