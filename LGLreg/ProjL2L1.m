function [x,pk] = ProjL2L1(g,u,c,L)
%% Solves the proximal operator of g(x) = indicator{||x||_1 <= c}  w.r.t l2-norm
% x \in argmin g^Tx + L/2 ||x - u||_2^2 s.t  ||x||_1 <= c
%     = argmin 0.5 ||x - (u - g/L)||_2^2 + lbd/L ||x||_1
% pk \in subdiff( g + indicator{||x||_1 <= c}) && subdiff(-L/2 ||x - u||_2^2)
%%
z = u - g/L;

sz   = sort(abs(nonzeros(z)), 'descend');
csz  = cumsum(sz);
nidz = find( csz - (1:numel(sz))'.*[sz(2:end); 0] >= c ...
     + 2*eps(c),1);
if ~isempty(nidz)
    dz   = ( csz(nidz) - c ) /nidz;
    x = z.*( 1 - dz./ max(abs(z), dz) );
else
    x = z;
end
    
pk = -L*(x - u);
end