function [x,pk] = ProxL1L1(g,u,lbd,L)
%% Solves the proximal operator of g(x) = lbd ||x||_1 w.r.t l1-norm
% x \in argmin g^Tx + lbd ||x||_1 + L/2 ||x - u||_1^2
% pk \in subdiff( g + lbd ||x||_1) && subdiff(-L/2 ||x - u||_1^2)
% TO DO: clean up the code
%%
p=length(u);
x_new = u;
s = sign(u);
ii= s==0; 

s(ii) = -(abs(g(ii))>= lbd).*sign(g(ii));
%s(ii) = - sign(SoftTh(g(ii), lbd));

strue = zeros(p,1);
c_used = 0;
w = g + lbd*s;
count = 0;

for k=1:2*p
    count = count + 1;
    [slope,i]=max(abs(s.*w));
    %slope = - abs(s(i)*w(i));
    c_left = slope/L - c_used;
    if (c_left <= 0)
        break
    end
    x = x_new;
    if (sign(s(i)*w(i)) > 0)
        x_new(i) = s(i)*max(abs(x(i))-c_left,0);
        strue(i) = s(i);
        c_used = c_used -abs(x_new(i)) + abs(x(i));
        if x_new(i)==0
            %s(i) = - sign(SoftTh(g(i), lbd));
            s(i) = -(abs(g(i))>= lbd).*sign(g(i));
            strue(i) = s(i);
            w(i) = g(i) + lbd*s(i);
        else 
             break
        end
    else
        x_new(i)= s(i)*(abs(x(i))+c_left);
        strue(i) = s(i);
        c_used = c_used + abs(x_new(i)) - abs(x(i));
        break
    end
end

c = c_used;
x = x_new;
 
pk = s.^2.*(g + lbd*s); 
% if xi~= 0, the above choice is obligatory 
% if xi=0 & ui=0, the above choice is a valid one
% if xi=0 & xi-ui~=0, the choice pki = -L*c*sign(x-u) is obligatory
if 0
% when multiple valid choices are possible choose a random one
I = (x==0) & (x==u);
a = zeros(p,1);
b = a;

a(I) = max(g(I)-lbd,-L*c);
b(I) = min(g(I)+lbd,L*c);

pk(I) = a(I) + (b(I) -a(I)).*rand(nnz(I),1);
end
%[s1,s2] = Signs(x-u);

%pk_l = min(-L*c*s1,-L*c*s2);
%pk_u = max(-L*c*s1,-L*c*s2);
%ind = ((x==0)+(pk_l==pk_u)==2);
%pk(ind) = pk_l(ind);
ind = ((x==0) + (x~= u) ==2);
pk(ind)=-L*c*sign(x(ind)-u(ind));

% cond = (pk >= pk_l- 0.0001) + (pk <= pk_u + 0.0001) + (pk <= g+lbd+0.0001) + (pk>=g-lbd-0.001);
% if(sum(cond)~=4*p)
%     display('pk is wrong!')
%     keyboard
% end
end