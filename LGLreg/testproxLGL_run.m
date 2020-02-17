function testproxLGL_run(sizeGrp)

pvalues = 2^5; %2.^(5:1:9);  %l2-prox of Linf-LGL is quite slow
tols = [1e-2,1e-3,1e-4,1e-5];
[P,T] = meshgrid(pvalues, tols);

ntols = length(tols);
nps = length(pvalues);
N = 10;
randomGrps = 1;

parfor ind = 1:ntols*nps
    % set seed, to have same examples for each tols and p values
    rng(1)

    %%
    tol = T(ind);
    p = P(ind);
    fprintf('Running testproxLGL with p = %d and tol = %2.5f \n\n',p,tol);
    [info_cyclic(ind),info_bslp_linf(ind),info_bslp_l1(ind)] = testproxLGL(p,sizeGrp,tol,N,randomGrps)
end

fsave = sprintf('Results/testproxLGL_sizeGrp%d',sizeGrp);
save(fsave)
end