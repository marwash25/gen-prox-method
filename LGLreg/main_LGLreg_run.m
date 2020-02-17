function main_LGLreg_run(p,n,s,sizeGrp,randomGrps,maxit,adapt_tol,tol)

sigma = 1e-3;
tol_GPM = tol;
tol_FISTA = tol;

algos.runGPM = 0;
algos.runaccGPM = 1;
algos.runaccGPM_ls = 0;
algos.runGPM_1 = 0;
algos.runaccGPM_1 = 0;
algos.runaccGPM_ls_1 = 0;
algos.runFISTA = 1;

fprintf('Running main_LGLreg_cluster with p%d_s%d_sizeGrp%d_randomGrps%d_maxit%d \n\n',p,s,sizeGrp,randomGrps,maxit);

[info_fista, info_gpm_inf,info_agpm_inf,info_agpm_inf_ls,info_gpm_1,info_agpm_1,info_agpm_1_ls] = main_LGLreg_parallel(p,n,s,maxit,sigma,algos,sizeGrp,randomGrps,tol_GPM, tol_FISTA,adapt_tol);

tol_str = strrep(num2str(tol), '.', '');
fsave = sprintf('Results/main_LGLreg_p%d_n%d_s%d_sizeGrp%d_randomGrps%d_adaptol%d_tol%s',p,n,s,sizeGrp,randomGrps,adapt_tol,tol_str);
fprintf('saving to %s \n', fsave);
save(fsave);
end