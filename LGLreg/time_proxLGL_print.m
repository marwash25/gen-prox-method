close all
clear all
clc
%%
sizeGrp = 10;
fload = sprintf('Results/TimingProxResults/testproxLGL_sizeGrp%d',sizeGrp)
load(fload)
%%
display('Running time for l2-prox of LGL')
fprintf('%s\n', repmat('-', 1, 79));
fprintf('tol \\ p |')
fprintf('      %d     |      %d     |     %d     |     %d     |     %d     |\n',pvalues)
fprintf('%s\n', repmat('-', 1, 79));

for ind_tol = 1:ntols
    fprintf('  %2.0e |',tols(ind_tol))
    for ind_p = 1:nps    
        ind  = (ind_p-1)*ntols + ind_tol;
        fprintf('   %7.3f   |',info_cyclic(ind).time)
    end
    fprintf('\n')
end
fprintf('\n')
%%
display('Running time for linf-prox of LGL + pk')
fprintf('%s\n', repmat('-', 1, 79));
fprintf('tol \\ p |')
fprintf('      %d     |      %d     |     %d     |     %d     |     %d     |\n',pvalues)
fprintf('%s\n', repmat('-', 1, 79));

for ind_tol = 1:ntols
    fprintf('  %2.0e |',tols(ind_tol))
    for ind_p = 1:nps    
        ind  = (ind_p-1)*ntols + ind_tol;
        fprintf('%5.3f + %5.3f|',info_bslp_linf(ind).time, info_bslp_linf(ind).time_pk)
    end
    fprintf('\n')
end