clear all;
close all;

p = 1000;
n = 400;
s = 50;
maxit = 3000;
noise_sigma = 1e-4;
synthetic = 1;
plot_ls = 1;

fload = sprintf('Results/main_l1reg_p%d_n%d_s%d_maxit%d_noisesig%1.0e_synthetic%d',p,n,s,maxit,noise_sigma, synthetic)
load(fload)
%% Plot figures

figure
loglog(info_gpm1.ObjErr,'r-','linewidth',2); hold on
loglog(info_ista.ObjErr,'k-','linewidth',2); hold on
loglog(info_fista.ObjErr,'m-','linewidth',2); hold on
loglog(info_agpm1_b.ObjErr,'b-','linewidth',2); hold on
if plot_ls
hold on, loglog(info_gpm1_ls.ObjErr,'r--','linewidth',2); 
hold on, loglog(info_ista_ls.ObjErr,'k--','linewidth',2);

h = legend('$\ell_1$-GPM','ISTA','FISTA','$\ell_1$-accGPM-$\tau = 0$','$\ell_1$-GPM-ls','ISTA-ls','Location', 'Best');
else
h = legend('$\ell_1$-GPM','ISTA','FISTA','$\ell_1$-accGPM-$\tau = 0$','Location', 'Best');
end
grid on;
set(h,'Interpreter','latex')
title('$F(x_k) - F^*$','Interpreter','Latex')
set(gca,'fontsize',25)
%%
figure
loglog(info_gpm1.EstErr2,'r-','linewidth',2); hold on
loglog(info_ista.EstErr2,'k-','linewidth',2); hold on
loglog(info_fista.EstErr2,'m-','linewidth',2); hold on
loglog(info_agpm1_b.EstErr2,'b-','linewidth',2); hold on
if plot_ls
hold on, loglog(info_gpm1_ls.EstErr2,'r--','linewidth',2); 
hold on, loglog(info_ista_ls.EstErr2,'k--','linewidth',2);
h = legend('$\ell_1$-GPM','ISTA','FISTA','$\ell_1$-accGPM-$\tau = 0$','$\ell_1$-GPM-ls','ISTA-ls','Location', 'Best');
else
h = legend('$\ell_1$-GPM','ISTA','FISTA','$\ell_1$-accGPM-$\tau = 0$','Location', 'Best');
end
grid on;
set(h,'Interpreter','latex')
title('$\|x_k - x^\natural\|_2$','Interpreter','Latex')
set(gca,'fontsize',25)

