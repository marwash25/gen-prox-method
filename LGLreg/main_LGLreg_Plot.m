close all
clear all
%% Plots
p = 100;
n = 50;
s = 2;
sizeGrp = 10;
randomGrps = 1;
adapt_tol = 1;
tol = 1e-5;
tol_str = strrep(num2str(tol), '.', '');
fload = sprintf('Results/p100_n50_s2_sizeGrp10_randomGrps1_maxit5000/main_LGLreg_p%d_n%d_s%d_sizeGrp%d_randomGrps%d_adaptol%d_tol%s',p,n,s,sizeGrp,randomGrps,adapt_tol,tol_str)
load(fload)

runCGD = algos.runCGD;
runCGD_ls = algos.runCGD_ls;
runGPM = algos.runGPM;
runaccGPM = algos.runaccGPM;
runaccGPM_ls = algos.runaccGPM_ls;
runGPM_1 = algos.runGPM_1;
runaccGPM_1 = algos.runaccGPM_1;
runaccGPM_ls_1 = algos.runaccGPM_ls_1;
runFISTA = algos.runFISTA;

print_flag = 0;
count = 1;

figure
if runGPM
legendStr{count} = '$\ell_\infty$-GPM';
loglog(info_gpm_inf.ObjErr,'r--','linewidth',2); hold on
count = count + 1;
end

if runaccGPM
legendStr{count} = '$\ell_\infty$-accGPM';
loglog(info_agpm_inf.ObjErr,'r-','linewidth',2); hold on
count = count + 1;
end

if runaccGPM_ls
legendStr{count} = '$\ell_\infty$-accGPM-ls';
loglog(info_agpm_inf_ls.ObjErr,'r-.','linewidth',2); hold on
count = count + 1;
end

if runGPM_1
legendStr{count} = '$\ell_1$-GPM';
loglog(info_gpm_1.ObjErr,'g--','linewidth',2); hold on
count = count + 1;
end

if runaccGPM_1
legendStr{count} = '$\ell_1$-accGPM';
loglog(info_agpm_1.ObjErr,'g-','linewidth',2); hold on
count = count + 1;
end

if runaccGPM_ls_1
legendStr{count} = '$\ell_1$-accGPM-ls';
loglog(info_agpm_1_ls.ObjErr,'g-.','linewidth',2); hold on
count = count + 1;
end

if runFISTA
legendStr{count} = 'FISTA';
loglog(info_fista.ObjErr,'b-','linewidth',2); hold on
count = count + 1;
end

if runCGD
legendStr{count} =   'CGD';
loglog(info_cg.ObjErr,'b-','linewidth',2); hold on
count = count + 1;
end


if runCGD_ls
legendStr{count} =   'CGD';
loglog(info_cg_ls.ObjErr,'b--','linewidth',2); hold on
count = count + 1;
end

h = legend(legendStr,'Location', 'Best');
grid on;
set(h,'Interpreter','latex')
xlabel('iterations')
title('$F(x_k) - F^*$','Interpreter','Latex')
set(gca,'fontsize',25)
shg
%fname = sprintf('ObjErrIter_p%d_s%d_sizeGrp%d_randomGrps%d',p,s,sizeGrp,randomGrps);
fname = sprintf('ObjErrIter_p%d_n%d_s%d_sizeGrp%d_randomGrps%d_adaptol%d_tol%s',p,n,s,sizeGrp,randomGrps,adapt_tol,tol_str);
if print_flag
    print(gcf,'-dpdf','-r600',fname)
end
%%
figure

if runGPM
loglog(info_gpm_inf.time/3,info_gpm_inf.ObjErr,'r--','linewidth',2); hold on
end

if runaccGPM
loglog(info_agpm_inf.time/3,info_agpm_inf.ObjErr,'r-','linewidth',2); hold on
end

if runaccGPM_ls
loglog(info_agpm_inf_ls.time,info_agpm_inf.ObjErr,'r-.','linewidth',2); hold on
end

if runGPM_1
loglog(info_gpm_1.time/3,info_gpm_1.ObjErr,'g--','linewidth',2); hold on
end

if runaccGPM_1
loglog(info_agpm_1.time/3,info_agpm_1.ObjErr,'g-','linewidth',2); hold on
end

if runaccGPM_ls_1
loglog(info_agpm_1_ls.time,info_agpm_1.ObjErr,'g-.','linewidth',2); hold on
end

if runFISTA
loglog(info_fista.time,info_fista.ObjErr,'b-','linewidth',2); hold on
end

if runCGD
loglog(info_cg.time,info_cg.ObjErr,'b-','linewidth',2); hold on
end

if runCGD_ls
loglog(info_cg_ls.time,info_cg_ls.ObjErr,'b--','linewidth',2); hold on
end

h = legend(legendStr,'Location', 'Best');
grid on;
set(h,'Interpreter','latex')
xlabel('time (sec)')
title('$F(x_k) - F^*$','Interpreter','Latex')
set(gca,'fontsize',25)
shg

fname = sprintf('ObjErrTime_p%d_n%d_s%d_sizeGrp%d_randomGrps%d_adaptol%d_tol%s',p,n,s,sizeGrp,randomGrps,adapt_tol,tol_str)
if print_flag
    print(gcf,'-dpdf','-r600',fname)
end
%%
figure

if runGPM
loglog(info_gpm_inf.OptErr2,'r--','linewidth',2); hold on
end

if runaccGPM
loglog(info_agpm_inf.OptErr2,'r-','linewidth',2); hold on
end

if runaccGPM_ls
loglog(info_agpm_inf_ls.OptErr2,'r-.','linewidth',2); hold on
end

if runGPM_1
loglog(info_gpm_1.OptErr2,'g--','linewidth',2); hold on
end

if runaccGPM_1
loglog(info_agpm_1.OptErr2,'g-','linewidth',2); hold on
end

if runaccGPM_ls_1
loglog(info_agpm_1_ls.OptErr2,'g-.','linewidth',2); hold on
end


if runFISTA
loglog(info_fista.OptErr2,'b-','linewidth',2); hold on
end

if runCGD
loglog(info_cg.OptErr2,'b-','linewidth',2); hold on
end

if runCGD_ls
loglog(info_cg_ls.OptErr2,'b--','linewidth',2); hold on
end

h = legend(legendStr,'Location', 'Best');
grid on;
set(h,'Interpreter','latex')
xlabel('iterations')
title('$\|x_k - x_{cvx}\|_2$','Interpreter','Latex')
set(gca,'fontsize',25)
shg
fname = sprintf('OptErrIter_p%d_n%d_s%d_sizeGrp%d_randomGrps%d_adaptol%d_tol%s',p,n,s,sizeGrp,randomGrps,adapt_tol,tol_str);
if print_flag
    print(gcf,'-dpdf','-r600',fname)
end
%%
figure

if runGPM
loglog(info_gpm_inf.time/3,info_gpm_inf.OptErr2,'r--','linewidth',2); hold on
end

if runaccGPM
loglog(info_agpm_inf.time/3,info_agpm_inf.OptErr2,'r-','linewidth',2); hold on
end

if runaccGPM_ls
loglog(info_agpm_inf_ls.time,info_agpm_inf.OptErr2,'r-.','linewidth',2); hold on
end

if runGPM_1
loglog(info_gpm_1.time/3,info_gpm_inf.OptErr2,'g--','linewidth',2); hold on
end

if runaccGPM_1
loglog(info_agpm_1.time/3,info_agpm_1.OptErr2,'g-','linewidth',2); hold on
end

if runaccGPM_ls_1
loglog(info_agpm_1_ls.time,info_agpm_1.OptErr2,'g-.','linewidth',2); hold on
end

if runFISTA
loglog(info_fista.time,info_fista.OptErr2,'b-','linewidth',2); hold on
end

if runCGD
loglog(info_cg.time,info_cg.OptErr2,'b-','linewidth',2); hold on
end

if runCGD_ls
loglog(info_cg_ls.time,info_cg_ls.OptErr2,'b--','linewidth',2); hold on
end

h = legend(legendStr,'Location', 'Best');
grid on;
set(h,'Interpreter','latex')
xlabel('time (sec)')
title('$\|x_k - x_{cvx}\|_2$','Interpreter','Latex')
set(gca,'fontsize',25)
shg
fname = sprintf('OptErrTime_p%d_n%d_s%d_sizeGrp%d_randomGrps%d_adaptol%d_tol%s',p,n,s,sizeGrp,randomGrps,adapt_tol,tol_str);
if print_flag
    print(gcf,'-dpdf','-r600',fname)
end
%%
figure; plot(info_agpm_inf.itertime,'r-'); hold on
%plot(info_agpm_1.itertime,'g-'); hold on
plot(info_fista.itertime,'b-')
%plot(info_cg.itertime,'b-')
h = legend(legendStr,'Location', 'Best');
grid on;
set(h,'Interpreter','latex')
xlabel('iterations')
title('Iteration time','Interpreter','Latex')
set(gca,'fontsize',25)
shg

%%

