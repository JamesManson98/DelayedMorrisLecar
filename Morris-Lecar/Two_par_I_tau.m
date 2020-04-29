%% code to produce 2 parameter continuations (in I and Tau) and Hopf bifurcations 
clear
addpath( genpath('../../../dde_biftool_v3.1.1'));
myfuncs
load('Two_par_I_Kappa_vars.mat', 'ind_I','ind_Kappa','ind_Tau','stst_br_kappa')
%% Varying Tau
parameters_tau=stst_br_kappa.point(105).parameter;
x2=stst_br_kappa.point(105).x; %initial guess to find an equilibrium
[stst_br_tau,suc]=SetupStst(funcs,...
    'parameter',parameters_tau,'x',x2,...
    'contpar',ind_Tau,'step',0.01);
stst_br_tau.method.continuation.plot=0;
stst_br_tau.parameter.max_step=[17 1];
stst_br_tau.parameter.min_bound=[17 0];
stst_br_tau.parameter.max_bound=[17 100];
[stst_br_tau,s,f,r]=br_contn(funcs,stst_br_tau,200);
stst_br_tau=br_rvers(stst_br_tau);
[stst_br_tau,s,f,r]=br_contn(funcs,stst_br_tau,90);
stst_br_tau=br_stabl(funcs, stst_br_tau,0,0);
%% Varying I
k=1;
for i=10:20:140
    parameters_cont=stst_br_tau.point(i).parameter;
    x2=stst_br_tau.point(i).x; %initial guess to find an equilibrium
    cont_par=13;
    [stst_br_I_tau,suc]=SetupStst(funcs,...
        'parameter',parameters_cont,'x',x2,...
        'contpar',cont_par,'step',0.001,'plot_progress',0);
    stst_br_I_tau.method.stability.minimal_real_part=-3;
    stst_br_I_tau.parameter.max_step=[13 0.25; 17 1];
    stst_br_I_tau.parameter.min_bound=[13 -60; 17 0];
    stst_br_I_tau.parameter.max_bound=[13 60; 17 100];
    [stst_br_I_tau,s,f,r]=br_contn(funcs,stst_br_I_tau,500);
    stst_br_I_tau=br_rvers(stst_br_I_tau);
    [stst_br_I_tau,s,f,r]=br_contn(funcs,stst_br_I_tau,1000);
    stst_br_I_tau=br_stabl(funcs, stst_br_I_tau,0,0);
    [xm,ym]=df_measr(0,stst_br_I_tau);
    branchesIT(k).br=stst_br_I_tau;
    k=k+1;
end     
%% 2 parameter continuation of hopf for i and tau
nunst=GetStability(branchesIT(1).br,'funcs',funcs);
indhopf=find(abs(diff(nunst))==2);
[hopf_br_IT,suc]=SetupHopf(funcs,branchesIT(1).br,indhopf,'contpar',[13,17],...
    'dir',17,'step',0.001,'plot',0);
hopf_br_IT.parameter.max_step=[13 0.2; 17 1];    
hopf_br_IT.parameter.min_bound=[13 0; 17 0];
hopf_br_IT.parameter.max_bound=[13 100; 17 100];
hopf_br_IT.method.stability.minimal_real_part=-2;
[hopf_br_IT,s,f,r]=br_contn(funcs,hopf_br_IT,10000);
hopf_br_IT=br_rvers(hopf_br_IT);
[hopf_br_IT,s,f,r]=br_contn(funcs,hopf_br_IT,7200);
%% 3d plot
figure(5)
for j=1:4
    [xm,ym]=df_measr(0, branchesIT(1).br);
    zm=ym;
    ym=xm;
    ym.col=ind_Tau;
    subplot(2,2,j)
    br_plot3(hopf_br_IT,xm,ym,zm,'k')
    hold on
    plot3(-70,40,0,'kx')
    plot3(-70,40,0,'ko','MarkerFaceColor', 'k')
    xlim([-60,60])
    zlim([-90,10])
    for i=1:length(branchesIT)
        br_splot3(branchesIT(i).br,xm,ym,zm)
    end
    if j==2
        view(2)
    elseif j==3
        view([0 0]);
    elseif j==4
        view([90 0]);
    end
    xlabel('I', 'FontSize', 18);ylabel('\tau', 'FontSize', 18);zlabel('V', 'FontSize', 18)
end