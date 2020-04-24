clear
addpath( genpath('../../../dde_biftool_v3.1.1'));
myfuncs
load('Two_par_I_Kappa_vars.mat', 'stst_br2')
%% Varying Tau vs. I 
parameters_tau=stst_br2.point(105).parameter;
x2=stst_br2.point(105).x; %initial guess to find an equilibrium
ind_Tau=17;
[stst_br4,suc]=SetupStst(funcs,...
    'parameter',parameters_tau,'x',x2,...
    'contpar',ind_Tau,'step',0.01);
stst_br4.method.continuation.plot=0;
stst_br4.parameter.max_step=[17 1];
stst_br4.parameter.min_bound=[17 0];
stst_br4.parameter.max_bound=[17 100];
[stst_br4,s,f,r]=br_contn(funcs,stst_br4,200);
stst_br4=br_rvers(stst_br4)
[stst_br4,s,f,r]=br_contn(funcs,stst_br4,90);
stst_br4=br_stabl(funcs, stst_br4,0,0)
%%
k=1;
for i=10:20:140
    parameters_cont=stst_br4.point(i).parameter;
    x2=stst_br4.point(i).x; %initial guess to find an equilibrium
    cont_par=13;
    [stst_br5,suc]=SetupStst(funcs,...
        'parameter',parameters_cont,'x',x2,...
        'contpar',cont_par,'step',0.001,'plot_progress',0);
    stst_br5.method.stability.minimal_real_part=-3;
    stst_br5.parameter.max_step=[13 0.25; 17 1];
    stst_br5.parameter.min_bound=[13 -60; 17 0];
    stst_br5.parameter.max_bound=[13 60; 17 100];
    figure(7)
    [stst_br5,s,f,r]=br_contn(funcs,stst_br5,500);
    stst_br5=br_rvers(stst_br5)
    [stst_br5,s,f,r]=br_contn(funcs,stst_br5,1000);
    stst_br5=br_stabl(funcs, stst_br5,0,0)
    [xm,ym]=df_measr(0,stst_br5);
    br_splot(stst_br5,xm,ym)
    branchesIT(k).br=stst_br5;
    k=k+1;
end     
title('Bifurcation diagram in the (I,V) plane for 0<=\tau<=100', 'FontSize', 18)
xlabel('I', 'FontSize', 18);ylabel('V', 'FontSize', 18);
%% 2 parameter continuation of hopf for i and tau
nunst=GetStability(branchesIT(1).br,'funcs',funcs);
indhopf=find(abs(diff(nunst))==2);
[branch2,suc]=SetupHopf(funcs,branchesIT(1).br,indhopf,'contpar',[13,17],...
    'dir',17,'step',0.001,'plot',0);
branch2.parameter.max_step=[13 0.2; 17 1];    
branch2.parameter.min_bound=[13 0; 17 0];
branch2.parameter.max_bound=[13 100; 17 100];
branch2.method.stability.minimal_real_part=-2;
[branch2,s,f,r]=br_contn(funcs,branch2,10000);
branch2=br_rvers(branch2)
[branch2,s,f,r]=br_contn(funcs,branch2,7200);
%%
figure(8)
for j=1:4
    [xm,ym]=df_measr(0, branchesIT(1).br);
    zm=ym;
    ym=xm;
    ym.col=17;
    subplot(2,2,j)
    br_plot3(branch2,xm,ym,zm,'k')
    hold on
    plot3(-70,40,0,'kx')
    plot3(-70,40,0,'ko','MarkerFaceColor', 'k')
    xlim([-60,60])
    for i=1:length(branchesIT)
        br_splot3(branchesIT(i).br,xm,ym,zm)
    end
    if j==2
        view(2)
    elseif j==3
        view([0 0]);
    end
    xlabel('I', 'FontSize', 18);ylabel('\tau', 'FontSize', 18);zlabel('V', 'FontSize', 18)
end