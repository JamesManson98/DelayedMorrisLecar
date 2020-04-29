%% code to produce 2 parameter continuations (in I and Kappa) and Hopf bifurcations 
clear
addpath( genpath('../../../dde_biftool_v3.1.1'));
myfuncs
load('One_par_I_vars.mat','ind_I','ind_Kappa','ind_Tau', 'stst_br_I')
%% Varying Kappa
parameters_kappa=stst_br_I.point(21).parameter;
x1=stst_br_I.point(21).x; %initial guess to find an equilibrium
[stst_br_kappa,suc]=SetupStst(funcs,...
    'parameter',parameters_kappa,'x',x1,...
    'contpar',ind_Kappa,'step',0.01);
stst_br_kappa.method.continuation.plot=0;              
stst_br_kappa.parameter.max_step=[ind_Kappa 1];
stst_br_kappa.parameter.min_bound=[ind_Kappa 0];
stst_br_kappa.parameter.max_bound=[ind_Kappa 100];
[stst_br_kappa,s,f,r]=br_contn(funcs,stst_br_kappa,200);
stst_br_kappa=br_rvers(stst_br_kappa)
[stst_br_kappa,s,f,r]=br_contn(funcs,stst_br_kappa,90);
stst_br_kappa=br_stabl(funcs, stst_br_kappa,0,0)
%% Varying I
k=1;
for i=10:15:120
    parameters_cont=stst_br_kappa.point(i).parameter;
    x1=stst_br_kappa.point(i).x; %initial guess to find an equilibrium
    cont_par=ind_I;
    [stst_br_I_kappa,suc]=SetupStst(funcs,...
        'parameter',parameters_cont,'x',x1,...
        'contpar',cont_par,'step',0.001,'plot_progress',0);
    stst_br_I_kappa.parameter.max_step=[ind_I 0.5];
    stst_br_I_kappa.parameter.min_bound=[ind_I -60];
    stst_br_I_kappa.parameter.max_bound=[ind_I 60];
    stst_br_I_kappa.method.continuation.plot=0;
    [stst_br_I_kappa,s,f,r]=br_contn(funcs,stst_br_I_kappa,500);
    stst_br_I_kappa=br_rvers(stst_br_I_kappa)
    [stst_br_I_kappa,s,f,r]=br_contn(funcs,stst_br_I_kappa,500);
    stst_br_I_kappa=br_stabl(funcs, stst_br_I_kappa,0,0)
    [xm,ym]=df_measr(0,stst_br_I_kappa);
    branches_I_kappa(k).br=stst_br_I_kappa;
    k=k+1;
end
xlabel('I', 'FontSize', 18);ylabel('V', 'FontSize', 18);
axis([-60,60,-100,20])
%% 2 parameter continuation of hopf
nunst=GetStability(branches_I_kappa(3).br,'funcs',funcs);
indhopf=find(abs(diff(nunst))==2);
[hopf_br_I_kappa(1).br,suc]=SetupHopf(funcs,branches_I_kappa(3).br,indhopf(2),'contpar',[ind_I,ind_Kappa],...
    'dir',13,'step',0.1,'plot',0);
hopf_br_I_kappa(1).br.parameter.max_step=[ind_I 0.25; ind_Kappa 1];
hopf_br_I_kappa(1).br.parameter.min_bound=[ind_I -60; ind_Kappa -2];
hopf_br_I_kappa(1).br.parameter.max_bound=[ind_I 60; ind_Kappa 100];
hopf_br_I_kappa(1).br.method.stability.minimal_real_part=-2;
[hopf_br_I_kappa(1).br,s,f,r]=br_contn(funcs,hopf_br_I_kappa(1).br,1200);
hopf_br_I_kappa(1).br=br_rvers(hopf_br_I_kappa(1).br);
[hopf_br_I_kappa(1).br,s,f,r]=br_contn(funcs,hopf_br_I_kappa(1).br,1200);
% we find the other 2d branch of hopf orbits
[hopf_br_I_kappa(2).br,suc]=SetupHopf(funcs,branches_I_kappa(3).br,indhopf(3),'contpar',[ind_I,ind_Kappa],...
    'dir',13,'step',0.1,'plot',0);
hopf_br_I_kappa(2).br.parameter.max_step=[ind_I 0.25; ind_Kappa 1];
hopf_br_I_kappa(2).br.parameter.min_bound=[ind_I -60; ind_Kappa 0];
hopf_br_I_kappa(2).br.parameter.max_bound=[ind_I 60; ind_Kappa 100];
hopf_br_I_kappa(2).br.method.stability.minimal_real_part=-2;
[hopf_br_I_kappa(2).br,s,f,r]=br_contn(funcs,hopf_br_I_kappa(2).br,1200);
hopf_br_I_kappa(2).br=br_rvers(hopf_br_I_kappa(2).br);
[hopf_br_I_kappa(2).br,s,f,r]=br_contn(funcs,hopf_br_I_kappa(2).br,77);
%% 3d plot
figure(4)
for j=1:4
    [xm,ym]=df_measr(0, branches_I_kappa(1).br);
    zm=ym;
    ym=xm;
    ym.col=ind_Kappa;
    subplot(2,2,j)
    br_plot3(hopf_br_I_kappa(1).br,xm,ym,zm,'k')
    hold on
    plot3(-70,40,0,'kx')
    plot3(-70,40,0,'ko','MarkerFaceColor', 'k')
    xlim([-60,60])
    zlim([-90,10])
    br_plot3(hopf_br_I_kappa(2).br,xm,ym,zm,'k')
    for i=1:length(branches_I_kappa)
        br_splot3(branches_I_kappa(i).br,xm,ym,zm)
    end
    if j==2
        view(2)
    elseif j==3
        view([0 0]);
    elseif j==4
        view([90 0]);
    end
    xlabel('I', 'FontSize', 18);ylabel('\kappa', 'FontSize', 18);zlabel('V', 'FontSize', 18)
end