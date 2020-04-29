%% code to produce 2 parameter continuations (in Kappa and Tau) and Hopf bifurcations 
clear
addpath( genpath('../../../dde_biftool_v3.1.1'));
myfuncs
load('Two_par_I_Kappa_vars.mat', 'stst_br_kappa','ind_I','ind_Kappa','ind_Tau')
%% Varying Tau 
parameters_tau_kappa=stst_br_kappa.point(105).parameter;
x3=stst_br_kappa.point(105).x; %initial guess to find an equilibrium
[stst_br_tau,suc]=SetupStst(funcs,...
    'parameter',parameters_tau_kappa,'x',x3,...
    'contpar',ind_Tau,'step',0.01);
stst_br_tau.method.continuation.plot=0;
stst_br_tau.parameter.max_step=[ind_Tau 0.5];
stst_br_tau.parameter.min_bound=[ind_Tau 0];
stst_br_tau.parameter.max_bound=[ind_Tau 100];
[stst_br_tau,s,f,r]=br_contn(funcs,stst_br_tau,200);
stst_br_tau=br_rvers(stst_br_tau)
[stst_br_tau,s,f,r]=br_contn(funcs,stst_br_tau,90);
stst_br_tau=br_stabl(funcs, stst_br_tau,0,0)
%% Varying Kappa
k=1;
for i=10:30:198
    parameters_cont=stst_br_tau.point(i).parameter;
    x4=stst_br_tau.point(i).x; %initial guess to find an equilibrium
    [stst_br_tk,suc]=SetupStst(funcs,...
        'parameter',parameters_cont,'x',x4,...
        'contpar',ind_Kappa,'step',0.001,'plot_progress',0);
    stst_br_tk.method.stability.minimal_real_part=-3;
    stst_br_tk.parameter.max_step=[ind_Kappa 0.1; ind_Tau 1];    
    stst_br_tk.parameter.min_bound=[ind_Kappa 0; ind_Tau 0];
    stst_br_tk.parameter.max_bound=[ind_Kappa 100; ind_Tau 100];
    stst_br_tk.method.continuation.plot=0;
    [stst_br_tk,s,f,r]=br_contn(funcs,stst_br_tk,1500);
    stst_br_tk=br_rvers(stst_br_tk)
    [stst_br_tk,s,f,r]=br_contn(funcs,stst_br_tk,1500);
    stst_br_tk=br_stabl(funcs, stst_br_tk,0,0)
    [xm,ym]=df_measr(0,stst_br_tk);
    branchesI_tk(k).br=stst_br_tk;
    k=k+1;
end
%% 2 parameter continuation of hopf for tau and kappa
nunst=GetStability(branchesI_tk(1).br,'funcs',funcs);
indhopf=find(abs(diff(nunst))==2);
[hopf_br_tk,suc]=SetupHopf(funcs,branchesI_tk(1).br,indhopf(1),'contpar',[ind_Kappa,ind_Tau],...
    'dir',14,'step',0.001,'plot',0);
hopf_br_tk.parameter.max_step=[ind_Tau 0.1; ind_Kappa 0.1];    
hopf_br_tk.parameter.min_bound=[ind_Kappa 0; ind_Tau 0];
hopf_br_tk.parameter.max_bound=[ind_Kappa 100; ind_Tau 200];%150
hopf_br_tk.method.stability.minimal_real_part=-2;
[hopf_br_tk,s,f,r]=br_contn(funcs,hopf_br_tk,15000);
hopf_br_tk=br_rvers(hopf_br_tk)
[hopf_br_tk,s,f,r]=br_contn(funcs,hopf_br_tk,15000);
%% 3d plot
figure(6)
for j=1:4
    [xm,ym]=df_measr(0, branchesI_tk(1).br);
    zm=ym;
    ym=xm;
    ym.col=ind_Tau;
    subplot(2,2,j)
    br_plot3(hopf_br_tk,xm,ym,zm,'k')
    hold on
    plot3(-70,40,0,'kx')
    plot3(-70,40,0,'ko','MarkerFaceColor', 'k')
    xlim([0,100]);
    ylim([0,100]);
    zlim([0,8]);
    xlabel('\kappa', 'FontSize', 18);ylabel('\tau', 'FontSize', 18);zlabel('V', 'FontSize', 18);
    for i=1:length(branchesI_tk)
        branchesI_tk(i).br.method.continuation.plot=0;
        branchesI_tk(i).br=br_rvers(branchesI_tk(i).br);
        branchesI_tk(i).br=br_contn(funcs,branchesI_tk(i).br,1000);
        branchesI_tk(i).br=br_stabl(funcs,branchesI_tk(i).br,0,0);
        br_splot3(branchesI_tk(i).br,xm,ym,zm)
    end
    if j==2
        view(gca,[90 270]);
    elseif j==3
        view([0 0]);
    elseif j==4
        view([90 0]);
    end
end
%% Hopf-Hopf and Generalised-Hopf points
figure(7)
xm.field='parameter';
xm.col=17;
ym.field='parameter';
ym.col=14;
br_plot(hopf_br_tk,xm, ym,'k')
hold on
hopf_br_tk = br_stabl(funcs,hopf_br_tk,0,0);
hopf_br_tk = br_bifdet(funcs,hopf_br_tk);
flagged_hopfs=br_getflags(hopf_br_tk);
GH_points=flagged_hopfs(5,:);
GH_points(GH_points==0)=[];
HH_points=flagged_hopfs(6,:);
HH_points(HH_points==0)=[];
for i=1:length(GH_points)
    plot(hopf_br_tk.point(GH_points(i)).parameter(ind_Tau),hopf_br_tk.point(GH_points(i)).parameter(ind_Kappa),'Marker','x','color', 'r')
end
for i=1:length(HH_points)
    plot(hopf_br_tk.point(HH_points(i)).parameter(ind_Tau),hopf_br_tk.point(HH_points(i)).parameter(ind_Kappa),'Marker','o','color', 'r')
end
plot(hopf_br_tk.point(4700).parameter(ind_Tau) , hopf_br_tk.point(4700).parameter(ind_Kappa),'o', 'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10)
plot(hopf_br_tk.point(4850).parameter(ind_Tau) , hopf_br_tk.point(4850).parameter(ind_Kappa),'o', 'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',10)
plot(hopf_br_tk.point(4950).parameter(ind_Tau) , hopf_br_tk.point(4950).parameter(ind_Kappa),'o', 'MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',10)
plot(hopf_br_tk.point(5050).parameter(ind_Tau) , hopf_br_tk.point(5050).parameter(ind_Kappa),'o', 'MarkerEdgeColor','k',...
                'MarkerFaceColor','y',...
                'MarkerSize',10)
h(1)=plot(1000,1000, 'k');
h(2)=plot(1000,1000, 'Marker', 'x', 'color', 'r', 'linestyle', 'none');
h(3)=plot(1000,1000, 'Marker', 'o', 'color', 'r', 'linestyle', 'none');
xlim([0,100]);
ylim([0,100]);
xlabel('\tau', 'FontSize', 18);ylabel('\kappa', 'FontSize', 18);
legend(h, {'Hopf continuation', 'GH points', 'HH points'},'FontSize', 18)