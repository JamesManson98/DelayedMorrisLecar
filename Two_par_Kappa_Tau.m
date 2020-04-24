clear
addpath( genpath('../../../dde_biftool_v3.1.1'));
myfuncs
load('Two_par_I_Kappa_vars.mat', 'stst_br2')
%% Varying Tau vs. Kappa 
parameters_tau_kappa=stst_br2.point(105).parameter;
x3=stst_br2.point(105).x; %initial guess to find an equilibrium
ind_Tau=17;
[stst_br_tk,suc]=SetupStst(funcs,...
    'parameter',parameters_tau_kappa,'x',x3,...
    'contpar',ind_Tau,'step',0.01);
stst_br_tk.method.continuation.plot=0;
stst_br_tk.parameter.max_step=[17 0.5];
stst_br_tk.parameter.min_bound=[17 0];
stst_br_tk.parameter.max_bound=[17 100];
[stst_br_tk,s,f,r]=br_contn(funcs,stst_br_tk,200);
stst_br_tk=br_rvers(stst_br_tk)
[stst_br_tk,s,f,r]=br_contn(funcs,stst_br_tk,90);
stst_br_tk=br_stabl(funcs, stst_br_tk,0,0)
%%
k=1;
figure(9)
for i=10:30:198
    parameters_cont=stst_br_tk.point(i).parameter;
    x4=stst_br_tk.point(i).x; %initial guess to find an equilibrium
    cont_par=14;
    [stst_br5,suc]=SetupStst(funcs,...
        'parameter',parameters_cont,'x',x4,...
        'contpar',cont_par,'step',0.001,'plot_progress',0);
    stst_br5.method.stability.minimal_real_part=-3;
    stst_br5.parameter.max_step=[14 0.1; 17 1];    
    stst_br5.parameter.min_bound=[14 0; 17 0];
    stst_br5.parameter.max_bound=[14 100; 17 100];
    [stst_br5,s,f,r]=br_contn(funcs,stst_br5,1500);
    stst_br5=br_rvers(stst_br5)
    [stst_br5,s,f,r]=br_contn(funcs,stst_br5,1500);
    stst_br5=br_stabl(funcs, stst_br5,0,0)
    [xm,ym]=df_measr(0,stst_br5);
    br_splot(stst_br5,xm,ym)
    branchesI_tk(k).br=stst_br5;
    k=k+1;
end
xlabel('\kappa', 'FontSize', 18);ylabel('V', 'FontSize', 18);
title('Bifurcation diagram in the (\kappa,V) plane for 0<=\tau<=100', 'FontSize', 18)
%% 2 parameter continuation of hopf for tau and kappa
nunst=GetStability(branchesI_tk(1).br,'funcs',funcs);
indhopf=find(abs(diff(nunst))==2);
[branch2,suc]=SetupHopf(funcs,branchesI_tk(1).br,indhopf(1),'contpar',[14,17],...
    'dir',14,'step',0.001,'plot',0);
branch2.parameter.max_step=[17 0.1; 14 0.1];    
%or remove by hand e.g branch2.point(i)=[]
branch2.parameter.min_bound=[14 0; 17 0];
branch2.parameter.max_bound=[14 100; 17 200];%150
branch2.method.stability.minimal_real_part=-2;
[branch2,s,f,r]=br_contn(funcs,branch2,15000);
branch2=br_rvers(branch2)
[branch2,s,f,r]=br_contn(funcs,branch2,15000);
%%
figure(10)
for j=1:4
    [xm,ym]=df_measr(0, branchesI_tk(1).br);
    zm=ym;
    ym=xm;
    ym.col=17;
    subplot(2,2,j)
    br_plot3(branch2,xm,ym,zm,'k')
    hold on
    plot3(-70,40,0,'kx')
    plot3(-70,40,0,'ko','MarkerFaceColor', 'k')
    xlim([0,100]);
    ylim([0,100]);
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
    end
end
%% Criticality of hopf (Found HH, how to find DH and lines?)
figure(11)
br_plot(branch2,ym, xm,'k')
hold on
branch2 = br_stabl(funcs,branch2,0,0);
branch2 = br_bifdet(funcs,branch2);
flagged_hopfs=br_getflags(branch2);
GH_points=flagged_hopfs(5,:);
GH_points(GH_points==0)=[];
HH_points=flagged_hopfs(6,:);
HH_points(HH_points==0)=[];
for i=1:length(GH_points)
    plot(branch2.point(GH_points(i)).parameter(17),branch2.point(GH_points(i)).parameter(14),'Marker','x','color', 'r')
end
for i=1:length(HH_points)
    plot(branch2.point(HH_points(i)).parameter(17),branch2.point(HH_points(i)).parameter(14),'Marker','o','color', 'r')
end
plot(branch2.point(4700).parameter(17) , branch2.point(4700).parameter(14),'o', 'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10)
plot(branch2.point(4850).parameter(17) , branch2.point(4850).parameter(14),'o', 'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',10)
plot(branch2.point(4950).parameter(17) , branch2.point(4950).parameter(14),'o', 'MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',10)
plot(branch2.point(5050).parameter(17) , branch2.point(5050).parameter(14),'o', 'MarkerEdgeColor','k',...
                'MarkerFaceColor','y',...
                'MarkerSize',10)
h(1)=plot(1000,1000, 'k');
h(2)=plot(1000,1000, 'Marker', 'x', 'color', 'r', 'linestyle', 'none');
h(3)=plot(1000,1000, 'Marker', 'o', 'color', 'r', 'linestyle', 'none');
xlim([0,100]);
ylim([0,100]);
xlabel('\tau', 'FontSize', 18);ylabel('\kappa', 'FontSize', 18);
legend(h, {'Hopf continuation', 'GH points', 'HH points'},'FontSize', 18)