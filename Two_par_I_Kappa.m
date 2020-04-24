clear
addpath( genpath('../../../dde_biftool_v3.1.1'));
myfuncs
load('One_par_I_vars.mat','ind_I', 'stst_br')
%% Varying Kappa
parameters_kappa=stst_br.point(21).parameter;
x1=stst_br.point(21).x; %initial guess to find an equilibrium
ind_Kappa=14;
[stst_br2,suc]=SetupStst(funcs,...
    'parameter',parameters_kappa,'x',x1,...
    'contpar',ind_Kappa,'step',0.01);
stst_br2.method.continuation.plot=0;              
stst_br2.parameter.max_step=[14 1];
stst_br2.parameter.min_bound=[14 0];
stst_br2.parameter.max_bound=[14 100];
[stst_br2,s,f,r]=br_contn(funcs,stst_br2,200);
stst_br2=br_rvers(stst_br2)
[stst_br2,s,f,r]=br_contn(funcs,stst_br2,90);
stst_br2=br_stabl(funcs, stst_br2,0,0)
%%
k=1;
for i=10:15:120
    parameters_cont=stst_br2.point(i).parameter;
    x1=stst_br2.point(i).x; %initial guess to find an equilibrium
    cont_par=13;
    [stst_br3,suc]=SetupStst(funcs,...
        'parameter',parameters_cont,'x',x1,...
        'contpar',cont_par,'step',0.001,'plot_progress',0);
    stst_br3.parameter.max_step=[13 0.5];
    stst_br3.parameter.min_bound=[13 -60];
    stst_br3.parameter.max_bound=[13 60];
    figure(5)
    stst_br3.method.continuation.plot=0;
    [stst_br3,s,f,r]=br_contn(funcs,stst_br3,500);
    stst_br3=br_rvers(stst_br3)
    [stst_br3,s,f,r]=br_contn(funcs,stst_br3,500);
    stst_br3=br_stabl(funcs, stst_br3,0,0)
    [xm,ym]=df_measr(0,stst_br3);
    branchesI(k).br=stst_br3;
    k=k+1;
    br_splot(stst_br3,xm,ym)
end
xlabel('I', 'FontSize', 18);ylabel('V', 'FontSize', 18);
axis([-60,60,-100,20])
%% 2 parameter continuation of hopf
nunst=GetStability(branchesI(3).br,'funcs',funcs);
indhopf=find(abs(diff(nunst))==2)
[branch2,suc]=SetupHopf(funcs,branchesI(3).br,indhopf(2),'contpar',[ind_I,14],...
    'dir',13,'step',0.1,'plot',0);
branch2.parameter.max_step=[13 0.25; 14 1];
branch2.parameter.min_bound=[13 -60; 14 -2];
branch2.parameter.max_bound=[13 60; 14 100];
branch2.method.stability.minimal_real_part=-2;
[branch2,s,f,r]=br_contn(funcs,branch2,1200);
branch2=br_rvers(branch2);
[branch2,s,f,r]=br_contn(funcs,branch2,1200);
% we find the other 2d branch of hopf orbits
[branch2a,suc]=SetupHopf(funcs,branchesI(3).br,indhopf(3),'contpar',[ind_I,14],...
    'dir',13,'step',0.1,'plot',0);
branch2a.parameter.max_step=[13 0.25; 14 1];
branch2a.parameter.min_bound=[13 -60; 14 0];
branch2a.parameter.max_bound=[13 60; 14 100];
branch2a.method.stability.minimal_real_part=-2;
[branch2a,s,f,r]=br_contn(funcs,branch2a,1200);
branch2a=br_rvers(branch2a);
[branch2a,s,f,r]=br_contn(funcs,branch2a,77);
%%
figure(6)
for j=1:4
    [xm,ym]=df_measr(0, branchesI(1).br);
    zm=ym;
    ym=xm;
    ym.col=14;
    subplot(2,2,j)
    br_plot3(branch2,xm,ym,zm,'k')
    hold on
    plot3(-70,40,0,'kx')
    plot3(-70,40,0,'ko','MarkerFaceColor', 'k')
    xlim([-60,60])
    br_plot3(branch2a,xm,ym,zm,'k')
    for i=1:length(branchesI)
        br_splot3(branchesI(i).br,xm,ym,zm)
    end
    if j==2
        view(2)
    elseif j==3
        view([0 0]);
    end
    xlabel('I', 'FontSize', 18);ylabel('\kappa', 'FontSize', 18);zlabel('V', 'FontSize', 18)
end
%% periodic solutions:
% intervals=18;
% degree=3;
% for i=1:100:301
%     [psol,stepcond]=p_topsol(funcs,branch2a.point(i),1e-2,degree,intervals);
%     %correct periodic solution guess:
%     method=df_mthod(funcs,'psol');
%     [psol,success]=p_correc(funcs,psol,4,stepcond,method.point);
%     branch4=df_brnch(funcs,ind_I,'psol'); % empty branch:
%     % make degenerate periodic solution with amplitude zero at hopf point:
%     deg_psol=p_topsol(funcs,branch2a.point(i),0,degree,intervals);
%     % use deg_psol and psol as first two points on branch:
%     deg_psol.mesh=[];
%     branch4.point=deg_psol;
%     psol.mesh=[];
%     branch4.point(2)=psol;
%     [branch4,s,f,r]=br_contn(funcs,branch4,43); 
%     figure(6)
%     hold on
%     [xm,zm]=df_measr(0,branch4);
%     zm.col='max';
%     ym.col=14;
%     ym.field='parameter';
%     br_plot3(branch4,xm,ym,zm,'k')
%     ym.col='min';
%     br_plot3(branch4,xm,ym,zm,'k')
% end