%%generates branches to initialise fold of POs continuation
clear all
close all
addpath( genpath('../../../dde_biftool_v3.1.1'));
myfuncs
load('Two_par_Kappa_Tau_vars','hopf_br_tk','xm','ym','GH_points','HH_points','ind_I','ind_Kappa','ind_Tau')
figure(14)
xlabel('\tau', 'FontSize', 18);ylabel('\kappa', 'FontSize', 18);
hold on
hopf_br_tk_unstable_tau=NaN(1,length(hopf_br_tk.point));
hopf_br_tk_stable_tau=NaN(1,length(hopf_br_tk.point));
hopf_br_tk_unstable_kappa=NaN(1,length(hopf_br_tk.point));
hopf_br_tk_stable_kappa=NaN(1,length(hopf_br_tk.point));
for i=1:length(hopf_br_tk.point)
    hopf_br_tk.point(i)=nmfm_hopf(funcs, hopf_br_tk.point(i)); % normal form computation for hopf    
    if hopf_br_tk.point(i).nmfm.L1>=0
        hopf_br_tk_unstable_tau(i)=hopf_br_tk.point(i).parameter(ind_Tau);
        hopf_br_tk_unstable_kappa(i)=hopf_br_tk.point(i).parameter(ind_Kappa);
    else
        hopf_br_tk_stable_tau(i)=hopf_br_tk.point(i).parameter(ind_Tau);
        hopf_br_tk_stable_kappa(i)=hopf_br_tk.point(i).parameter(ind_Kappa);
    end
end
plot(hopf_br_tk_stable_tau(:),hopf_br_tk_stable_kappa(:),'color','g')
plot(hopf_br_tk_unstable_tau(:),hopf_br_tk_unstable_kappa(:),'color','r')
for i=1:length(GH_points)
    plot(hopf_br_tk.point(GH_points(i)).parameter(ind_Tau),hopf_br_tk.point(GH_points(i)).parameter(ind_Kappa),'Marker','x','color', 'k')
end
for i=1:length(HH_points)
    plot(hopf_br_tk.point(HH_points(i)).parameter(ind_Tau),hopf_br_tk.point(HH_points(i)).parameter(ind_Kappa),'Marker','o','color', 'k')
end
legend('Stable Hopf continuation','Unstable Hopf continuation','Generalised Hopf Points', 'HH points')
%% periodic solutions:
figure()
xlabel('\tau', 'FontSize', 18);ylabel('\kappa', 'FontSize', 18);
hold on
plot(hopf_br_tk_stable_tau(:),hopf_br_tk_stable_kappa(:),'color','g')
plot(hopf_br_tk_unstable_tau(:),hopf_br_tk_unstable_kappa(:),'color','r')
for i=1:length(GH_points)
    plot(hopf_br_tk.point(GH_points(i)).parameter(ind_Tau),hopf_br_tk.point(GH_points(i)).parameter(ind_Kappa),'Marker','x','color', 'k')
end
for i=1:length(HH_points)
    plot(hopf_br_tk.point(HH_points(i)).parameter(ind_Tau),hopf_br_tk.point(HH_points(i)).parameter(ind_Kappa),'Marker','o','color', 'k')
end
legend('Stable Hopf continuation','Unstable Hopf continuation','Generalised Hopf Points', 'HH points')
intervals=18;
degree=3;
for i=[301,1501,GH_points(3)+5,GH_points(4)+10,GH_points(5)+10,GH_points(6)+30,GH_points(7)+20,GH_points(8)+20,GH_points(10)+30,GH_points(12)+40,GH_points(14)+51]
    [psol,stepcond]=p_topsol(funcs,hopf_br_tk.point(i),1e-2,degree,intervals);
    %correct periodic solution guess:
    method=df_mthod(funcs,'psol');
    [psol,success]=p_correc(funcs,psol,ind_Kappa,stepcond,method.point);
    temp_branch=df_brnch(funcs,ind_Tau,'psol'); % empty branch:
    % make degenerate periodic solution with amplitude zero at hopf point:
    deg_psol=p_topsol(funcs,hopf_br_tk.point(i),0,degree,intervals);
    deg_psol.parameter(14)=psol.parameter(ind_Kappa);
    % use deg_psol and psol as first two points on branch:
%     deg_psol.mesh=[];
    temp_branch.point=deg_psol;
%     psol.mesh=[];
    temp_branch.point(2)=psol;
    temp_branch.method.continuation.plot=0;
    if i==GH_points(3)+5||i==GH_points(3)+10||i==GH_points(3)+15
        temp_branch.parameter.min_bound=[ind_Tau hopf_br_tk.point(i).parameter(ind_Tau)-10];
        temp_branch.parameter.max_bound=[ind_Tau hopf_br_tk.point(i).parameter(ind_Tau)+10];
    else
        temp_branch.parameter.min_bound=[ind_Tau hopf_br_tk.point(i).parameter(ind_Tau)-1];
        temp_branch.parameter.max_bound=[ind_Tau hopf_br_tk.point(i).parameter(ind_Tau)+1];
    end
    [temp_branch,s,f,r]=br_contn(funcs,temp_branch,500);
    br_rvers(temp_branch);
    [temp_branch,s,f,r]=br_contn(funcs,temp_branch,500);
    temp_branch=br_stabl(funcs,temp_branch,0,0);
    [xm,ym]=df_measr(0,temp_branch);
    ym.col=ind_Kappa;
    ym.field='parameter';
    xm.col=ind_Tau;
    br_splot(temp_branch,xm,ym,0);
    if i==301
        init_PO_fold(1).br=temp_branch;
    elseif i==1501
        init_PO_fold(2).br=temp_branch;
    elseif i==GH_points(4)+10
        init_PO_fold(3).br=temp_branch;        
    elseif i==GH_points(6)+30
        init_PO_fold(4).br=temp_branch;
    elseif i==GH_points(3)+5
        init_PO_fold(5).br=temp_branch;
    elseif i==GH_points(5)+10
        init_PO_fold(6).br=temp_branch;
    elseif i==GH_points(7)+20
        init_PO_fold(7).br=temp_branch;
    elseif i==GH_points(8)+20
        init_PO_fold(8).br=temp_branch;
    elseif i==GH_points(10)+30
        init_PO_fold(9).br=temp_branch;
    elseif i==GH_points(12)+40
        init_PO_fold(10).br=temp_branch;
    elseif i==GH_points(14)+51
        init_PO_fold(11).br=temp_branch;
    end
end