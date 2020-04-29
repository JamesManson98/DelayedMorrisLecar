%%generates branches to initialise torus continuation
clear all
close all
addpath( genpath('../../../dde_biftool_v3.1.1'));
myfuncs
load('Two_par_Kappa_Tau_vars','hopf_br_tk','xm','ym','GH_points','HH_points','ind_I','ind_Kappa','ind_Tau')
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
%% periodic solutions:
intervals=18;
degree=3;
for i=[7811,HH_points(1)-10, HH_points(3)+10, HH_points(1)-5, HH_points(2)-5, HH_points(5)-5, HH_points(10)+10, HH_points(6)+10, HH_points(7)+10, HH_points(9), HH_points(11)+10, HH_points(4)+10, HH_points(27)+10, GH_points(4)+200, HH_points(19)+10, HH_points(29)+10, HH_points(13)+10, HH_points(15)+10, HH_points(16)+10, HH_points(18)+10, HH_points(30)+10]
%     figure()
    [psol,stepcond]=p_topsol(funcs,hopf_br_tk.point(i),1e-2,degree,intervals);
    %correct periodic solution guess:
    method=df_mthod(funcs,'psol');
    [psol,success]=p_correc(funcs,psol,14,stepcond,method.point);
    temp_branch=df_brnch(funcs,ind_Tau,'psol'); % empty branch:
    % make degenerate periodic solution with amplitude zero at hopf point:
    deg_psol=p_topsol(funcs,hopf_br_tk.point(i),0,degree,intervals);
    deg_psol.parameter(14)=psol.parameter(14);
    % use deg_psol and psol as first two points on branch:
%     deg_psol.mesh=[];
    temp_branch.point=deg_psol;
%     psol.mesh=[];
    temp_branch.point(2)=psol;
    temp_branch.parameter.max_step=[14,0.1; 17,0.1];
    temp_branch.method.continuation.plot=0;
    temp_branch.parameter.min_bound=[17 hopf_br_tk.point(i).parameter(17)-20];
    temp_branch.parameter.max_bound=[17 hopf_br_tk.point(i).parameter(17)+20];
    if i==7811||i==HH_points(1)-20
        temp_branch.parameter.min_bound=[17 hopf_br_tk.point(i).parameter(17)-50];
        temp_branch.parameter.max_bound=[17 hopf_br_tk.point(i).parameter(17)+50];
    end
    [temp_branch,s,f,r]=br_contn(funcs,temp_branch,500);
    br_rvers(temp_branch);
    [temp_branch,s,f,r]=br_contn(funcs,temp_branch,500);
    temp_branch=br_stabl(funcs,temp_branch,0,0);
    [xm,ym]=df_measr(0,temp_branch);
    %br_splot(branch4,xm,ym,0);
    ym.col=14;
    ym.field='parameter';
    xm.col=17;
%     br_splot(branch4,xm,ym,0);
    figure(1)
    br_splot(temp_branch,xm,ym,0);
    if i==HH_points(1)-10
        init_torus(1)=temp_branch;
    elseif i==HH_points(3)+10
        init_torus(3)=temp_branch;
    elseif i==HH_points(1)-5
        init_torus(2)=temp_branch;
    elseif i==HH_points(2)-5
        init_torus(4)=temp_branch;
    elseif i==HH_points(5)-5
        init_torus(5)=temp_branch;
    elseif i==HH_points(10)+10
        init_torus(6)=temp_branch;
    elseif i==HH_points(6)+10
        init_torus(7)=temp_branch;
    elseif i==HH_points(7)+10
        init_torus(8)=temp_branch;
    elseif i==HH_points(9)
        init_torus(9)=temp_branch;
    elseif i==HH_points(11)+10
        init_torus(10)=temp_branch;
    elseif i==HH_points(4)+10
        init_torus(11)=temp_branch;
    elseif i==HH_points(27)+10
        init_torus(12)=temp_branch;
    elseif i==GH_points(4)+200
        init_torus(13)=temp_branch;
    elseif i==HH_points(19)+10
        init_torus(14)=temp_branch;
    elseif i==HH_points(29)+10
        init_torus(15)=temp_branch;
    elseif i==HH_points(13)+10
        init_torus(16)=temp_branch;
    elseif i==HH_points(15)+10
        init_torus(17)=temp_branch;
    elseif i==HH_points(16)+10
        init_torus(18)=temp_branch;
    elseif i==HH_points(18)+10
        init_torus(19)=temp_branch;
    elseif i==HH_points(30)+10
        init_torus(20)=temp_branch;
    elseif i==7811
        init_torus(21)=temp_branch;
    end
end