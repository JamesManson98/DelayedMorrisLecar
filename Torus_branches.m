clear all
close all
addpath( genpath('../../../dde_biftool_v3.1.1'));
myfuncs
load('Two_par_Kappa_Tau_vars','branch2','xm','ym','GH_points','HH_points')
ind_Tau=17;
branch2positivetau=NaN(1,length(branch2.point));
branch2negativetau=NaN(1,length(branch2.point));
branch2positivekappa=NaN(1,length(branch2.point));
branch2negativekappa=NaN(1,length(branch2.point));
for i=1:length(branch2.point)
    branch2.point(i)=nmfm_hopf(funcs, branch2.point(i)); % normal form computation for hopf    
    if branch2.point(i).nmfm.L1>=0;
        branch2positivetau(i)=branch2.point(i).parameter(17);
        branch2positivekappa(i)=branch2.point(i).parameter(14);
    else
        branch2negativetau(i)=branch2.point(i).parameter(17);
        branch2negativekappa(i)=branch2.point(i).parameter(14);
    end
end
%% periodic solutions:
figure()
xlabel('\tau', 'FontSize', 18);ylabel('\kappa', 'FontSize', 18);
hold on
plot(branch2negativetau(:),branch2negativekappa(:),'color','g')
plot(branch2positivetau(:),branch2positivekappa(:),'color','r')
for i=1:length(GH_points)
    plot(branch2.point(GH_points(i)).parameter(17),branch2.point(GH_points(i)).parameter(14),'Marker','x','color', 'k')
end
for i=1:length(HH_points)
    plot(branch2.point(HH_points(i)).parameter(17),branch2.point(HH_points(i)).parameter(14),'Marker','o','color', 'k')
end
legend('Stable Hopf continuation','Unstable Hopf continuation','Generalised Hopf Points', 'HH points')
intervals=18;
degree=3;
%
for i=[HH_points(1)-10, HH_points(3)+10, HH_points(1)-5, HH_points(2)-5, HH_points(5)-5, HH_points(10)+10, HH_points(6)+10, HH_points(7)+10, HH_points(9), HH_points(11)+10, HH_points(4)+10, HH_points(27)+10, HH_points(1)-20, HH_points(19)+10, HH_points(29)+10, HH_points(13)+10, HH_points(15)+10, HH_points(16)+10, HH_points(18)+10, HH_points(30)+10]
%     figure()
    [psol,stepcond]=p_topsol(funcs,branch2.point(i),1e-2,degree,intervals);
    %correct periodic solution guess:
    method=df_mthod(funcs,'psol');
    [psol,success]=p_correc(funcs,psol,14,stepcond,method.point);
    branch4=df_brnch(funcs,ind_Tau,'psol'); % empty branch:
    % make degenerate periodic solution with amplitude zero at hopf point:
    deg_psol=p_topsol(funcs,branch2.point(i),0,degree,intervals);
    deg_psol.parameter(14)=psol.parameter(14);
    % use deg_psol and psol as first two points on branch:
%     deg_psol.mesh=[];
    branch4.point=deg_psol;
%     psol.mesh=[];
    branch4.point(2)=psol;
    branch4.parameter.max_step=[14,0.1; 17,0.1];
    branch4.method.continuation.plot=0;
    branch4.parameter.min_bound=[17 branch2.point(i).parameter(17)-20];
    branch4.parameter.max_bound=[17 branch2.point(i).parameter(17)+20];
    [branch4,s,f,r]=br_contn(funcs,branch4,500);
    br_rvers(branch4);
    [branch4,s,f,r]=br_contn(funcs,branch4,500);
    branch4=br_stabl(funcs,branch4,0,0);
    [xm,ym]=df_measr(0,branch4);
    %br_splot(branch4,xm,ym,0);
    ym.col=14;
    ym.field='parameter';
    xm.col=17;
%     br_splot(branch4,xm,ym,0);
    figure(1)
    br_splot(branch4,xm,ym,0);
    if i==HH_points(1)-10
        torus_branch(1).br=branch4;
    elseif i==HH_points(3)+10
        torus_branch(3).br=branch4;
    elseif i==HH_points(1)-5
        BRANCH_TORUS_2=branch4;
    elseif i==HH_points(2)-5
        BRANCH_TORUS=branch4;
    elseif i==HH_points(5)-5
        BRANCH_TORUS_3=branch4;
    elseif i==HH_points(10)+10
        TEST_BRANCH_5=branch4;
    elseif i==HH_points(6)+10
        test_5_branch=branch4;
    elseif i==HH_points(7)+10
        more_test_br=branch4;
    elseif i==HH_points(9)
        test_branch=branch4;
    elseif i==HH_points(11)+10
        test_branch_11=branch4;
    elseif i==HH_points(4)+10
        test_branch_4=branch4;
    elseif i==HH_points(27)+10
        test_branch_27=branch4;
    elseif i==HH_points(1)-20
        HH_points1_minus20_branch=branch4;
    elseif i==HH_points(19)+10
        test_branch_19=branch4;
    elseif i==HH_points(29)+10
        test_branch_29=branch4;
    elseif i==HH_points(13)+10
        test_branch_13=branch4;
    elseif i==HH_points(15)+10
        test_branch_15=branch4;
    elseif i==HH_points(16)+10
        test_branch_16=branch4;
    elseif i==HH_points(18)+10
        test_branch_18=branch4;
    elseif i==HH_points(30)+10
        test_branch_30=branch4;
    end
end
%%
for j=27:27
    j
    for i=799%HH_points(j)-10:2:HH_points(j)+10%+10
        figure()
        [psol,stepcond]=p_topsol(funcs,branch2.point(i),1e-2,degree,intervals);
        %correct periodic solution guess:
        method=df_mthod(funcs,'psol');
        [psol,success]=p_correc(funcs,psol,14,stepcond,method.point);
        branch4=df_brnch(funcs,ind_Tau,'psol'); % empty branch:
        % make degenerate periodic solution with amplitude zero at hopf point:
        deg_psol=p_topsol(funcs,branch2.point(i),0,degree,intervals);
        deg_psol.parameter(14)=psol.parameter(14);
        % use deg_psol and psol as first two points on branch:
    %     deg_psol.mesh=[];
        branch4.point=deg_psol;
    %     psol.mesh=[];
        branch4.point(2)=psol;
        branch4.method.continuation.plot=0;
        branch4.parameter.min_bound=[17 branch2.point(i).parameter(17)-20];
        branch4.parameter.max_bound=[17 branch2.point(i).parameter(17)+20];
        branch4.parameter.max_step=[14,0.01; 17,0.01];

        [branch4,s,f,r]=br_contn(funcs,branch4,500);
        br_rvers(branch4);
        [branch4,s,f,r]=br_contn(funcs,branch4,500);
        branch4=br_stabl(funcs,branch4,0,0);
        [xm,ym]=df_measr(0,branch4);
        br_splot(branch4,xm,ym,0);
        % change max step for BRANCH_TORUS_2 AND BRANCH_TORUS_3 (0.1)
        figure(1)
        ym.col=14;
        ym.field='parameter';
        xm.col=17;
        br_splot(branch4,xm,ym,0);
    end
end