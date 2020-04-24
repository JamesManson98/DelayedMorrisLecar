clear all
close all
addpath( genpath('../../../dde_biftool_v3.1.1'));
myfuncs
load('Two_par_Kappa_Tau_vars','branch2','xm','ym','GH_points','HH_points')
ind_Tau=17;
figure(14)
xlabel('\tau', 'FontSize', 18);ylabel('\kappa', 'FontSize', 18);
hold on
branch2positivetau=NaN(1,length(branch2.point));
branch2negativetau=NaN(1,length(branch2.point));
branch2positivekappa=NaN(1,length(branch2.point));
branch2negativekappa=NaN(1,length(branch2.point));
for i=1:length(branch2.point)
    branch2.point(i)=nmfm_hopf(funcs, branch2.point(i)); % normal form computation for hopf    
    if branch2.point(i).nmfm.L1>=0;
        branch2positivetau(i)=branch2.point(i).parameter(17);
        branch2positivekappa(i)=branch2.point(i).parameter(14);
        %indexes_positive=[indexes_positive,i];
    else
        %indexes_negative=[
        branch2negativetau(i)=branch2.point(i).parameter(17);
        branch2negativekappa(i)=branch2.point(i).parameter(14);
    end
end
plot(branch2negativetau(:),branch2negativekappa(:),'color','g')
plot(branch2positivetau(:),branch2positivekappa(:),'color','r')
for i=1:length(GH_points)
    plot(branch2.point(GH_points(i)).parameter(17),branch2.point(GH_points(i)).parameter(14),'Marker','x','color', 'k')
end
for i=1:length(HH_points)
    plot(branch2.point(HH_points(i)).parameter(17),branch2.point(HH_points(i)).parameter(14),'Marker','o','color', 'k')
end
legend('Stable Hopf continuation','Unstable Hopf continuation','Generalised Hopf Points', 'HH points')
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
% for i=[4700 4850 4950 5100]
for i=[301,1501,GH_points(3)+5,GH_points(3)+10,GH_points(3)+15,GH_points(4)+10,GH_points(5)+10,GH_points(6)+30,GH_points(7)+20,GH_points(8)+20,GH_points(10)+30,GH_points(12)+40,GH_points(14)+51]
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
    if i==GH_points(3)+5||i==GH_points(3)+10||i==GH_points(3)+15
        branch4.parameter.min_bound=[17 branch2.point(i).parameter(17)-10];
        branch4.parameter.max_bound=[17 branch2.point(i).parameter(17)+10];
    else
        branch4.parameter.min_bound=[17 branch2.point(i).parameter(17)-1];
        branch4.parameter.max_bound=[17 branch2.point(i).parameter(17)+1];
    end
    [branch4,s,f,r]=br_contn(funcs,branch4,500);
    br_rvers(branch4);
    [branch4,s,f,r]=br_contn(funcs,branch4,500);
    branch4=br_stabl(funcs,branch4,0,0);
    [xm,ym]=df_measr(0,branch4);
    ym.col=14;
    ym.field='parameter';
    xm.col=17;
    br_splot(branch4,xm,ym,0);
    if i==301
        branch4a=branch4;
    elseif i==1501
        branch4b=branch4;
    elseif i==GH_points(6)+30
        branch4d=branch4;
    elseif i==GH_points(3)+5
        branch4e=branch4;
    elseif i==GH_points(3)+10
        branch4f=branch4;
    elseif i==GH_points(3)+15
        branch4g=branch4;
    elseif i==GH_points(4)+10
        branch4c=branch4;
    elseif i==GH_points(5)+10
        branch4h=branch4;
    elseif i==GH_points(7)+20
        branch4i=branch4;
    elseif i==GH_points(8)+20
        branch4j=branch4;
    elseif i==GH_points(10)+30
        branch4k=branch4;
    elseif i==GH_points(12)+40
        branch4l=branch4;
    elseif i==GH_points(14)+51
        branch4m=branch4;
    end
end
%%
plot(branch2.point(4700).parameter(17) , branch2.point(4700).parameter(14),'o', 'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10)
plot(branch2.point(4850).parameter(17) , branch2.point(4850).parameter(14),'o', 'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',10)
plot(branch2.point(4950).parameter(17) , branch2.point(4950).parameter(14),'o', 'MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',10)
plot(branch2.point(5100).parameter(17) , branch2.point(5100).parameter(14),'o', 'MarkerEdgeColor','k',...
                'MarkerFaceColor','k',...
                'MarkerSize',10)