clear
close all
addpath( genpath('../../../dde_biftool_v3.1.1'));
myfuncs
load('Two_par_Kappa_Tau_vars','branch2','xm','ym','GH_points','HH_points')
ind_Tau=17;
figure(14)
xlabel('\tau', 'FontSize', 18);ylabel('\kappa', 'FontSize', 18);
hold on
for i=1:length(GH_points)
    plot(branch2.point(GH_points(i)).parameter(17),branch2.point(GH_points(i)).parameter(14),'Marker','x','color', 'r')
end
for i=1:length(HH_points)
    plot(branch2.point(HH_points(i)).parameter(17),branch2.point(HH_points(i)).parameter(14),'Marker','o','color', 'r')
end
br_plot(branch2,ym, xm,'k')
xlim([0,100]);
ylim([0,100]);
%% periodic solutions:
intervals=21;
degree=3;
x=1;
colors=['g','r','b','y'];
for i=[4700 4850 4950 5050]
    [psol,stepcond]=p_topsol(funcs,branch2.point(i),1e-2,degree,intervals);
    %correct periodic solution guess:
    method=df_mthod(funcs,'psol');
    [psol,success]=p_correc(funcs,psol,14,stepcond,method.point);
    branch4=df_brnch(funcs,ind_Tau,'psol'); % empty branch:
    % make degenerate periodic solution with amplitude zero at hopf point:
    deg_psol=p_topsol(funcs,branch2.point(i),0,degree,intervals);
    deg_psol.parameter(14)=psol.parameter(14);
    % use deg_psol and psol as first two points on branch:
    branch4.point=deg_psol;
    branch4.point(2)=psol;
    branch4.parameter.max_step=[17,0.1];
    branch4.method.continuation.plot=0;
    [branch4,s,f,r]=br_contn(funcs,branch4,1000);
    if i==5050
        [branch4,s,f,r]=br_contn(funcs,branch4,150);
    end
    branch4=br_stabl(funcs,branch4,0,0);
    [xm,ym]=df_measr(1,branch4);
    xm.col=17;
    figure(1)
    subplot(2,2,x)
    ylim([-35,25])
%     subplot(121)
%     br_plot(branch4,[],ym,'k')
%     xlabel('Point in continuation')
%     ylabel('Modulus of largest eigenvalue')
%     subplot(122)
    plot(1,10,'o', 'MarkerEdgeColor','k',...
         'MarkerFaceColor',colors(x),...
         'MarkerSize',10)
    hold on
    [xm,ym]=df_measr(0,branch4);
    ym.col='max';
    %branch4=br_stabl(funcs,branch4,0,0);
    br_plot(branch4,xm,ym,'k')
    ym.col='min';
    br_plot(branch4,xm,ym,'k')    
    xlabel('\tau')
    ylabel('V')    
    [stst_br,suc]=SetupStst(funcs,...
    'parameter',branch2.point(i).parameter,'x',branch2.point(i).x,...
    'contpar',17,'step',0.001);
    stst_br.parameter.max_step=[17,0.1];
    stst_br.parameter.max_bound=[17 100];
    stst_br.parameter.min_bound=[17 0];
    stst_br.method.continuation.plot=0;

    [stst_br,s,f,r]=br_contn(funcs,stst_br,10000);
    stst_br=br_rvers(stst_br)
    [stst_br,s,f,r]=br_contn(funcs,stst_br,10000);
    stst_br=br_stabl(funcs, stst_br,0,0);

    [xm,ym]=df_measr(0,stst_br);
    br_splot(stst_br,xm,ym)
    figure(14)
    ym.col=14;
    ym.field='parameter';
    xm.col=17;
    br_splot(branch4,xm,ym,0);
    x=x+1;
end