clear all
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

%% periodic solutions
intervals=21;
degree=3;
k=1;
colors=['g','r','b','y'];
for i=[4700 4850 4950 5050]
    [psol,stepcond]=p_topsol(funcs,branch2.point(i),1e-2,degree,intervals);
    method=df_mthod(funcs,'psol');
    [psol,success]=p_correc(funcs,psol,17,stepcond,method.point);
    br4(k).pr=df_brnch(funcs,ind_Tau,'psol');
    deg_psol=p_topsol(funcs,branch2.point(i),0,degree,intervals);
    %deg_psol.parameter(14)=psol.parameter(14);
    br4(k).pr.point=deg_psol;
    br4(k).pr.point(2)=psol;
    br4(k).pr.parameter.max_step=[17,0.1];
    br4(k).pr.method.continuation.plot=0;
    [br4(k).pr,s,f,r]=br_contn(funcs,br4(k).pr,1000);
    if i==5050
        [br4(k).pr,s,f,r]=br_contn(funcs,br4(k).pr,150);
    end
    br4(k).pr=br_stabl(funcs,br4(k).pr,0,0);
    k=k+1;
end
%% steady states
k=1;
for i=[4700 4850 4950 5050]
    [br4(k).st,suc]=SetupStst(funcs,...
    'parameter',branch2.point(i).parameter,'x',branch2.point(i).x,...
    'contpar',17,'corpar',17,'step',0.001);
    br4(k).st.parameter.max_step=[17,0.1];
    br4(k).st.parameter.max_bound=[17 100];
    br4(k).st.parameter.min_bound=[17 0];
    br4(k).st.method.continuation.plot=0;

    [br4(k).st,s,f,r]=br_contn(funcs,br4(k).st,10000);
    br4(k).st=br_rvers(br4(k).st)
    [br4(k).st,s,f,r]=br_contn(funcs,br4(k).st,10000);
    br4(k).st=br_stabl(funcs,br4(k).st,0,0);
    k=k+1;
end
%%
[xmp,ymp]=df_measr(0,br4(1).pr);
zmp=xmp;
zmp.col=14;

[xms,yms]=df_measr(0,br4(1).st);
zms=xms;
zms.col=14;
colors=['g','r','b','y'];

for k=1:4
    subplot(2,2,k)
    plot(1,10,'o', 'MarkerEdgeColor','k',...
     'MarkerFaceColor',colors(k),...
     'MarkerSize',10)
    hold on
    ymp.col='max';
    br_plot(br4(k).pr,xmp,ymp,'k')
    ymp.col='min';
    br_plot(br4(k).pr,xmp,ymp,'k')    
    br_splot_l1(br4(k).st,xms,yms);
    ylim([-35,25])
    xlabel('\tau')
    ylabel('V')
end

%%
figure()
[hbranch,suc]=SetupHopf(funcs,br4(2).st,511,'contpar',[14 17],'correc',1,'dir',[14])
hbranch.parameter.max_step=[14,0.1];
hbranch=br_contn(funcs,hbranch,2000);
hbranch=br_rvers(hbranch);
hbranch=br_contn(funcs,hbranch,2000);

br_plot(branch2,ym, xm,'k')
br_plot(hbranch,ym, xm,'r')
%%
figure()
for i=1:1042
    s(i)=p_dststb(br4(2).st.point(i));
end
s2=GetStability(br4(2).st);
[xmp,ymp]=df_measr(1,br4(2).st);
ymp.subfield='l0';
br_plot(br4(2).st,[],ymp)

[xmp,ymp]=df_measr(1,br4(2).st);
ymp.subfield='l1';
br_plot(br4(2).st,[],ymp)
plot(s-2)
plot(s2-1)
plot([0 1042],[0 0],'k')
ylim([-0.1 0.1])
