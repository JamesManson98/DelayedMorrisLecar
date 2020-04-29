%%checking the criticality of the Hopf points
clear all
close all
addpath( genpath('../../../dde_biftool_v3.1.1'));
myfuncs
load('Two_par_Kappa_Tau_vars','hopf_br_tk','xm','ym','GH_points','HH_points')
ind_Tau=17;
%% plotting sub/supercritical Hopfs
figure(14)
xlabel('\tau', 'FontSize', 18);ylabel('\kappa', 'FontSize', 18);
hold on
for i=1:length(GH_points)
    plot(hopf_br_tk.point(GH_points(i)).parameter(17),hopf_br_tk.point(GH_points(i)).parameter(14),'Marker','x','color', 'r')
end
for i=1:length(HH_points)
    plot(hopf_br_tk.point(HH_points(i)).parameter(17),hopf_br_tk.point(HH_points(i)).parameter(14),'Marker','o','color', 'r')
end
br_plot(hopf_br_tk,xm, ym,'k')
xlim([0,100]);
ylim([0,100]);

%% periodic solutions
intervals=21;
degree=3;
k=1;
colors=['g','r','b','y'];
for i=[4700 4850 4950 5050]
    [psol,stepcond]=p_topsol(funcs,hopf_br_tk.point(i),1e-2,degree,intervals);
    method=df_mthod(funcs,'psol');
    [psol,success]=p_correc(funcs,psol,17,stepcond,method.point);
    branch(k).pr=df_brnch(funcs,ind_Tau,'psol');
    deg_psol=p_topsol(funcs,hopf_br_tk.point(i),0,degree,intervals);
    %deg_psol.parameter(14)=psol.parameter(14);
    branch(k).pr.point=deg_psol;
    branch(k).pr.point(2)=psol;
    branch(k).pr.parameter.max_step=[17,0.1];
    branch(k).pr.method.continuation.plot=0;
    [branch(k).pr,s,f,r]=br_contn(funcs,branch(k).pr,1000);
    if i==5050
        [branch(k).pr,s,f,r]=br_contn(funcs,branch(k).pr,150);
    end
    branch(k).pr=br_stabl(funcs,branch(k).pr,0,0);
    k=k+1;
end
%% steady states
k=1;
for i=[4700 4850 4950 5050]
    [branch(k).st,suc]=SetupStst(funcs,...
    'parameter',hopf_br_tk.point(i).parameter,'x',hopf_br_tk.point(i).x,...
    'contpar',17,'corpar',17,'step',0.001);
    branch(k).st.parameter.max_step=[17,0.1];
    branch(k).st.parameter.max_bound=[17 100];
    branch(k).st.parameter.min_bound=[17 0];
    branch(k).st.method.continuation.plot=0;

    [branch(k).st,s,f,r]=br_contn(funcs,branch(k).st,10000);
    branch(k).st=br_rvers(branch(k).st)
    [branch(k).st,s,f,r]=br_contn(funcs,branch(k).st,10000);
    branch(k).st=br_stabl(funcs,branch(k).st,0,0);
    k=k+1;
end
%%
[xmp,ymp]=df_measr(0,branch(1).pr);
zmp=xmp;
zmp.col=14;

[xms,yms]=df_measr(0,branch(1).st);
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
    br_splot_l1(branch(k).pr,xmp,ymp,'k')
    ymp.col='min';
    br_splot_l1(branch(k).pr,xmp,ymp,'k')    
    br_splot_l1(branch(k).st,xms,yms);
    ylim([-35,25])
    xlabel('\tau')
    ylabel('V')
end

%%
figure()
[hbranch,suc]=SetupHopf(funcs,branch(2).st,511,'contpar',[14 17],'correc',1,'dir',[14])
hbranch.parameter.max_step=[14,0.1];
hbranch=br_contn(funcs,hbranch,2000);
hbranch=br_rvers(hbranch);
hbranch=br_contn(funcs,hbranch,2000);

br_plot(hopf_br_tk,xm, ym,'k')
br_plot(hbranch,xm, ym,'r')
%%
figure()
for i=1:1042
    s(i)=p_dststb(branch(2).st.point(i));
end
s2=GetStability(branch(2).st);
[xmp,ymp]=df_measr(1,branch(2).st);
ymp.subfield='l0';
br_plot(branch(2).st,[],ymp)

[xmp,ymp]=df_measr(1,branch(2).st);
ymp.subfield='l1';
br_plot(branch(2).st,[],ymp)
plot(s-2)
plot(s2-1)
plot([0 1042],[0 0],'k')
ylim([-0.1 0.1])
