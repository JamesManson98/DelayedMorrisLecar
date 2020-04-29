%% code to produce 1 parameter continuations and bifurcations in I
clear all
close all
addpath( genpath('../../../dde_biftool_v3.1.1'));
myfuncs
%%
figure(1)
hold on
h(1)=plot(linspace(-200,-199,2),linspace(-200,-199,2),'k');
h(2)=plot(-200,-200,'kx');
h(3)=plot(-200,-200,'ko','MarkerSize', 6,'MarkerFaceColor', 'k');
parameters_I=[0.125, 4.0, 2.0, 8.0, -2.8, 26.0, 12.0, 17.4, -60.0, 120.0, -91.893652410657, 20.0, -60.0, 0.0, 0.0, 5.0, 10.0];
% [phi, gca, gl, gk, v1, v2, v3, v4, vl, vca, vk, C, I, Kappa, v_s, Kappa_s, Tau] 

x0=[-100;0]; %initial guess to find an equilibrium

ind_I=13;
ind_Kappa=14;
ind_Tau=17;
[stst_br_I,suc]=SetupStst(funcs,...
    'parameter',parameters_I,'x',x0,...
    'contpar',ind_I,'step',0.01);
% stst_br_I.parameter.max_step=[13,1]; for more accurate computation set
% a small max step and a longer continuation 
% stst_br_I.parameter.max_bound=[13 62];
% stst_br_I.parameter.min_bound=[13 -62]; 
stst_br_I.method.continuation.plot=0;

[stst_br_I,s,f,r]=br_contn(funcs,stst_br_I,90);
stst_br_I=br_rvers(stst_br_I);
[stst_br_I,s,f,r]=br_contn(funcs,stst_br_I,30);
stst_br_I=br_stabl(funcs, stst_br_I,0,0);

[xm,ym]=df_measr(0,stst_br_I);
br_splot(stst_br_I,xm,ym)
xlabel('I');ylabel('V');
%this plots hopf at x,y=21.9479,1.7006
%% finding the hopf point:
nunst=GetStability(stst_br_I,'funcs',funcs);
indhopf=find(abs(diff(nunst))==2);
hopf=p_tohopf(funcs,stst_br_I.point(indhopf)); 
method=df_mthod(funcs,'hopf',1); % get hopf calculation method parameters:
method.stability.minimal_real_part=-1;
[hopf,success]=p_correc(funcs,hopf,ind_I,[],method.point); % correct hopf
first_hopf=hopf;                    % store hopf point in other variable for later use
hopf.stability=p_stabil(funcs,hopf,method.stability); % compute stability of hopf point
%% periodic solutions:
intervals=18;
degree=3;
[psol,stepcond]=p_topsol(funcs,first_hopf,1e-2,degree,intervals);
%correct periodic solution guess:
method=df_mthod(funcs,'psol');
[psol,success]=p_correc(funcs,psol,ind_I,stepcond,method.point);
PO_branch_I=df_brnch(funcs,ind_I,'psol'); % empty branch:
% make degenerate periodic solution with amplitude zero at hopf point:
deg_psol=p_topsol(funcs,first_hopf,0,degree,intervals);
% use deg_psol and psol as first two points on branch:
% deg_psol.mesh=[];
PO_branch_I.point=deg_psol;
% psol.mesh=[];
PO_branch_I.point(2)=psol;
figure(2)
PO_branch_I.parameter.max_step=[13,1];
[PO_branch_I,s,f,r]=br_contn(funcs,PO_branch_I,59); 
xlabel('I');ylabel('Amplitude');
%%
figure(1)
hold on
[xm,ym]=df_measr(0,PO_branch_I);
ym.col='max';
br_plot(PO_branch_I,xm,ym,'k')
ym.col='min';
br_plot(PO_branch_I,xm,ym,'k')
axis([-60,60,-100,20])
names={'Branch of periodic orbits','Fold bifurcation','Hopf bifurcation'};
legend(h,names,'location','SouthEast','FontSize',18)
figure(3)
subplot(2,2,1)
for i=1:40
    hold on 
    if i~=33
        plot(PO_branch_I.point(i).profile(1,:),PO_branch_I.point(i).profile(2,:),'k'),
    else
        plot(PO_branch_I.point(i).profile(1,:),PO_branch_I.point(i).profile(2,:),'r'),
    end
end
xlabel('V');ylabel('W');title('(a)');
time=PO_branch_I.point(33).mesh*PO_branch_I.point(33).period;
time=[time,PO_branch_I.point(33).period+time,2*PO_branch_I.point(33).period+time];
set(gca,'FontSize',18)
subplot(2,2,4)
plot(time,[PO_branch_I.point(33).profile(2,:),PO_branch_I.point(33).profile(2,:),PO_branch_I.point(33).profile(2,:)],'r');
xlim([0,75]);
xlabel('Time (ms)');ylabel('W');title('(d)');
set(gca,'FontSize',18)
subplot(2,2,3)
plot(time,[PO_branch_I.point(33).profile(1,:),PO_branch_I.point(33).profile(1,:),PO_branch_I.point(33).profile(1,:)],'r');
xlim([0,75]);
xlabel('Time (ms)');ylabel('V');title('(c)');
set(gca,'FontSize',18)
subplot(2,2,2)
ym.col=1;
ym.field='period';
br_plot(PO_branch_I,xm,ym,'k')
ylim([20,120])
xlabel('I');ylabel('Period');title('(b)');
set(gca,'FontSize',18)