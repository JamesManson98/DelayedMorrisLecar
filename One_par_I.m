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
parameters_I=[0.125, 4.0, 2.0, 8.0, -2.8, 26.0, 12.0, 17.4, -60.0, 120.0, -91.893652410657, 20.0, -60.0, 0.0, 0.0, 5.0, 10.0];  %coursework last year parameters
% [phi, gca, gl, gk, v1, v2, v3, v4, vl, vca, vk, C, I, Kappa, v_s, Kappa_s, Tau] 

x0=[-100;0]; %initial guess to find an equilibrium

ind_I=13;
[stst_br,suc]=SetupStst(funcs,...
    'parameter',parameters_I,'x',x0,...
    'contpar',ind_I,'step',0.01);
stst_br.parameter.max_step=[13,1];
stst_br.parameter.max_bound=[13 62];
stst_br.parameter.min_bound=[13 -62];
stst_br.method.continuation.plot=0;

[stst_br,s,f,r]=br_contn(funcs,stst_br,400);
stst_br=br_rvers(stst_br)
[stst_br,s,f,r]=br_contn(funcs,stst_br,30);
stst_br=br_stabl(funcs, stst_br,0,0);

[xm,ym]=df_measr(0,stst_br);
br_splot(stst_br,xm,ym)
xlabel('I');ylabel('V');
%this plots hopf at x,y=21.9479,1.7006
%% finding the hopf point:
nunst=GetStability(stst_br,'funcs',funcs);
indhopf=find(abs(diff(nunst))==2);
hopf=p_tohopf(funcs,stst_br.point(indhopf)); 
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
[psol,success]=p_correc(funcs,psol,4,stepcond,method.point);
branch4=df_brnch(funcs,ind_I,'psol'); % empty branch:
% make degenerate periodic solution with amplitude zero at hopf point:
deg_psol=p_topsol(funcs,first_hopf,0,degree,intervals);
% use deg_psol and psol as first two points on branch:
% deg_psol.mesh=[];
branch4.point=deg_psol;
% psol.mesh=[];
branch4.point(2)=psol;
figure(2)
branch4.parameter.max_step=[13,1];
[branch4,s,f,r]=br_contn(funcs,branch4,59); 
xlabel('I');ylabel('Amplitude');
%%
figure(1)
hold on
[xm,ym]=df_measr(0,branch4);
ym.col='max';
br_plot(branch4,xm,ym,'k')
ym.col='min';
br_plot(branch4,xm,ym,'k')
axis([-60,60,-100,20])
names={'Branch of periodic orbits','Fold bifurcation','Hopf bifurcation'};
legend(h,names,'location','SouthEast','FontSize',18)
figure(3)
subplot(2,2,1)
for i=1:40
    hold on 
    if i~=33
        plot(branch4.point(i).profile(1,:),branch4.point(i).profile(2,:),'k'),
    else
        plot(branch4.point(i).profile(1,:),branch4.point(i).profile(2,:),'r'),
    end
end
xlabel('V');ylabel('W');title('(a)');
time=branch4.point(33).mesh*branch4.point(33).period;
time=[time,branch4.point(33).period+time,2*branch4.point(33).period+time];
set(gca,'FontSize',18)
subplot(2,2,4)
plot(time,[branch4.point(33).profile(2,:),branch4.point(33).profile(2,:),branch4.point(33).profile(2,:)],'r');
xlim([0,75]);
xlabel('Time (ms)');ylabel('W');title('(d)');
set(gca,'FontSize',18)
subplot(2,2,3)
plot(time,[branch4.point(33).profile(1,:),branch4.point(33).profile(1,:),branch4.point(33).profile(1,:)],'r');
xlim([0,75]);
xlabel('Time (ms)');ylabel('V');title('(c)');
set(gca,'FontSize',18)
subplot(2,2,2)
ym.col=1;
ym.field='period';
br_plot(branch4,xm,ym,'k')
ylim([20,120])
xlabel('I');ylabel('Period');title('(b)');
set(gca,'FontSize',18)