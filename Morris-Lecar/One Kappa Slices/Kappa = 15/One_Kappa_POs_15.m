%%Periodic orbit branches emerging from the Hopf points close to kappa = 15
clear all
close all
addpath(genpath('../../../../../dde_biftool_v3.1.1'));
load('Two_par_Kappa_Tau_vars','hopf_br_tk','ind_Tau')
myfuncs
%%
kappa=15; %15, 45, 85
ind=find(diff(sign(arrayfun(@(x) x.parameter(14), hopf_br_tk.point)-kappa))~=0);
figure
hold on
intervals=18;
degree=3;
ind(end)=[];
for i=1:length(ind)
    hopf_br_tk.point(ind(i)).parameter(14)=kappa; %we fix kappa
    method=df_mthod(funcs,'hopf',1); % get hopf calculation method parameters:
    method.stability.minimal_real_part=-1;
    [hopf_br_tk.point(ind(i)),~]=p_correc(funcs,hopf_br_tk.point(ind(i)),17,[],method.point); % correct hopf
    [psol,stepcond]=p_topsol(funcs,hopf_br_tk.point(ind(i)),1e-2,degree,intervals);
    %correct periodic solution guess:
    method=df_mthod(funcs,'psol');
    [psol,success]=p_correc(funcs,psol,17,stepcond,method.point);
    branches_PO_kappa_15(i).branch=df_brnch(funcs,ind_Tau,'psol'); % empty branch:
    % make degenerate periodic solution with amplitude zero at hopf point:
    deg_psol=p_topsol(funcs,hopf_br_tk.point(ind(i)),0,degree,intervals);
    % use deg_psol and psol as first two points on branch:
%     deg_psol.mesh=[];
    branches_PO_kappa_15(i).branch.point=deg_psol;
%     psol.mesh=[];
    branches_PO_kappa_15(i).branch.point(2)=psol;
    branches_PO_kappa_15(i).branch.method.continuation.plot=0;    
    branches_PO_kappa_15(i).branch.parameter.min_bound=[17 -20];
    branches_PO_kappa_15(i).branch.parameter.max_bound=[17 200];
    branches_PO_kappa_15(i).branch.parameter.max_step=[17 0.1; 14 0.1];
    [branches_PO_kappa_15(i).branch,s,f,r]=br_contn(funcs,branches_PO_kappa_15(i).branch,1000);
    br_rvers(branches_PO_kappa_15(i).branch);
    [branches_PO_kappa_15(i).branch,s,f,r]=br_contn(funcs,branches_PO_kappa_15(i).branch,1000);
    %branch5=br_stabl(funcs,branch5,0,0);
    [xm,ym]=df_measr(0,branches_PO_kappa_15(i).branch);
    %we find a low value of amplitude where there is a turning point
    %(allows a lower stepsize)
    same_dim_sol = branches_PO_kappa_15(i).branch.point(1:end-2);
    z=find(diff(sign(diff(arrayfun(@(x) max(x.profile(1,:)), branches_PO_kappa_15(i).branch.point)-arrayfun(@(x) min(x.profile(1,:)), branches_PO_kappa_15(i).branch.point))))==2 & arrayfun(@(x) max(x.profile(1,:)), same_dim_sol)-arrayfun(@(x) min(x.profile(1,:)), same_dim_sol)<5,1,'first');
    if ~isempty(z)
        branches_PO_kappa_15(i).branch.point=branches_PO_kappa_15(i).branch.point(1:z+1);
    end
    xm.col=17;
    ym.col='ampl';
    branches_PO_kappa_15(i).branch = br_stabl(funcs,branches_PO_kappa_15(i).branch,0,0);
    br_splot(branches_PO_kappa_15(i).branch,xm, ym,'b')
end
xlabel('\tau');ylabel('Amplitude');
xlim([-15,200]);