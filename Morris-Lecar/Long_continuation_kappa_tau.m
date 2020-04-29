clear all
close all
load('Two_par_Kappa_Tau_vars','hopf_br_tk','xm','ym')
myfuncs
%% 2 parameter continuation of hopf for tau and kappa
kappa_600_branch=hopf_br_tk;
kappa_600_branch.parameter.max_step=[17 0.1; 14 0.1];    
kappa_600_branch.parameter.min_bound=[14 0; 17 0];
kappa_600_branch.parameter.max_bound=[14 100; 17 600];
kappa_600_branch.method.stability.minimal_real_part=-2;
kappa_600_branch=br_rvers(kappa_600_branch);
[kappa_600_branch,s,f,r]=br_contn(funcs,kappa_600_branch,100000);
kappa_600_branch = br_stabl(funcs,kappa_600_branch,0,0);
hopf_br_tk_unstable_tau=NaN(1,length(kappa_600_branch.point));
hopf_br_tk_stable_tau=NaN(1,length(kappa_600_branch.point));
hopf_br_tk_unstable_kappa=NaN(1,length(kappa_600_branch.point));
hopf_br_tk_stable_kappa=NaN(1,length(kappa_600_branch.point));
for i=1:length(kappa_600_branch.point)
    B(i)=nmfm_hopf(funcs, kappa_600_branch.point(i)); % normal form computation for hopf    
    if B(i).nmfm.L1>=0
        hopf_br_tk_unstable_tau(i)=kappa_600_branch.point(i).parameter(17);
        hopf_br_tk_unstable_kappa(i)=kappa_600_branch.point(i).parameter(14);
    else
        hopf_br_tk_stable_tau(i)=kappa_600_branch.point(i).parameter(17);
        hopf_br_tk_stable_kappa(i)=kappa_600_branch.point(i).parameter(14);
    end
end
plot(hopf_br_tk_stable_tau(:),hopf_br_tk_stable_kappa(:),'color','g')
hold on
plot(hopf_br_tk_unstable_tau(:),hopf_br_tk_unstable_kappa(:),'color','r')
xlim([0 300])
xlabel('\tau');ylabel('\kappa');
figure()
plot(hopf_br_tk_stable_tau(:),hopf_br_tk_stable_kappa(:),'color','g')
hold on
plot(hopf_br_tk_unstable_tau(:),hopf_br_tk_unstable_kappa(:),'color','r')
axis([180 280 1.5 3])
xlabel('\tau');ylabel('\kappa');
