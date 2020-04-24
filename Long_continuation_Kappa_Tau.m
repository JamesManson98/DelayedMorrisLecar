clear all
close all
load('WORKSPACE','branch2','xm','ym')
myfuncs
%% 2 parameter continuation of hopf for tau and kappa
kappa_600_branch=branch2;
kappa_600_branch.parameter.max_step=[17 0.1; 14 0.1];    
kappa_600_branch.parameter.min_bound=[14 0; 17 0];
kappa_600_branch.parameter.max_bound=[14 100; 17 600];
kappa_600_branch.method.stability.minimal_real_part=-2;
kappa_600_branch=br_rvers(kappa_600_branch);
[kappa_600_branch,s,f,r]=br_contn(funcs,kappa_600_branch,100000);
kappa_600_branch = br_stabl(funcs,kappa_600_branch,0,0);
branch2positivetau=NaN(1,length(kappa_600_branch.point));
branch2negativetau=NaN(1,length(kappa_600_branch.point));
branch2positivekappa=NaN(1,length(kappa_600_branch.point));
branch2negativekappa=NaN(1,length(kappa_600_branch.point));
for i=1:length(kappa_600_branch.point)
    B(i)=nmfm_hopf(funcs, kappa_600_branch.point(i)); % normal form computation for hopf    
    if B(i).nmfm.L1>=0
        branch2positivetau(i)=kappa_600_branch.point(i).parameter(17);
        branch2positivekappa(i)=kappa_600_branch.point(i).parameter(14);
    else
        branch2negativetau(i)=kappa_600_branch.point(i).parameter(17);
        branch2negativekappa(i)=kappa_600_branch.point(i).parameter(14);
    end
end
plot(branch2negativetau(:),branch2negativekappa(:),'color','g')
hold on
plot(branch2positivetau(:),branch2positivekappa(:),'color','r')
xlim([0 300])
xlabel('\tau');ylabel('\kappa');
figure()
plot(branch2negativetau(:),branch2negativekappa(:),'color','g')
hold on
plot(branch2positivetau(:),branch2positivekappa(:),'color','r')
axis([180 280 1.5 3])
xlabel('\tau');ylabel('\kappa');
