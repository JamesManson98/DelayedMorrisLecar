clear all
close all
addpath( genpath('../../../dde_biftool_v3.1.1'));
myfuncs
load('Two_par_Kappa_Tau_vars','xm','ym','GH_points','HH_points')
load('continuation_kappa_600','kappa_600_branch')
load('Two_par_Kappa_Tau_vars','branchesI_tk','indhopf')
ind_Tau=17;
figure()
xlabel('\tau', 'FontSize', 18);ylabel('\kappa', 'FontSize', 18);
hold on
% [kappa_600_branch,suc]=SetupHopf(funcs,branchesI_tk(1).br,indhopf(1),'contpar',[14,17],...
%     'dir',14,'step',0.001,'plot',0);
% kappa_600_branch = br_stabl(funcs,kappa_600_branch,0,0);
branch2positivetau=NaN(1,length(kappa_600_branch.point));
branch2negativetau=NaN(1,length(kappa_600_branch.point));
branch2positivekappa=NaN(1,length(kappa_600_branch.point));
branch2negativekappa=NaN(1,length(kappa_600_branch.point));
for i=1:length(kappa_600_branch.point)
    A(i)=nmfm_hopf(funcs, kappa_600_branch.point(i)); % normal form computation for hopf    
    if A(i).nmfm.L1>=0
        branch2positivetau(i)=kappa_600_branch.point(i).parameter(17);
        branch2positivekappa(i)=kappa_600_branch.point(i).parameter(14);
    else
        branch2negativetau(i)=kappa_600_branch.point(i).parameter(17);
        branch2negativekappa(i)=kappa_600_branch.point(i).parameter(14);
    end
end
for i=1:2
    subplot(2,1,i)
    plot(branch2negativetau(:),branch2negativekappa(:),'color','g')
    hold on
    plot(branch2positivetau(:),branch2positivekappa(:),'color','r')
    xlabel('\tau', 'FontSize', 18);ylabel('\kappa', 'FontSize', 18);
end