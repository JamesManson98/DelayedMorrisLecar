%%checking what happens at a swallowtail
close all
clear all
addpath( genpath('../../../dde_biftool_v3.1.1'));
myfuncs
load('Fold_of_POs','PO_fold','ym','xm')
inds=[100,300,650];
for i=1:length(inds)
    point_near_swallowtail=PO_fold(5).point(inds(i));
    point_near_swallowtail.parameter=point_near_swallowtail.parameter(1:17);
    point_near_swallowtail.profile=point_near_swallowtail.profile(1:2,:);

    point_near_swallowtail_2=point_near_swallowtail;
    point_near_swallowtail_2.parameter(17)=point_near_swallowtail.parameter(17)+0.001;
    branch=df_brnch(funcs,17,'psol');
    branch.point=point_near_swallowtail;
    branch.point(2)=point_near_swallowtail_2;
    method=df_mthod(funcs,'psol',1); % get hopf calculation method parameters:
    %method.stability.minimal_real_part=-1;
    [branch.point(2),~]=p_correc(funcs,branch.point(2),17,[],method.point); % correct hopf

    branch.method.continuation.plot=0;
    [branch,s,f,r]=br_contn(funcs,branch,100);
    sw(i).br=branch;
    [xm,ym]=df_measr(0, branch);
    zm=ym;
    ym.field='parameter';
    ym.col=14;
    br_plot3(branch,xm,ym,zm,'k')
    xlabel('\tau');ylabel('\kappa');zlabel('V');
end
br_plot3(PO_fold(5),xm, ym, zm, 'b') %HOW TO PLOT IN 3D A 2D BRANCH??
