%% Torus bifurcations
%clear
% close all
addpath( genpath('../../../dde_biftool_v3.1.1'));
myfuncs
load('Two_par_Kappa_Tau_vars','branch2','xm','ym','GH_points','HH_points')
load('initialpoint_torus','torus_branch','branch2negativekappa','branch2negativetau','branch2positivekappa','branch2positivetau')
% load('TEST','BRANCH_TORUS')
% load('TEST_2','BRANCH_TORUS_2')
% load('TEST_3','BRANCH_TORUS_3')
% load('Random_test','TEST_BRANCH_5')
% load('Random_test_2','test_5_branch')
% load('Random_test_3','more_test_br','test_branch','test_branch_11')
% load('test_br_4','test_branch_4')
% load('test_br_19','test_branch_19')
% load('test_br_29','test_branch_29')
% load('test_br_13','test_branch_13')
% load('test_br_15','test_branch_15')
% load('test_br_16','test_branch_16')
% load('test_br_18','test_branch_18')
% load('test_br_27','test_branch_27')
% load('test_br_30','test_branch_30')
% load('test_br_1','HH_points1_minus20_branch')
load('Torus_branches','HH_points1_minus20_branch','test_branch_4','test_branch_27','test_branch_30','test_branch_18','test_branch_19','test_branch_29','test_branch_13','test_branch_15','test_branch_16','BRANCH_TORUS','torus_branch','BRANCH_TORUS_2','BRANCH_TORUS_3','TEST_BRANCH_5','test_5_branch','more_test_br','test_branch','test_branch_11')

xlabel('\tau', 'FontSize', 18);ylabel('\kappa', 'FontSize', 18);
hold on
plot(branch2negativetau(:),branch2negativekappa(:),'color','g')
plot(branch2positivetau(:),branch2positivekappa(:),'color','r')
for i=1:length(GH_points)
    plot(branch2.point(GH_points(i)).parameter(17),branch2.point(GH_points(i)).parameter(14),'Marker','x','color', 'k')
end
for i=1:length(HH_points)
    plot(branch2.point(HH_points(i)).parameter(17),branch2.point(HH_points(i)).parameter(14),'Marker','o','color', 'k')
end
xlim([0,100]);
ylim([0,100]);
% h(1)=plot(1000,1000, 'g');
% h(2)=plot(1000,1000, 'r');
% h(3)=plot(1000,1000, 'b');
% h(4)=plot(1000,1000, 'k');
% h(5)=plot(1000,1000, 'Marker', 'x', 'color', 'k', 'linestyle', 'none');
% h(6)=plot(1000,1000, 'Marker', 'o', 'color', 'k', 'linestyle', 'none');
xm.field='parameter';
xm.col=17;
ym.field='parameter';
ym.col=14;
% legend(h, {'Stable Hopf continuation','Unstable Hopf continuation','Fold of POs', 'Torus Bifurcations', 'Generalised Hopf Points', 'HH points'},'FontSize', 18)
%%
[foldfuncs,branchtorus(1).br]=SetupTorusBifurcation(funcs, torus_branch(1).br,18,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus(1).br.method.point.print_residual_info=0;
branchtorus(1).br.method.point.adapt_mesh_after_correct=1;
branchtorus(1).br=br_contn(foldfuncs,branchtorus(1).br,50);
branchtorus(1).br=br_rvers(branchtorus(1).br);
branchtorus(1).br=br_contn(foldfuncs,branchtorus(1).br,4000);
br_plot(branchtorus(1).br,xm, ym,'k')
%%
[foldfuncs,branchtorus_2]=SetupTorusBifurcation(funcs, BRANCH_TORUS_2,105,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus_2.method.point.print_residual_info=0;
branchtorus_2.method.point.adapt_mesh_after_correct=1;
branchtorus_2=br_contn(foldfuncs,branchtorus_2,1000);
branchtorus_2=br_rvers(branchtorus_2);
branchtorus_2=br_contn(foldfuncs,branchtorus_2,1000);
br_plot(branchtorus_2,xm, ym,'k')
%% KEEP
[foldfuncs,branchtorus_3]=SetupTorusBifurcation(funcs, torus_branch(3).br,24,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus_3.method.point.print_residual_info=0;
branchtorus_3.method.point.adapt_mesh_after_correct=1;
branchtorus_3=br_contn(foldfuncs,branchtorus_3,50);
branchtorus_3=br_rvers(branchtorus_3);
branchtorus_3=br_contn(foldfuncs,branchtorus_3,2000);
br_plot(branchtorus_3,xm, ym,'m')
%%
[foldfuncs,branchtorus_4]=SetupTorusBifurcation(funcs, BRANCH_TORUS,23,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus_4.method.point.print_residual_info=0;
branchtorus_4.method.point.adapt_mesh_after_correct=1;
branchtorus_4=br_contn(foldfuncs,branchtorus_4,50);
branchtorus_4=br_rvers(branchtorus_4);
branchtorus_4=br_contn(foldfuncs,branchtorus_4,1000);
br_plot(branchtorus_4,xm, ym,'k')

%%
%19 20 91 193 194
for i=[91,193]
    [foldfuncs,branchtorus_5]=SetupTorusBifurcation(funcs, BRANCH_TORUS_3,i,'contpar',[14,17],...
        'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
        'max_step',[14,0.1; 17,0.1]);
    branchtorus_5.method.point.print_residual_info=0;
    branchtorus_5.method.point.adapt_mesh_after_correct=1;
    branchtorus_5=br_contn(foldfuncs,branchtorus_5,2000);
    branchtorus_5=br_rvers(branchtorus_5);
    branchtorus_5=br_contn(foldfuncs,branchtorus_5,200);
    br_plot(branchtorus_5,xm, ym,'k')
    if i==91
        branchtorus_5_1=branchtorus_5;
    end
end
%%
[foldfuncs,branchtorus_6]=SetupTorusBifurcation(funcs, TEST_BRANCH_5,90,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus_6.method.point.print_residual_info=0;
branchtorus_6.method.point.adapt_mesh_after_correct=1;
branchtorus_6=br_contn(foldfuncs,branchtorus_6,1000);
branchtorus_6=br_rvers(branchtorus_6);
branchtorus_6=br_contn(foldfuncs,branchtorus_6,1000);
br_plot(branchtorus_6,xm, ym,'k')
%%
[foldfuncs,branchtorus_7]=SetupTorusBifurcation(funcs, test_5_branch,121,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus_7.method.point.print_residual_info=0;
branchtorus_7.method.point.adapt_mesh_after_correct=1;
branchtorus_7=br_contn(foldfuncs,branchtorus_7,600);
branchtorus_7=br_rvers(branchtorus_7);
branchtorus_7=br_contn(foldfuncs,branchtorus_7,900);
br_plot(branchtorus_7,xm, ym,'k')
%%
[foldfuncs,branchtorus_8]=SetupTorusBifurcation(funcs, more_test_br,135,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus_8.method.point.print_residual_info=0;
branchtorus_8.method.point.adapt_mesh_after_correct=1;
branchtorus_8=br_contn(foldfuncs,branchtorus_8,1000);
branchtorus_8=br_rvers(branchtorus_8);
branchtorus_8=br_contn(foldfuncs,branchtorus_8,1000);
br_plot(branchtorus_8,xm, ym,'k')
%% KEEP
[foldfuncs,branchtorus_9]=SetupTorusBifurcation(funcs, test_branch,112,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,150],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus_9.method.point.print_residual_info=0;
branchtorus_9.method.point.adapt_mesh_after_correct=1;
branchtorus_9=br_contn(foldfuncs,branchtorus_9,200);
branchtorus_9=br_rvers(branchtorus_9);
branchtorus_9=br_contn(foldfuncs,branchtorus_9,200); %increase this
br_plot(branchtorus_9,xm, ym,'m')
%% 17 165 84/85
for i=[17,84,165]
    [foldfuncs,branchtorus_10]=SetupTorusBifurcation(funcs, test_branch_11,i,'contpar',[14,17],...
        'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
        'max_step',[14,0.1; 17,0.1]);
    branchtorus_10.method.point.print_residual_info=0;
    branchtorus_10.method.point.adapt_mesh_after_correct=1;
    branchtorus_10=br_contn(foldfuncs,branchtorus_10,400);
    branchtorus_10=br_rvers(branchtorus_10);
    branchtorus_10=br_contn(foldfuncs,branchtorus_10,800);
    br_plot(branchtorus_10,xm, ym,'k')
    if i==17
        branchtorus_10_1=branchtorus_10;
    elseif i==84
        branchtorus_10_2=branchtorus_10;
    else
        branchtorus_10_3=branchtorus_10;
    end
end
%%
[foldfuncs,branchtorus_11]=SetupTorusBifurcation(funcs, test_branch_4,10,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus_11.method.point.print_residual_info=0;
branchtorus_11.method.point.adapt_mesh_after_correct=1;
branchtorus_11=br_contn(foldfuncs,branchtorus_11,200);
branchtorus_11=br_rvers(branchtorus_11);
branchtorus_11=br_contn(foldfuncs,branchtorus_11,1000);
br_plot(branchtorus_11,xm, ym,'k')
%%
[foldfuncs,branchtorus_12]=SetupTorusBifurcation(funcs, test_branch_19,15,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus_12.method.point.print_residual_info=0;
branchtorus_12.method.point.adapt_mesh_after_correct=1;
branchtorus_12=br_contn(foldfuncs,branchtorus_12,200);
branchtorus_12=br_rvers(branchtorus_12);
branchtorus_12=br_contn(foldfuncs,branchtorus_12,800);
br_plot(branchtorus_12,xm, ym,'k')
%%
[foldfuncs,branchtorus_13]=SetupTorusBifurcation(funcs, test_branch_29,31,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,150],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus_13.method.point.print_residual_info=0;
branchtorus_13.method.point.adapt_mesh_after_correct=1;
branchtorus_13=br_contn(foldfuncs,branchtorus_13,200);
branchtorus_13=br_rvers(branchtorus_13);
branchtorus_13=br_contn(foldfuncs,branchtorus_13,800);
br_plot(branchtorus_13,xm, ym,'k')
%%
[foldfuncs,branchtorus_14]=SetupTorusBifurcation(funcs, test_branch_13,86,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,150],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus_14.method.point.print_residual_info=0;
branchtorus_14.method.point.adapt_mesh_after_correct=1;
branchtorus_14=br_contn(foldfuncs,branchtorus_14,200);
branchtorus_14=br_rvers(branchtorus_14);
branchtorus_14=br_contn(foldfuncs,branchtorus_14,200);
br_plot(branchtorus_14,xm, ym,'k')
%%
[foldfuncs,branchtorus_15]=SetupTorusBifurcation(funcs, test_branch_15,61,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,150],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus_15.method.point.print_residual_info=0;
branchtorus_15.method.point.adapt_mesh_after_correct=1;
branchtorus_15=br_contn(foldfuncs,branchtorus_15,400);
branchtorus_15=br_rvers(branchtorus_15);
branchtorus_15=br_contn(foldfuncs,branchtorus_15,400);
br_plot(branchtorus_15,xm, ym,'m')
%%
[foldfuncs,branchtorus_16]=SetupTorusBifurcation(funcs, test_branch_16,59,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,150],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus_16.method.point.print_residual_info=0;
branchtorus_16.method.point.adapt_mesh_after_correct=1;
branchtorus_16=br_contn(foldfuncs,branchtorus_16,1600);
branchtorus_16=br_rvers(branchtorus_16);
branchtorus_16=br_contn(foldfuncs,branchtorus_16,200);
br_plot(branchtorus_16,xm, ym,'k')
%%
[foldfuncs,branchtorus_17]=SetupTorusBifurcation(funcs, test_branch_18,82,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,150],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus_17.method.point.print_residual_info=0;
branchtorus_17.method.point.adapt_mesh_after_correct=1;
branchtorus_17=br_contn(foldfuncs,branchtorus_17,400);
branchtorus_17=br_rvers(branchtorus_17);
branchtorus_17=br_contn(foldfuncs,branchtorus_17,1200);
br_plot(branchtorus_17,xm, ym,'k')
%%
[foldfuncs,branchtorus_18]=SetupTorusBifurcation(funcs, test_branch_27,81,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,150],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus_18.method.point.print_residual_info=0;
branchtorus_18.method.point.adapt_mesh_after_correct=1;
branchtorus_18=br_contn(foldfuncs,branchtorus_18,1000);
branchtorus_18=br_rvers(branchtorus_18);
branchtorus_18=br_contn(foldfuncs,branchtorus_18,2000);
br_plot(branchtorus_18,xm, ym,'k')
%%
[foldfuncs,branchtorus_19]=SetupTorusBifurcation(funcs, HH_points1_minus20_branch,1549,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,150],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus_19.method.point.print_residual_info=0;
branchtorus_19.method.point.adapt_mesh_after_correct=1;
branchtorus_19=br_contn(foldfuncs,branchtorus_19,3000);
branchtorus_19=br_rvers(branchtorus_19);
branchtorus_19=br_contn(foldfuncs,branchtorus_19,400);
br_plot(branchtorus_19,xm, ym,'m')%?
%%
[foldfuncs,branchtorus_20]=SetupTorusBifurcation(funcs, test_branch_30,29,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,150],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus_20.method.point.print_residual_info=0;
branchtorus_20.method.point.adapt_mesh_after_correct=1;
branchtorus_20=br_contn(foldfuncs,branchtorus_20,200);
branchtorus_20=br_rvers(branchtorus_20);
branchtorus_20=br_contn(foldfuncs,branchtorus_20,1200);
br_plot(branchtorus_20,xm, ym,'k')
%% ???
% load('One_Kappa_POs_45','branches_PO_kappa_45')
% [foldfuncs,branchtorus_21]=SetupTorusBifurcation(funcs, branches_PO_kappa_45(2).branch,264,'contpar',[14,17],...
%     'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,150],...
%     'max_step',[14,0.1; 17,0.1]);
% branchtorus_21.method.point.print_residual_info=0;
% branchtorus_21.method.point.adapt_mesh_after_correct=1;
% branchtorus_21=br_contn(foldfuncs,branchtorus_21,200);
% branchtorus_21=br_rvers(branchtorus_21);
% branchtorus_21=br_contn(foldfuncs,branchtorus_21,400);
% br_plot(branchtorus_21,xm, ym,'m')
%% ???
% load('One_Kappa_POs_85','branches_PO_kappa_85')
% [foldfuncs,branchtorus_22]=SetupTorusBifurcation(funcs, branches_PO_kappa_85(9).branch,583,'contpar',[14,17],...
%     'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,150],...
%     'max_step',[14,0.01; 17,0.01]);
% branchtorus_22.method.point.print_residual_info=0;
% branchtorus_22.method.point.adapt_mesh_after_correct=1;
% branchtorus_22=br_contn(foldfuncs,branchtorus_22,200);
% branchtorus_22=br_rvers(branchtorus_22);
% branchtorus_22=br_contn(foldfuncs,branchtorus_22,800);
% br_plot(branchtorus_22,xm, ym,'m')
%% same as one before but joins up better
load('branches_PO_kappa_95','branches_PO_kappa_85')
[foldfuncs,branchtorus_25]=SetupTorusBifurcation(funcs, branches_PO_kappa_85(7).branch,621,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,150],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus_25.method.point.print_residual_info=0;
branchtorus_25.method.point.adapt_mesh_after_correct=1;
branchtorus_25=br_contn(foldfuncs,branchtorus_25,200);
branchtorus_25=br_rvers(branchtorus_25);
branchtorus_25=br_contn(foldfuncs,branchtorus_25,1200);
br_plot(branchtorus_25,xm, ym,'m')
