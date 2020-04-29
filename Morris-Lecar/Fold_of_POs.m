%% Continuing the branches of fold of POs
%loading in variables
clear
addpath( genpath('../../../dde_biftool_v3.1.1'));
myfuncs
load('Two_par_Kappa_Tau_vars','hopf_br_tk','xm','ym','GH_points','HH_points')
load('initialpoint_pofold','init_PO_fold','hopf_br_tk_unstable_tau','hopf_br_tk_stable_tau','hopf_br_tk_unstable_kappa','hopf_br_tk_stable_kappa')
%plotting sub/supercritical Hopfs
figure(8)
xlabel('\tau', 'FontSize', 18);ylabel('\kappa', 'FontSize', 18);
hold on
for i=1:length(GH_points)
    plot(hopf_br_tk.point(GH_points(i)).parameter(17),hopf_br_tk.point(GH_points(i)).parameter(14),'Marker','x','color', 'k')
end
for i=1:length(HH_points)
    plot(hopf_br_tk.point(HH_points(i)).parameter(17),hopf_br_tk.point(HH_points(i)).parameter(14),'Marker','o','color', 'k')
end
plot(hopf_br_tk_stable_tau(:),hopf_br_tk_stable_kappa(:),'color','g')
plot(hopf_br_tk_unstable_tau(:),hopf_br_tk_unstable_kappa(:),'color','r')
xlim([0,100]);
ylim([0,100]);
h(1)=plot(1000,1000, 'g');
h(2)=plot(1000,1000, 'r');
h(1)=plot(1000,1000, 'Marker', 'x', 'color', 'k', 'linestyle', 'none');
h(2)=plot(1000,1000, 'Marker', 'o', 'color', 'k', 'linestyle', 'none');
%% initialising, continuting and plotting branches of fold bifurcations 
[foldfuncs,PO_fold(1)]=SetupPOfold(funcs,init_PO_fold(1).br,33,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
PO_fold(1).method.point.print_residual_info=0;
PO_fold(1).method.point.adapt_mesh_after_correct=1;
PO_fold(1)=br_contn(foldfuncs,PO_fold(1),650);
PO_fold(1)=br_rvers(PO_fold);
PO_fold(1)=br_contn(foldfuncs,PO_fold(1),400);
PO_fold(1)=br_rvers(PO_fold);
PO_fold(1).parameter.max_step=[14,0.01; 17,0.01];
PO_fold(1)=br_contn(foldfuncs,PO_fold(1),350);
br_plot(PO_fold(1),xm, ym,'b')
%%
[foldfuncs,PO_fold(2)]=SetupPOfold(funcs,init_PO_fold(2).br,28,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
PO_fold(2).method.point.print_residual_info=0;
PO_fold(2).method.point.adapt_mesh_after_correct=1;
PO_fold(2)=br_contn(foldfuncs,PO_fold(2),370);
PO_fold(2)=br_rvers(PO_fold(2));
PO_fold(2)=br_contn(foldfuncs,PO_fold(2),500);
br_plot(PO_fold(2),xm, ym,'b')
%%
[foldfuncs,PO_fold(3)]=SetupPOfold(funcs,init_PO_fold(3).br,17,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
PO_fold(3).method.point.print_residual_info=0;
PO_fold(3).method.point.adapt_mesh_after_correct=1;
PO_fold(3)=br_contn(foldfuncs,PO_fold(3),50);
PO_fold(3)=br_rvers(PO_fold(3));
PO_fold(3)=br_contn(foldfuncs,PO_fold(3),900);
br_plot(PO_fold(3),xm, ym,'b')
%%
[foldfuncs,PO_fold(4)]=SetupPOfold(funcs,init_PO_fold(4).br,20,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
PO_fold(4).method.point.print_residual_info=0;
PO_fold(4).method.point.adapt_mesh_after_correct=1;
PO_fold(4)=br_contn(foldfuncs,PO_fold(4),50);
PO_fold(4)=br_rvers(PO_fold(4));
PO_fold(4)=br_contn(foldfuncs,PO_fold(4),800);
br_plot(PO_fold(4),xm, ym,'b')
%%
[foldfuncs,PO_fold(5)]=SetupPOfold(funcs,init_PO_fold(5).br,25,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,100],...
    'max_step',[14,0.01; 17,0.01]);
PO_fold(5).method.point.print_residual_info=0;
PO_fold(5).method.point.adapt_mesh_after_correct=1;
PO_fold(5)=br_contn(foldfuncs,PO_fold(5),100);
PO_fold(5)=br_rvers(PO_fold(5));
PO_fold(5)=br_contn(foldfuncs,PO_fold(5),500);
PO_fold(5).parameter.max_step=[14,0.1; 17,0.1];
PO_fold(5)=br_contn(foldfuncs,PO_fold(5),2000); 
br_plot(PO_fold(5),xm, ym,'b')
%%
[foldfuncs,PO_fold(6)]=SetupPOfold(funcs,init_PO_fold(6).br,25,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
PO_fold(6).method.point.print_residual_info=0;
PO_fold(6).method.point.adapt_mesh_after_correct=1;
PO_fold(6)=br_contn(foldfuncs,PO_fold(6),50);
PO_fold(6)=br_rvers(PO_fold(6));
PO_fold(6)=br_contn(foldfuncs,PO_fold(6),800); 
br_plot(PO_fold(6),xm, ym,'b')
%%
[foldfuncs,PO_fold(7)]=SetupPOfold(funcs,init_PO_fold(7).br,27,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
PO_fold(7).method.point.print_residual_info=0;
PO_fold(7).method.point.adapt_mesh_after_correct=1;
PO_fold(7)=br_contn(foldfuncs,PO_fold(7),50);
PO_fold(7)=br_rvers(PO_fold(7));
PO_fold(7)=br_contn(foldfuncs,PO_fold(7),400);
br_plot(PO_fold(7),xm, ym,'b')
%%
[foldfuncs,PO_fold(8)]=SetupPOfold(funcs,init_PO_fold(8).br,16,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
PO_fold(8).method.point.print_residual_info=0;
PO_fold(8).method.point.adapt_mesh_after_correct=1;
PO_fold(8)=br_contn(foldfuncs,PO_fold(8),50);
PO_fold(8)=br_rvers(PO_fold(8));
PO_fold(8)=br_contn(foldfuncs,PO_fold(8),800);
br_plot(PO_fold(8),xm, ym,'b')
%%
[foldfuncs,PO_fold(9)]=SetupPOfold(funcs,init_PO_fold(9).br,17,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
PO_fold(9).method.point.print_residual_info=0;
PO_fold(9).method.point.adapt_mesh_after_correct=1;
PO_fold(9)=br_contn(foldfuncs,PO_fold(9),50);
PO_fold(9)=br_rvers(PO_fold(9));
PO_fold(9)=br_contn(foldfuncs,PO_fold(9),800);
br_plot(PO_fold(9),xm, ym,'b')
%%
[foldfuncs,PO_fold(10)]=SetupPOfold(funcs,init_PO_fold(10).br,15,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
PO_fold(10).method.point.print_residual_info=0;
PO_fold(10).method.point.adapt_mesh_after_correct=1;
PO_fold(10)=br_contn(foldfuncs,PO_fold(10),50);
PO_fold(10)=br_rvers(PO_fold(10));
PO_fold(10)=br_contn(foldfuncs,PO_fold(10),800);
br_plot(PO_fold(10),xm, ym,'b')
%%
[foldfuncs,PO_fold(11)]=SetupPOfold(funcs,init_PO_fold(11).br,14,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,150],...
    'max_step',[14,0.1; 17,0.1]);
PO_fold(11).method.point.print_residual_info=0;
PO_fold(11).method.point.adapt_mesh_after_correct=1;
PO_fold(11)=br_contn(foldfuncs,PO_fold(11),50);
PO_fold(11)=br_rvers(PO_fold(11));
PO_fold(11)=br_contn(foldfuncs,PO_fold(11),800);
br_plot(PO_fold(11),xm, ym,'b')