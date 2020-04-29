%% Torus bifurcations
clear
addpath( genpath('../../../dde_biftool_v3.1.1'));
myfuncs
load('Two_par_Kappa_Tau_vars','hopf_br_tk','xm','ym','GH_points','HH_points')
load('initialpoint_torus','init_torus','hopf_br_tk_unstable_tau','hopf_br_tk_stable_tau','hopf_br_tk_unstable_kappa','hopf_br_tk_stable_kappa')

%plotting sub/supercritical Hopfs
figure(9)
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
xm.field='parameter';
xm.col=17;
ym.field='parameter';
ym.col=14;
%% initialising, continuing and plotting branches of Torus bifurcations
[foldfuncs,branchtorus(1)]=SetupTorusBifurcation(funcs, init_torus(1),18,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus(1).method.point.print_residual_info=0;
branchtorus(1).method.point.adapt_mesh_after_correct=1;
branchtorus(1)=br_contn(foldfuncs,branchtorus(1),50);
branchtorus(1)=br_rvers(branchtorus(1));
branchtorus(1)=br_contn(foldfuncs,branchtorus(1),4000);
br_plot(branchtorus(1),xm, ym,'k')
%%
[foldfuncs,branchtorus(2)]=SetupTorusBifurcation(funcs, init_torus(2),105,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus(2).method.point.print_residual_info=0;
branchtorus(2).method.point.adapt_mesh_after_correct=1;
branchtorus(2)=br_contn(foldfuncs,branchtorus(2),1000);
branchtorus(2)=br_rvers(branchtorus(2));
branchtorus(2)=br_contn(foldfuncs,branchtorus(2),1000);
br_plot(branchtorus(2),xm, ym,'k')
%% KEEP
[foldfuncs,branchtorus(3)]=SetupTorusBifurcation(funcs, init_torus(3),24,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus(3).method.point.print_residual_info=0;
branchtorus(3).method.point.adapt_mesh_after_correct=1;
branchtorus(3)=br_contn(foldfuncs,branchtorus(3),50);
branchtorus(3)=br_rvers(branchtorus(3));
branchtorus(3)=br_contn(foldfuncs,branchtorus(3),2000);
br_plot(branchtorus(3),xm, ym,'m')
%%
[foldfuncs,branchtorus(4)]=SetupTorusBifurcation(funcs, init_torus(4),23,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus(4).method.point.print_residual_info=0;
branchtorus(4).method.point.adapt_mesh_after_correct=1;
branchtorus(4)=br_contn(foldfuncs,branchtorus(4),50);
branchtorus(4)=br_rvers(branchtorus(4));
branchtorus(4)=br_contn(foldfuncs,branchtorus(4),1000);
br_plot(branchtorus(4),xm, ym,'k')

%%
%19 20 91 193 194
for i=[91,193]
    [foldfuncs,branchtorus(5)]=SetupTorusBifurcation(funcs, init_torus(5),i,'contpar',[14,17],...
        'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
        'max_step',[14,0.1; 17,0.1]);
    branchtorus(5).method.point.print_residual_info=0;
    branchtorus(5).method.point.adapt_mesh_after_correct=1;
    branchtorus(5)=br_contn(foldfuncs,branchtorus(5),2000);
    branchtorus(5)=br_rvers(branchtorus(5));
    branchtorus(5)=br_contn(foldfuncs,branchtorus(5),200);
    br_plot(branchtorus(5),xm, ym,'k')
    if i==91
        branchtorus(22)=branchtorus(5);
    end
end
%%
[foldfuncs,branchtorus(6)]=SetupTorusBifurcation(funcs, init_torus(6),90,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus(6).method.point.print_residual_info=0;
branchtorus(6).method.point.adapt_mesh_after_correct=1;
branchtorus(6)=br_contn(foldfuncs,branchtorus(6),1000);
branchtorus(6)=br_rvers(branchtorus(6));
branchtorus(6)=br_contn(foldfuncs,branchtorus(6),1000);
br_plot(branchtorus(6),xm, ym,'k')
%%
[foldfuncs,branchtorus(7)]=SetupTorusBifurcation(funcs, init_torus(7),121,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus(7).method.point.print_residual_info=0;
branchtorus(7).method.point.adapt_mesh_after_correct=1;
branchtorus(7)=br_contn(foldfuncs,branchtorus(7),600);
branchtorus(7)=br_rvers(branchtorus(7));
branchtorus(7)=br_contn(foldfuncs,branchtorus(7),900);
br_plot(branchtorus(7),xm, ym,'k')
%%
[foldfuncs,branchtorus(8)]=SetupTorusBifurcation(funcs, init_torus(8),135,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus(8).method.point.print_residual_info=0;
branchtorus(8).method.point.adapt_mesh_after_correct=1;
branchtorus(8)=br_contn(foldfuncs,branchtorus(8),1000);
branchtorus(8)=br_rvers(branchtorus(8));
branchtorus(8)=br_contn(foldfuncs,branchtorus(8),1000);
br_plot(branchtorus(8),xm, ym,'k')
%% KEEP
[foldfuncs,branchtorus(9)]=SetupTorusBifurcation(funcs, init_torus(9),112,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,150],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus(9).method.point.print_residual_info=0;
branchtorus(9).method.point.adapt_mesh_after_correct=1;
branchtorus(9)=br_contn(foldfuncs,branchtorus(9),200);
branchtorus(9)=br_rvers(branchtorus(9));
branchtorus(9)=br_contn(foldfuncs,branchtorus(9),200); %increase this
br_plot(branchtorus(9),xm, ym,'m')
%% 17 165 84/85
for i=[17,84,165]
    [foldfuncs,branchtorus(10)]=SetupTorusBifurcation(funcs, init_torus(10),i,'contpar',[14,17],...
        'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
        'max_step',[14,0.1; 17,0.1]);
    branchtorus(10).method.point.print_residual_info=0;
    branchtorus(10).method.point.adapt_mesh_after_correct=1;
    branchtorus(10)=br_contn(foldfuncs,branchtorus(10),400);
    branchtorus(10)=br_rvers(branchtorus(10));
    branchtorus(10)=br_contn(foldfuncs,branchtorus(10),800);
    br_plot(branchtorus(10),xm, ym,'k')
    if i==17
        branchtorus(23)=branchtorus(10);
    elseif i==84
        branchtorus(24)=branchtorus(10);
    end
end
%%
[foldfuncs,branchtorus(11)]=SetupTorusBifurcation(funcs, init_torus(11),10,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus(11).method.point.print_residual_info=0;
branchtorus(11).method.point.adapt_mesh_after_correct=1;
branchtorus(11)=br_contn(foldfuncs,branchtorus(11),200);
branchtorus(11)=br_rvers(branchtorus(11));
branchtorus(11)=br_contn(foldfuncs,branchtorus(11),1000);
br_plot(branchtorus(11),xm, ym,'k')
%%
[foldfuncs,branchtorus(12)]=SetupTorusBifurcation(funcs, init_torus(12),81,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,150],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus(12).method.point.print_residual_info=0;
branchtorus(12).method.point.adapt_mesh_after_correct=1;
branchtorus(12)=br_contn(foldfuncs,branchtorus(12),1000);
branchtorus(12)=br_rvers(branchtorus(12));
branchtorus(12)=br_contn(foldfuncs,branchtorus(12),2000);
br_plot(branchtorus(12),xm, ym,'k')
%%
[foldfuncs,branchtorus(13)]=SetupTorusBifurcation(funcs, init_torus(13),46,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus(13).method.point.print_residual_info=0;
branchtorus(13).method.point.adapt_mesh_after_correct=1;
branchtorus(13)=br_contn(foldfuncs,branchtorus(13),3000);
branchtorus(13)=br_rvers(branchtorus(13));
branchtorus(13)=br_contn(foldfuncs,branchtorus(13),2000);
br_plot(branchtorus(13),xm, ym,'m')
%%
[foldfuncs,branchtorus(14)]=SetupTorusBifurcation(funcs, init_torus(14),11,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus(14).method.point.print_residual_info=0;
branchtorus(14).method.point.adapt_mesh_after_correct=1;
branchtorus(14)=br_contn(foldfuncs,branchtorus(14),200);
branchtorus(14)=br_rvers(branchtorus(14));
branchtorus(14)=br_contn(foldfuncs,branchtorus(14),800);
br_plot(branchtorus(14),xm, ym,'k')
%%
[foldfuncs,branchtorus(15)]=SetupTorusBifurcation(funcs, init_torus(15),31,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,150],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus(15).method.point.print_residual_info=0;
branchtorus(15).method.point.adapt_mesh_after_correct=1;
branchtorus(15)=br_contn(foldfuncs,branchtorus(15),200);
branchtorus(15)=br_rvers(branchtorus(15));
branchtorus(15)=br_contn(foldfuncs,branchtorus(15),800);
br_plot(branchtorus(15),xm, ym,'k')
%%
[foldfuncs,branchtorus(16)]=SetupTorusBifurcation(funcs, init_torus(16),86,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,150],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus(16).method.point.print_residual_info=0;
branchtorus(16).method.point.adapt_mesh_after_correct=1;
branchtorus(16)=br_contn(foldfuncs,branchtorus(16),200);
branchtorus(16)=br_rvers(branchtorus(16));
branchtorus(16)=br_contn(foldfuncs,branchtorus(16),200);
br_plot(branchtorus(16),xm, ym,'k')
%%
[foldfuncs,branchtorus(17)]=SetupTorusBifurcation(funcs, init_torus(17),61,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,150],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus(17).method.point.print_residual_info=0;
branchtorus(17).method.point.adapt_mesh_after_correct=1;
branchtorus(17)=br_contn(foldfuncs,branchtorus(17),400);
branchtorus(17)=br_rvers(branchtorus(17));
branchtorus(17)=br_contn(foldfuncs,branchtorus(17),400);
br_plot(branchtorus(17),xm, ym,'m')
%%
[foldfuncs,branchtorus(18)]=SetupTorusBifurcation(funcs, init_torus(18),59,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,150],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus(18).method.point.print_residual_info=0;
branchtorus(18).method.point.adapt_mesh_after_correct=1;
branchtorus(18)=br_contn(foldfuncs,branchtorus(18),1600);
branchtorus(18)=br_rvers(branchtorus(18));
branchtorus(18)=br_contn(foldfuncs,branchtorus(18),200);
br_plot(branchtorus(18),xm, ym,'k')
%%
[foldfuncs,branchtorus(19)]=SetupTorusBifurcation(funcs, init_torus(19),82,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,150],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus(19).method.point.print_residual_info=0;
branchtorus(19).method.point.adapt_mesh_after_correct=1;
branchtorus(19)=br_contn(foldfuncs,branchtorus(19),400);
branchtorus(19)=br_rvers(branchtorus(19));
branchtorus(19)=br_contn(foldfuncs,branchtorus(19),1200);
br_plot(branchtorus(19),xm, ym,'k')
%%
[foldfuncs,branchtorus(20)]=SetupTorusBifurcation(funcs, init_torus(20),29,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,150],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus(20).method.point.print_residual_info=0;
branchtorus(20).method.point.adapt_mesh_after_correct=1;
branchtorus(20)=br_contn(foldfuncs,branchtorus(20),200);
branchtorus(20)=br_rvers(branchtorus(20));
branchtorus(20)=br_contn(foldfuncs,branchtorus(20),1200);
br_plot(branchtorus(20),xm, ym,'k')
%%
[foldfuncs,branchtorus(21)]=SetupTorusBifurcation(funcs, init_torus(21),621,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,150; 17,150],...
    'max_step',[14,0.1; 17,0.1]);
branchtorus(21).method.point.print_residual_info=0;
branchtorus(21).method.point.adapt_mesh_after_correct=1;
branchtorus(21)=br_contn(foldfuncs,branchtorus(21),200);
branchtorus(21)=br_rvers(branchtorus(21));
branchtorus(21)=br_contn(foldfuncs,branchtorus(21),1200);
br_plot(branchtorus(21),xm, ym,'m')
