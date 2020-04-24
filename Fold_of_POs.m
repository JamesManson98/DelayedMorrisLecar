%% Fold of POs
clear
close all
addpath( genpath('../../../dde_biftool_v3.1.1'));
myfuncs
load('Two_par_Kappa_Tau_vars','branch2','xm','ym','GH_points','HH_points')
load('initialpoint_pofold','branch4a','branch4b','branch4c','branch4d','branch4e','branch4f','branch4g','branch4h','branch4i','branch4j','branch4k','branch4l','branch4m','branch2negativekappa','branch2negativetau','branch2positivekappa','branch2positivetau')
xlabel('\tau', 'FontSize', 18);ylabel('\kappa', 'FontSize', 18);
hold on
for i=1:length(GH_points)
    plot(branch2.point(GH_points(i)).parameter(17),branch2.point(GH_points(i)).parameter(14),'Marker','x','color', 'k')
end
for i=1:length(HH_points)
    plot(branch2.point(HH_points(i)).parameter(17),branch2.point(HH_points(i)).parameter(14),'Marker','o','color', 'k')
end
plot(branch2negativetau(:),branch2negativekappa(:),'color','g')
plot(branch2positivetau(:),branch2positivekappa(:),'color','r')
xlim([0,100]);
ylim([0,100]);
h(1)=plot(1000,1000, 'g');
h(2)=plot(1000,1000, 'r');
% h(3)=plot(1000,1000, 'b');
h(1)=plot(1000,1000, 'Marker', 'x', 'color', 'k', 'linestyle', 'none');
h(2)=plot(1000,1000, 'Marker', 'o', 'color', 'k', 'linestyle', 'none');
% a=legend(h, { 'GH points', 'HH points'});%'Stable Hopf continuation','Unstable Hopf continuation','Fold of POs'},'FontSize', 18)
% a.Position(1)=0.85;
% a.Position(2)=0.7;
%%
[foldfuncs,branch7a]=SetupPOfold(funcs,branch4a,33,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branch7a.method.point.print_residual_info=0;
branch7a.method.point.adapt_mesh_after_correct=1;
branch7a=br_contn(foldfuncs,branch7a,650);
branch7a=br_rvers(branch7a);
branch7a=br_contn(foldfuncs,branch7a,400);
branch7a=br_rvers(branch7a);
branch7a.parameter.max_step=[14,0.01; 17,0.01];
branch7a=br_contn(foldfuncs,branch7a,350);
br_plot(branch7a,ym, xm,'b')
%%
[foldfuncs,branch7b]=SetupPOfold(funcs,branch4b,28,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branch7b.method.point.print_residual_info=0;
branch7b.method.point.adapt_mesh_after_correct=1;
branch7b=br_contn(foldfuncs,branch7b,370);
branch7b=br_rvers(branch7b);
branch7b=br_contn(foldfuncs,branch7b,500);
br_plot(branch7b,ym, xm,'b')
%%
[foldfuncs,branch7c]=SetupPOfold(funcs,branch4c,17,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branch7c.method.point.print_residual_info=0;
branch7c.method.point.adapt_mesh_after_correct=1;
branch7c=br_contn(foldfuncs,branch7c,50);
branch7c=br_rvers(branch7c);
branch7c=br_contn(foldfuncs,branch7c,900);
br_plot(branch7c,ym, xm,'b')
%%
[foldfuncs,branch7e]=SetupPOfold(funcs,branch4e,25,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,100],...
    'max_step',[14,0.01; 17,0.01]);
branch7e.method.point.print_residual_info=0;
branch7e.method.point.adapt_mesh_after_correct=1;
branch7e=br_contn(foldfuncs,branch7e,100);
branch7e=br_rvers(branch7e);
branch7e=br_contn(foldfuncs,branch7e,500);
branch7e.parameter.max_step=[14,0.1; 17,0.1];
branch7e=br_contn(foldfuncs,branch7e,2000); 
br_plot(branch7e,ym, xm,'b')
%%
[foldfuncs,branch7h]=SetupPOfold(funcs,branch4h,25,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branch7h.method.point.print_residual_info=0;
branch7h.method.point.adapt_mesh_after_correct=1;
branch7h=br_contn(foldfuncs,branch7h,50);
branch7h=br_rvers(branch7h);
branch7h=br_contn(foldfuncs,branch7h,800); 
br_plot(branch7h,ym, xm,'b')
%%
[foldfuncs,branch7d]=SetupPOfold(funcs,branch4d,20,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branch7d.method.point.print_residual_info=0;
branch7d.method.point.adapt_mesh_after_correct=1;
branch7d=br_contn(foldfuncs,branch7d,50);
branch7d=br_rvers(branch7d);
branch7d=br_contn(foldfuncs,branch7d,800);
br_plot(branch7d,ym, xm,'b')
%%
[foldfuncs,branch7i]=SetupPOfold(funcs,branch4i,27,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branch7i.method.point.print_residual_info=0;
branch7i.method.point.adapt_mesh_after_correct=1;
branch7i=br_contn(foldfuncs,branch7i,50);
branch7i=br_rvers(branch7i);
branch7i=br_contn(foldfuncs,branch7i,400);
br_plot(branch7i,ym, xm,'b')
%%
[foldfuncs,branch7j]=SetupPOfold(funcs,branch4j,16,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branch7j.method.point.print_residual_info=0;
branch7j.method.point.adapt_mesh_after_correct=1;
branch7j=br_contn(foldfuncs,branch7j,50);
branch7j=br_rvers(branch7j);
branch7j=br_contn(foldfuncs,branch7j,800);
br_plot(branch7j,ym, xm,'b')
%%
[foldfuncs,branch7k]=SetupPOfold(funcs,branch4k,17,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branch7k.method.point.print_residual_info=0;
branch7k.method.point.adapt_mesh_after_correct=1;
branch7k=br_contn(foldfuncs,branch7k,50);
branch7k=br_rvers(branch7k);
branch7k=br_contn(foldfuncs,branch7k,800);
br_plot(branch7k,ym, xm,'b')
%%
[foldfuncs,branch7l]=SetupPOfold(funcs,branch4l,15,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,100],...
    'max_step',[14,0.1; 17,0.1]);
branch7l.method.point.print_residual_info=0;
branch7l.method.point.adapt_mesh_after_correct=1;
branch7l=br_contn(foldfuncs,branch7l,50);
branch7l=br_rvers(branch7l);
branch7l=br_contn(foldfuncs,branch7l,800);
br_plot(branch7l,ym, xm,'b')
%%
[foldfuncs,branch7m]=SetupPOfold(funcs,branch4m,14,'contpar',[14,17],...
    'dir',17,'step',0.001,'min_bound',[14,0; 17,0],'max_bound',[14,100; 17,150],...
    'max_step',[14,0.1; 17,0.1]);
branch7m.method.point.print_residual_info=0;
branch7m.method.point.adapt_mesh_after_correct=1;
branch7m=br_contn(foldfuncs,branch7m,50);
branch7m=br_rvers(branch7m);
branch7m=br_contn(foldfuncs,branch7m,800);
br_plot(branch7m,ym, xm,'b')