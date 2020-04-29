%% makes the complete bifurcation diagram
%loading in variables
clear all
load('Fold_of_POs')
load('Torus_bifurcations')

%plotting sub/supercritical hopfs
figure(10)
xlabel('\tau', 'FontSize', 18);ylabel('\kappa', 'FontSize', 18);
hold on
plot(hopf_br_tk_stable_tau(:),hopf_br_tk_stable_kappa(:),'color','g')
plot(hopf_br_tk_unstable_tau(:),hopf_br_tk_unstable_kappa(:),'color','r')
for i=1:length(GH_points)
    plot(hopf_br_tk.point(GH_points(i)).parameter(17),hopf_br_tk.point(GH_points(i)).parameter(14),'Marker','x','color', 'k')
end
for i=1:length(HH_points)
    plot(hopf_br_tk.point(HH_points(i)).parameter(17),hopf_br_tk.point(HH_points(i)).parameter(14),'Marker','o','color', 'k')
end
xlim([0,100]);
ylim([0,100]);
xm.field='parameter';
xm.col=17;
ym.field='parameter';
ym.col=14;
%% plotting fold bifurcations of periodic orbits
for i=1:length(PO_fold)
br_plot(PO_fold(i),xm,ym,'b')
end
%% plotting torus bifurcations that change stability
for i=[3,9,13,17,21]
    br_plot(branchtorus(i),xm, ym,'k')
end

%lines to show constant kappa slices

line([0 100], [85 85], 'LineStyle','--', 'color', 'k')
line([0 100], [45 45], 'LineStyle','--', 'color', 'k')
line([0 100], [15 15], 'LineStyle','--', 'color', 'k')
