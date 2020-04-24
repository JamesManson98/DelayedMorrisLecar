point=branches_PO_kappa_45(6).branch.point(270); %193
newmesh=linspace(0,1,201);
rm_point = p_remesh (point , 4 , newmesh );
period=rm_point.period;
tau=rm_point.parameter(17);
par=rm_point.parameter;
mesh=rm_point.mesh;
tspan=mesh*period;
history=rm_point.profile;
history_fun = @(t) interp1(tspan'-period, history', -mod(-t,period))';
sol = dde23(@(t,xx,yy)sys_rhs_test(t,xx,yy,par), [tau], history_fun, [0, 120*period]);
subplot(2,1,2)
plot(sol.x(:),sol.y(1,:))
xlim([0,2200])
xlabel('Time (ms)', 'FontSize', 18);ylabel('V', 'FontSize', 18);