%% Plotting phase portraits and time series at multistable taus
clear all
addpath(genpath('../../../../../dde_biftool_v3.1.1'));
load('One_Kappa_POs_15','branches_PO_kappa_15')
myfuncs

%% setting up tau/branches to be used
kappa_15(1).tau=34;
kappa_15(2).tau=65;
kappa_15(3).tau=99.5;

kappa_15(1).ind=[2,4];
kappa_15(2).ind=[5,7];
kappa_15(3).ind=[7,9];
%%
for i=1:length(kappa_15)
    c=1;
    tau=kappa_15(i).tau;
    subplot(3,3,i)
    figure(1)
    hold on
    for j=1:length(kappa_15(i).ind)
        br_ind=kappa_15(i).ind(j);
        k=find(diff(sign(arrayfun(@(x) x.parameter(17), branches_PO_kappa_15(br_ind).branch.point)-tau))~=0);
        for z=1:length(k)
            if p_dststb(branches_PO_kappa_15(br_ind).branch.point(k(z)))==1
                plot(branches_PO_kappa_15(br_ind).branch.point(k(z)).profile(1,:), branches_PO_kappa_15(br_ind).branch.point(k(z)).profile(2,:),'g')
                [M,I]=min(branches_PO_kappa_15(br_ind).branch.point(k(z)).profile(1,:));
                h(c)=plot(branches_PO_kappa_15(br_ind).branch.point(k(z)).profile(1,I),branches_PO_kappa_15(br_ind).branch.point(k(z)).profile(2,I),'x');
                figure(2)
                subplot(3,3,i)
                period=branches_PO_kappa_15(br_ind).branch.point(k(z)).period;
                time=branches_PO_kappa_15(br_ind).branch.point(k(z)).mesh*period;
                V=branches_PO_kappa_15(br_ind).branch.point(k(z)).profile(1,:);
                plot([time,time+period,time+2*period,time+3*period,time+4*period], [V,V,V,V,V]);
                hold on
                xlabel('Time (ms)');ylabel('V');
                xlim([0,100])
                figure(1)
            else
                plot(branches_PO_kappa_15(br_ind).branch.point(k(z)).profile(1,:), branches_PO_kappa_15(br_ind).branch.point(k(z)).profile(2,:),'r')
                [M,I]=min(branches_PO_kappa_15(br_ind).branch.point(k(z)).profile(1,:));
                h(c)=plot(branches_PO_kappa_15(br_ind).branch.point(k(z)).profile(1,I),branches_PO_kappa_15(br_ind).branch.point(k(z)).profile(2,I),'*');
            end
            string{c}=[num2str(round(branches_PO_kappa_15(br_ind).branch.point(k(z)).period,2))];
            c=c+1;
            xlabel('V');ylabel('W');
            par=[0.125, 4.0, 2.0, 8.0, -2.8, 26.0, 12.0, 17.4, -60.0, 120.0, -91.893652410657, 20.0, 15.0343886973906, 15.0, 0.0, 5.0, tau];
            fun = @(x)rhs(x,par);
            x0 = [2,0.23];
            x = fsolve(fun,x0); %solving for the fixed point using fsolve
            plot(x(1),x(2),'ko','MarkerFaceColor', 'k')
        end
    end
    a(i)=legend(h, string,'FontSize', 10,'location','NorthWest');
    a(i).Position(1)=a(i).Position(1)-0.00875;
    a(i).Position(2)=a(i).Position(2)+0.013;
end