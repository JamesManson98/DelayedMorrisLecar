%%
clear all
% close all
addpath(genpath('../../../../../dde_biftool_v3.1.1'));
load('One_Kappa_POs_85','branches_PO_kappa_85')
myfuncs

%%
kappa_85(1).tau=34;
kappa_85(3).tau=99.5;
kappa_85(2).tau=66;

kappa_85(1).ind=[2,4];
kappa_85(3).ind=[5,7,9,11,12];
kappa_85(2).ind=[4,6,8,12];
%%
for i=1:length(kappa_85)
    tau=kappa_85(i).tau;
    c=1;
    subplot(3,3,i+6)
    hold on
    for j=1:length(kappa_85(i).ind)
        br_ind=kappa_85(i).ind(j);
        k=find(diff(sign(arrayfun(@(x) x.parameter(17), branches_PO_kappa_85(br_ind).branch.point)-tau))~=0);
        for z=1:length(k)
            if p_dststb(branches_PO_kappa_85(br_ind).branch.point(k(z)))==1
                plot(branches_PO_kappa_85(br_ind).branch.point(k(z)).profile(1,:), branches_PO_kappa_85(br_ind).branch.point(k(z)).profile(2,:),'g')
                if i==2
                    [M,I]=min(branches_PO_kappa_85(br_ind).branch.point(k(z)).profile(2,:));
                else
                    [M,I]=min(branches_PO_kappa_85(br_ind).branch.point(k(z)).profile(1,:));
                end
                h(c)=plot(branches_PO_kappa_85(br_ind).branch.point(k(z)).profile(1,I),branches_PO_kappa_85(br_ind).branch.point(k(z)).profile(2,I),'x');
                figure(2)
                subplot(3,3,i+6)
                period=branches_PO_kappa_85(br_ind).branch.point(k(z)).period;
                time=branches_PO_kappa_85(br_ind).branch.point(k(z)).mesh*period;
                V=branches_PO_kappa_85(br_ind).branch.point(k(z)).profile(1,:);
                plot([time,time+period,time+2*period,time+3*period,time+4*period], [V,V,V,V,V]);
                hold on
                xlabel('Time (ms)');ylabel('V');
                xlim([0,100])
                figure(1)                
            else
                plot(branches_PO_kappa_85(br_ind).branch.point(k(z)).profile(1,:), branches_PO_kappa_85(br_ind).branch.point(k(z)).profile(2,:),'r')
                if i==2
                    [M,I]=min(branches_PO_kappa_85(br_ind).branch.point(k(z)).profile(2,:));
                else
                    [M,I]=min(branches_PO_kappa_85(br_ind).branch.point(k(z)).profile(1,:));
                end
                h(c)=plot(branches_PO_kappa_85(br_ind).branch.point(k(z)).profile(1,I),branches_PO_kappa_85(br_ind).branch.point(k(z)).profile(2,I),'*');
            end
            string{c}=[num2str(round(branches_PO_kappa_85(br_ind).branch.point(k(z)).period,2))];
            c=c+1;
            xlabel('V');ylabel('W');
            par=[0.125, 4.0, 2.0, 8.0, -2.8, 26.0, 12.0, 17.4, -60.0, 120.0, -91.893652410657, 20.0, 15.0343886973906, 85.0, 0.0, 5.0, tau];
            fun = @(x)rhs(x,par);
            x0 = [2,0.23];
            x = fsolve(fun,x0); %solving for the fixed point using fsolve
            plot(x(1),x(2),'ko','MarkerFaceColor', 'k')
        end
    end
    a3(i)=legend(h, string,'FontSize', 10,'location','NorthWest');
    a3(i).Position(1)=a3(i).Position(1)-0.00875;
    a3(i).Position(2)=a3(i).Position(2)+0.013;
end
%Fix bug in DDEBiftool coming off swallowtail and plots with 1d slices and
%the branches for the torus
%Try and find torus
%So if a change in stability is not at a fold
%Check floquet multipliers a pair of complex conjugates on unit circle