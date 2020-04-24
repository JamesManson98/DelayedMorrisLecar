%behaviour of f(V(t-tau))=(Kappa/2)*(1+tanh((V(t-tau)-Vs)/Kappa_s)) for a variety of
%Kappa and Kappa_s, (Vs=0)
Kappa_s=[0.1,1,10];
Kappa=[1,5];
V_delayed=-60:0.01:60;
colours=['r','b','g'];
for i=1:length(Kappa)
    subplot(1,length(Kappa),i)
    if i==1
        subplot(1,length(Kappa),1)
        title('(a)')
    else
        subplot(1,length(Kappa),2)
        title('(b)')
    end
    hold on
    for j=1:length(Kappa_s)
        plot(V_delayed,(Kappa(i)/2)*(1+tanh(V_delayed/Kappa_s(j))),colours(j))
    end
    xlabel('V(t-\tau)')
    ylabel('F(V(t-\tau))')
    legend({['\kappa_s = ' num2str(Kappa_s(1))], ...
        ['\kappa_s = ' num2str(Kappa_s(2))], ...
        ['\kappa_s = ' num2str(Kappa_s(3))]}, 'FontSize',18)
    set(gca,'FontSize',18)
end