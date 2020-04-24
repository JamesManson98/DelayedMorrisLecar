h=figure;
filename = 'Floquet2.gif';
theta = linspace(0,2*pi,300);
x = cos(theta);
y = sin(theta);
plot(x,y,'k')
hold on
xlabel('Real component')
ylabel('Imaginary component')
for i=1:length(branches_PO_kappa_45(6).branch.point)
    f=branches_PO_kappa_45(6).branch.point(i).profile;
    ampl=max(f(1,:))-min(f(1,:));
    p=plot(real(branches_PO_kappa_45(6).branch.point(i).stability.mu),imag(branches_PO_kappa_45(6).branch.point(i).stability.mu),'ro');
    %drawnow
    axis equal
    title(['Floquet Multiplier at (\tau, Amplitude)= (' num2str(round(branches_PO_kappa_45(6).branch.point(i).parameter(17))) ',' num2str(round(ampl)) ')'])
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if i == 1 
          imwrite(imind,cm,fullfile(cd,filename),'gif','DelayTime',0.01,'Loopcount',inf); 
      else 
          imwrite(imind,cm,fullfile(cd,filename),'gif','DelayTime',0.01,'WriteMode','append'); 
      end
    delete(p)
end
