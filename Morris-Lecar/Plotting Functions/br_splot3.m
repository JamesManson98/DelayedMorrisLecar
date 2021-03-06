function br_splot3(br,xm,ym,zm)
%% plot 3d branch stability
% function br_splot3(branch,x_measure,y_measure,z_measure,line_type)
% INPUT:
%	branch branch of points
%	x_measure scalar measure for x-coordinate
%	y_measure scalar measure for y-coordinate
%   z_measure scalar measure for z-coordinate

i_start=1;

while i_start<=length(br.point)
  s=p_dststb(br.point(i_start));
  i_end=i_start;
  while p_dststb(br.point(i_end))==s
    i_end=i_end+1;
    if i_end>length(br.point) 
        break; 
    end
  end
  if i_end>i_start+1
    b=br;
    b.point=b.point(i_start:i_end-1);
    if s==1, br_plot3(b,xm,ym,zm,'g'); else br_plot3(b,xm,ym,zm,'r'); end;
    hold on;
  end;
  i_start=i_end;
end;

for i=1:length(br.point)-1
  [s1,d1,c1]=p_dststb(br.point(i));
  [s2,d2,c2]=p_dststb(br.point(i+1));
  x1=p_measur(br.point(i),xm);
  x2=p_measur(br.point(i+1),xm);
  y1=p_measur(br.point(i),ym);
  y2=p_measur(br.point(i+1),ym);
  z1=p_measur(br.point(i),zm);
  z2=p_measur(br.point(i+1),zm);
  x=(abs(d1)*x2+abs(d2)*x1)/(abs(d1)+abs(d2));
  y=(abs(d1)*y2+abs(d2)*y1)/(abs(d1)+abs(d2));
  z=(abs(d1)*z2+abs(d2)*z1)/(abs(d1)+abs(d2));
  if s1~=s2
    if s1==1 
        plot3([x1 x],[y1 y],[z1 z],'g'); 
        plot3([x x2],[y y2],[z z2],'r');
    elseif s2==1
        plot3([x1 x],[y1 y],[z1 z],'r'); 
        plot3([x x2],[y y2],[z z2],'g'); 
    else
        plot3([x1 x2],[y1 y2],[z1 z2],'r'); 
    end
  end
  if s1~=s2
    if abs(d1)>abs(d2) 
        c1=c2; 
    end
    switch c1,
      case 0, plot3(x,y,z,'kx');
      case 1, plot3(x,y,z,'ko','MarkerFaceColor', 'k');
      case -1, plot3(x,y,z,'kd');
    end
  end
end
end
