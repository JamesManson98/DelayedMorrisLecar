function br_plot3(branch,x_m,y_m,z_m,line_type)
%% plot 3d branch
% function br_plot3(branch,x_measure,y_measure,z_measure,line_type)
% INPUT:
%	branch branch of points
%	x_measure scalar measure for x-coordinate
%	y_measure scalar measure for y-coordinate
%   z_measure scalar measure for z-coordinate
%	line_type line type to plot with

% (c) DDE-BIFTOOL v. 3.0(19), 11/04/2014
%
% 
%
if ~exist('line_type','var')
  line_type='';
end;
if isempty(line_type)
  line_type='';
end;

ll=length(branch.point);

if ll<=1
  error('BR_PLOT: ll=%d, branch contains too few points.',ll);
end;

if ~isempty(x_m)
  [ffx,llx]=br_measr(branch,x_m);
else
  ffx=1:ll;
  llx=ones(ll);
end;
if ~isempty(y_m)
  [ffy,lly]=br_measr(branch,y_m);
else
  ffy=1:ll;
  lly=ones(ll);
end;

  [ffz,llz]=br_measr(branch,z_m);



mx=max(llx);
my=max(lly);

if mx~=min(llx)
  for i=1:mx
    j_end=0;
    while 1
      llx(ll+1)=mx+1;
      j_start=j_end+1;
      while llx(j_start)<i
       j_start=j_start+1;
      end;
      if j_start==ll+1
        break;
      end;
      j_end=j_start+1;
      llx(ll+1)=0;
      while llx(j_end)>=i
        j_end=j_end+1;
      end;
      j_end=j_end-1;
      plot(ffx(j_start:j_end,i),ffy(j_start:j_end),line_type);
      hold on;
    end;
  end;
elseif my~=min(lly)
  for i=1:my
    j_end=0;
    while 1
      lly(ll+1)=my+1;
      j_start=j_end+1;
      while lly(j_start)<i
       j_start=j_start+1;
      end;
      if j_start==ll+1
        break;
      end;
      j_end=j_start+1;
      lly(ll+1)=0;
      while lly(j_end)>=i
        j_end=j_end+1;
      end;
      j_end=j_end-1;
      plot(ffx(j_start:j_end),ffy(j_start:j_end,i),line_type);
      hold on;
    end;
  end;
else
  plot3(ffx,ffy,ffz,line_type); 
  hold on;
end;

return;
