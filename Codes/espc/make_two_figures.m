close all
bdwidth = 5;
topbdwidth = 30;
set(0,'Units','pixels') 
scnsize = get(0,'ScreenSize');
pos1  = [3*bdwidth,... 
	.45*scnsize(4) + bdwidth,...
	scnsize(3)/2 - 5*bdwidth,...
	scnsize(4)/2 - (topbdwidth + bdwidth)];
pos2 = [pos1(1) + scnsize(3)/2,...
	pos1(2),...
	pos1(3),...
	pos1(4)];
figure('Position',pos1);
figure('Position',pos2);

%fig1_handle=figure(1);
%pause(.2);
%fig2_handle=figure(2);
%pause(.2);
%set(fig1_handle,'Position',pos1);
%set(fig2_handle,'Position',pos2);
