function [X,Y,Z]=Generate_Grid_3D_uniform(nx,ny,NO)

xmin=-1;
xmax=1.0;
ymin=-1;
ymax = 1;
% res=0.5;
% NO=4;
NO2=NO/2;
% nx = 64*res;
dx =(xmax-xmin)/nx;
% ny = 64*res;
dy = (ymax-ymin)/ny;
x=xmin:dx:xmax;
y=ymin:dy:ymax;
% plot(x,tan(x))
% tan(x(1))/tan(x(end));
% 
% scx = 30/tan(x(end));
% scy = 60/tan(y(end));

% x_grid = scx*tan(x);
% y_grid = scy*tan(y);

x_grid = x;
y_grid = y;

% add ghost cells 
x_ghost_left=0;
x_ghost_right=0;
count=1;
for i=NO2:-1:1
x_ghost_left(i) = 2*x_grid(1) - x_grid(1+count);
x_ghost_right(count) = 2*x_grid(end) - x_grid(end-count);
count=count+1;
end

count=1;
y_ghost_left=0;
y_ghost_right=0;
for i=NO2:-1:1
y_ghost_left(i) = 2*y_grid(1) - y_grid(1+count);
y_ghost_right(count) = 2*y_grid(end) - y_grid(end-count);
count=count+1;
end

x_grid_all = [x_ghost_left,x_grid,x_ghost_right];
y_grid_all = [y_ghost_left,y_grid,y_ghost_right];
z_grid_all = y_grid_all;

X=zeros(length(x_grid_all),length(y_grid_all),length(z_grid_all));
Y=X;
Z=X;

for i=1:length(x_grid_all)
    for j=1:length(y_grid_all)
        for k=1:length(z_grid_all)
            X(i,j,k) = x_grid_all(i);
            Y(i,j,k) = y_grid_all(j);
            Z(i,j,k) = z_grid_all(k);
        end
    end
end