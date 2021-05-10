% 1D1V Finite Volume Vlasov code for Landau Damping simulations
clear; 

NO = 8; % order of reconstrction / length of reconstruction stencil
NO2 = NO/2; % # of ghost cells on each side

% 1. generate a 2-D cartesian grid
nx = 64; % # of active cells in x space
ny = 64; % # of active cells in v space

L = 2*pi;
Vmax = 6;

deltax = 2/nx*L; % delta x
deltay = 2/ny*Vmax; % delta y

x1 = -L-NO2*deltax:deltax:L+NO2*deltax;
y1 = -Vmax-NO2*deltay:deltay:Vmax+NO2*deltay;

x=zeros(nx+1,ny+1); % 2-D mesh including corners
y=zeros(nx+1,ny+1);

% generate the corner grid including ghost cells
for i=1:length(x1)
    for j=1:length(y1)
            x(i,j) = x1(i);
            y(i,j) = y1(j);
    end
end
x = x+L;

[nx_total,ny_total]=size(x); 

% 2. index/metric calculations

% I,J are the indices for all cell centers
I = 1:nx_total-1;
J = 1:ny_total-1;
% Ip1,Jp1 are the indices for all cell corners
Ip1 = 1:nx_total;
Jp1 = 1:ny_total;

% index of active cell centers
Ic = NO2+1:NO2+nx;
Jc = NO2+1:NO2+ny;
% index of active face centers
If = NO2+1:NO2+nx+1;
Jf = NO2+1:NO2+ny+1;

% index of left-boudary for cell centers
Ic_lb = 1:NO2;
Jc_lb = 1:NO2;
% index of right-boundary for cell centers
Ic_rb = nx+NO2+1: nx+NO2+NO2;
Jc_rb = ny+NO2+1: ny+NO2+NO2;

% locations for cell centers
xc(I,J) = 0.25*( x(I,J) + x(I+1,J) + x(I,J+1) + x(I+1,J+1) );
yc(I,J) = 0.25*( y(I,J) + y(I+1,J) + y(I,J+1) + y(I+1,J+1) );

% locations for x-interfaces 
xi(Ip1,J) = 0.5*( x(Ip1,J) + x(Ip1,J+1) );
yi(Ip1,J) = 0.5*( y(Ip1,J) + y(Ip1,J+1) );

% locations for y-interfaces 
xj(I,Jp1) = 0.5*( x(I,Jp1) + x(I+1,Jp1) );
yj(I,Jp1) = 0.5*( y(I,Jp1) + y(I+1,Jp1) );

% cell length dx and dy
dx(I,J) = x(I+1,J) - x(I,J);
dy(I,J) = y(I,J+1) - y(I,J);

% calculate cell "volume" as dx times dy, which is used in reconstruction
volume(I,J) = dx(I,J).*dy(I,J);

% 3. Initial conditions
% two-stream (x - configuration space, y - velocity space)
A = 0.05;
k = 0.5;
% rho = 1/sqrt(2*pi).*yc.^2.*exp(-yc.^2/2).*(1+A.*cos(k.*xc));
rho = 1/sqrt(2*pi).*exp(-yc.^2/2).*(1+A.*cos(k.*xc));

% rho = 2/7/sqrt(2*pi).*(1+5*yc.^2).*exp(-yc.^2/2).*(1+A.*((cos(2*k.*xc)+cos(3*k.*xc))/1.2 + cos(k*xc)));
Ex = xc(:,1)*0;
qe = xc(:,1)*0;

% compute initial dt
CFL = 0.3;
PDMB = 4;

dtCFL = CFL ./(abs(yc(Ic,1))./deltax+abs(Ex(Ic))./deltay+eps);
dt = min(dtCFL(:));

% Adam-Bashforth predictor step
% Save the initial states for the first Adam-Bashforth time stepping
rho_p = rho; 
dt0 = dt;

Time = 0; % simulation time 

rho_interp = x*0;% interpolated value (for reconstruction)
rho_left = x*0;  % left state (after reconstruction)
rho_right = x*0; % right state (after reconstruction)
rho_flux_x = x*0;% x-interface flux (for finite volume update)
rho_flux_y = y*0;% y-interface flux (for finite volume update)
ex = xc*0;

count = 1;
% 4. loop over steps
for n=1:1000000 
    
    % step 0: update the time step based on the CFL condition
    dtCFL = CFL ./(abs(yc(Ic,1))./deltax+abs(Ex(Ic))./deltay+eps);
    dt = min(dtCFL(:));    
    
    % Step 1: predictor - get the half time step values
    rho_h = rho + dt./dt0./2.*(rho-rho_p); % 2nd-order extrapolation
    % save the current rho for the next AB time stepping
    rho_p = rho; 
    dt0 = dt;
    % advance simulation time
    Time = Time+dt;
    
    % compute electric field
    qe(Ic) = 1 - sum(rho_h(Ic,Jc).*deltay,2); % charge density
    
    Ex(:) = 0;
    for i = Ic(1):Ic(end)-1
        Ex(i+1) = Ex(i) + 0.5*(qe(i)+qe(i+1))*deltax; % integrate
    end

    Ex(Ic) = Ex(Ic) - mean(Ex(Ic)); % remove mean E field
    
    for j = Jc(1):Jf(end)
        ex(Ic,j) = - Ex(Ic);
    end
    
    rho_h0 = rho_h;        % density - for limiting

    if_act =If;
    jf_act = Jc;

%     % Step 2.1: reconstruct mass in the x-direction - 7th order
    % upwind interpolation and limiting for left-state
    rho_interp(if_act,jf_act) = -1/140*rho_h(if_act-4,jf_act)+5/84*rho_h(if_act-3,jf_act)-101/420*rho_h(if_act-2,jf_act)+319/420*rho_h(if_act-1,jf_act)+...
        107/210*rho_h(if_act,jf_act)-19/210*rho_h(if_act+1,jf_act)+1/105*rho_h(if_act+2,jf_act);
    %rho_interp(if_act,jf_act) = rho_interp(if_act,jf_act)./vol_jface(if_act,jf_act);
    rho_left(if_act,jf_act) = PDMU77(rho_h0(if_act-4,jf_act),rho_h0(if_act-3,jf_act),rho_h0(if_act-2,jf_act),rho_h0(if_act-1,jf_act),rho_h0(if_act,jf_act),rho_h0(if_act+1,jf_act),rho_h0(if_act+2,jf_act),rho_interp(if_act,jf_act),PDMB);
    
    % upwind interpolation and limiting for right-state
    rho_interp(if_act,jf_act) = -1/140*rho_h(if_act+3,jf_act)+5/84*rho_h(if_act+2,jf_act)-101/420*rho_h(if_act+1,jf_act)+319/420*rho_h(if_act,jf_act)+...
        107/210*rho_h(if_act-1,jf_act)-19/210*rho_h(if_act-2,jf_act)+1/105*rho_h(if_act-3,jf_act);
    %rho_interp(if_act,jf_act) = rho_interp(if_act,jf_act)./vol_jface(if_act,jf_act);
    rho_right(if_act,jf_act)= PDMU77(rho_h0(if_act+3,jf_act),rho_h0(if_act+2,jf_act),rho_h0(if_act+1,jf_act),rho_h0(if_act,jf_act),rho_h0(if_act-1,jf_act),rho_h0(if_act-2,jf_act),rho_h0(if_act-3,jf_act),rho_interp(if_act,jf_act),PDMB);
    
    % Step 2.2: mass flux in the x-direction; NOTE: since interface velocity is known (vxi), i'm just using a simple upwind flux
    rho_flux_x(If,Jc) = ( (1+sign(yc(If,Jc)))./2.* rho_left(If,Jc) + (1-sign(yc(If,Jc)))./2.* rho_right(If,Jc) ).*yc(If,Jc);

    if_act =Ic;
    jf_act = Jf;

    % upwind interpolation and limiting for left-state
    rho_interp(if_act,jf_act) = -1/140*rho_h(if_act,jf_act-4)+5/84*rho_h(if_act,jf_act-3)-101/420*rho_h(if_act,jf_act-2)+319/420*rho_h(if_act,jf_act-1)+...
        107/210*rho_h(if_act,jf_act)-19/210*rho_h(if_act,jf_act+1)+1/105*rho_h(if_act,jf_act+2);
    %rho_interp(if_act,jf_act) = rho_interp(if_act,jf_act)./vol_jface(if_act,jf_act);
    rho_left(if_act,jf_act) = PDMU77(rho_h0(if_act,jf_act-4),rho_h0(if_act,jf_act-3),rho_h0(if_act,jf_act-2),rho_h0(if_act,jf_act-1),rho_h0(if_act,jf_act),rho_h0(if_act,jf_act+1),rho_h0(if_act,jf_act+2),rho_interp(if_act,jf_act),PDMB);
    
    % upwind interpolation and limiting for right-state
    rho_interp(if_act,jf_act) = -1/140*rho_h(if_act,jf_act+3)+5/84*rho_h(if_act,jf_act+2)-101/420*rho_h(if_act,jf_act+1)+319/420*rho_h(if_act,jf_act)+...
        107/210*rho_h(if_act,jf_act-1)-19/210*rho_h(if_act,jf_act-2)+1/105*rho_h(if_act,jf_act-3);
    %rho_interp(if_act,jf_act) = rho_interp(if_act,jf_act)./vol_jface(if_act,jf_act);
    rho_right(if_act,jf_act)= PDMU77(rho_h0(if_act,jf_act+3),rho_h0(if_act,jf_act+2),rho_h0(if_act,jf_act+1),rho_h0(if_act,jf_act),rho_h0(if_act,jf_act-1),rho_h0(if_act,jf_act-2),rho_h0(if_act,jf_act-3),rho_interp(if_act,jf_act),PDMB);


    % Step 2.4: mass flux in the y-direction: use interface velocity vyj
    rho_flux_y(Ic,Jf) = ( (1+sign(ex(Ic,Jf)))./2.* rho_left(Ic,Jf) + (1-sign(ex(Ic,Jf)))./2.* rho_right(Ic,Jf) ).*ex(Ic,Jf);
    
    % step 3: finite volume update: delta(flux*face)/volume = delta(flux)/length
    rho(Ic,Jc) = rho_p(Ic,Jc) - dt.*( (rho_flux_x(Ic+1,Jc)-rho_flux_x(Ic,Jc))./dx(Ic,Jc) +  ...
                                      (rho_flux_y(Ic,Jc+1)-rho_flux_y(Ic,Jc))./dy(Ic,Jc) );

     % Step 4: boundary conditions 
     % --periodic boundary - x direction
    rho(Ic_rb,:) = rho(Ic(1):Ic(NO2),:);
    rho(Ic_lb,:) = rho(Ic(end-NO2+1):Ic(end),:);
    % --periodic boundary - y direction;
    rho(:,Jc_lb) = rho(:,Jc(end-NO2+1):Jc(end));
    rho(:,Jc_rb) = rho(:,Jc(1):Jc(NO2));   

    % end of calculations, now make some plots at run-time
    if(mod(n,20)==0)
        figure(1)
        contourf(xc(Ic,Jc)/pi,yc(Ic,Jc),rho(Ic,Jc),50,'linestyle','none');shading flat;%caxis([0 1])
        title(['f(x,v), Time = ',num2str(round(Time,2))]);
        colormap(turbo());
        xlabel('x')
        ylabel('v')
        set(gca,'fontsize',14);
        colorbar;
        pause(0.05);
        disp(['Loop ',num2str(n),', Sim Time = ',num2str(Time)]);
        
%         str = sprintf('%0*d',4,count);
%         image_fn = ['~/Data/solver/png/vlasov_128x256-',str];
%         disp(image_fn);
%         saveas(gcf,image_fn,'png'); % save image
        count = count+1;
        eee(count) = max(abs(Ex));
        tim(count) = Time;
    end 
    
    if(Time>50) % break if T>50
        break;
    end
end
%% figure compare the damping rate
figure
semilogy(tim,eee)
hold on
gam = 0.153;
plot(tim,A*exp(-gam*tim))

