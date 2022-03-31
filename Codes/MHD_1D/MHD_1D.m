%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                1-D HD/MHD solver on Cartesian Grid 
%
%    FEATURES:
%    1. Finite Volume
%    2. TVD/TVB limiter - PDM/WENO
%    2. Semi-conservative - only track plasma energy
%    3. high order reconstruction in space
%    4. 2nd order Adam-Bashforth scheme in time
%    5. Operator-splitting for the Lorentz force
%    6. Beam scheme for plasma and magnetic field 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

% Model Parameters
NO = 8;% order of reconstruction
NO2 = NO/2; % number of ghost cells on each side
gamma=5/3;
CFL = 0.3;
PDMB = 4;

% Grid information
res = 4;
nx = 256;
ny = 4;
nz = ny;
[x,y,z]=Generate_Grid_3D_uniform(nx,ny,NO);
[nx_total,ny_total,nz_total]=size(x);
x = x(:,1,1); % cell edge grid
I = 1:nx_total-1; % total number of cells

xc(I) = 0.5*( x(I) + x(I+1) ); % cell center grid
dx(I) = x(I+1) - x(I);         % dx of each cell

ic_act = NO2+1:NO2+nx;   % index of active cell centers
if_act = NO2+1:NO2+nx+1; % index of active cell faces: face = cell +1

ic_lb = 1:NO2;                % index of left boundary cells
ic_rb = nx+NO2+1: nx+NO2+NO2; % index of right boundar cells
%
% Define premitive Hydrodynamic variables at cell center
rho = zeros(size(xc));
vx = zeros(size(xc));
vy = zeros(size(xc));
vz = zeros(size(xc));
p = zeros(size(xc));

% Define magnetic fields at cell center - for 1-D MHD, the cell centered
% bx, by, bz are the evolved variables as long as bx=const. (which satisfies div B = 0)
bx = zeros(size(xc));
by = zeros(size(xc));
bz = zeros(size(xc));

% Initialize premitive Hydro variables - Brio-Wu shock tube simulation
rho(xc<0) = 1;
rho(xc>=0)= 0.125;
p(xc<0)=1;
p(xc>=0)=.1;
bx(:) = 0.75;
by(xc<0) = 1;
by(xc>=0) = -1;

% % quasi-linear advection
% rho(:)=1;
% rho(abs(xc)<0.1) = 2;
% p(:) = 1;
% vx(:) = 1;

% --symmetric boundary - x direction, left end
rho(ic_rb) = rho(2*nx+NO+1-ic_rb);
p(ic_rb) = p(2*nx+NO+1-ic_rb);
vx(ic_rb) = vx(2*nx+NO+1-ic_rb);
vy(ic_rb) = vy(2*nx+NO+1-ic_rb);
vz(ic_rb) = vz(2*nx+NO+1-ic_rb);
by(ic_rb) = by(2*nx+NO+1-ic_rb);
bz(ic_rb) = bz(2*nx+NO+1-ic_rb);

% --symmetric boundary - x direction, right end
rho(ic_lb) =rho(NO+1-ic_lb);
vx(ic_lb)  = vx(NO+1-ic_lb);
vy(ic_lb)  = vy(NO+1-ic_lb);
vz(ic_lb)  = vz(NO+1-ic_lb);
p(ic_lb)   =  p(NO+1-ic_lb);
by(ic_lb)  = by(NO+1-ic_lb);
bz(ic_lb)  = bz(NO+1-ic_lb);
% Get conserved Hydrodynamic variables
[rho,rhovx,rhovy,rhovz,eng] = getConservedVariables(rho,vx,vy,vz,p,gamma);

% get the first dt
dy=1e20;
dz=1e20;
dt0 = getDT(rho,vx,vy,vz,bx,by,bz,p,gamma,dx,dy,dz,CFL); 

% Save the first step states for Adam-Bashforth time stepping
rho_p = rho;
vx_p = vx;
vy_p = vy;
vz_p = vz;
bx_p = bx;
by_p = by;
bz_p = bz;
p_p= p;

% main loop
count = 0;
Nstep=1e6;
Time = 0;
figure

rho_left = zeros(size(vx));
vx_left = zeros(size(vx));
vy_left = zeros(size(vx));
vz_left = zeros(size(vx));
bx_left = zeros(size(vx));
by_left = zeros(size(vx));
bz_left = zeros(size(vx));
p_left = zeros(size(vx));

rho_right = zeros(size(vx));
vx_right = zeros(size(vx));
vy_right = zeros(size(vx));
vz_right = zeros(size(vx));
bx_right = zeros(size(vx));
by_right = zeros(size(vx));
bz_right = zeros(size(vx));
p_right = zeros(size(vx));

for n=1:Nstep
    
    tic;
    dt = getDT(rho,vx,vy,vz,bx,by,bz,p,gamma,dx,dy,dz,CFL);
    Time=Time+dt;
    [rho,rhovx,rhovy,rhovz,eng] = getConservedVariables(rho,vx,vy,vz,p,gamma);
    
    % Step 1: get the half time step values
    rho_h = rho + dt./dt0./2.*(rho-rho_p);
    vx_h = vx + dt./dt0./2.*(vx-vx_p);
    vy_h = vy + dt./dt0./2.*(vy-vy_p);
    vz_h = vz + dt./dt0./2.*(vz-vz_p);
    p_h = p + dt./dt0./2.*(p-p_p);    
    bx_h = bx + dt./dt0./2.*(bx-bx_p); % Bx is not evloved in 1-D MHD!
    by_h = by + dt./dt0./2.*(by-by_p); % Bx is not evloved in 1-D MHD!
    bz_h = bz + dt./dt0./2.*(bz-bz_p);
    
    rho_p = rho;
    vx_p = vx;
    vy_p = vy;
    vz_p = vz;
    p_p= p;
    by_p = by;
    bz_p = bz;
    dt0 = dt; 

    rho0 = rho;
    rhovx0 = rhovx;
    rhovy0 = rhovy;
    rhovz0 = rhovz;
    eng0 = eng;
    by0 = by;
    bz0 = bz;
    
    [rho_h,rhovx_h,rhovy_h,rhovz_h,eng_h] = ...
        getConservedVariables(rho_h,vx_h,vy_h,vz_h,p_h,gamma);
    
    % step 2: reconstruct the primitive variables to cell faces
    [rho_left, rho_right] = reconstruct_1D_x(rho_h,if_act,PDMB);    
    [vx_left, vx_right] = reconstruct_1D_x(vx_h,if_act,PDMB);   
    [vy_left, vy_right] = reconstruct_1D_x(vy_h,if_act,PDMB);       
    [vz_left, vz_right] = reconstruct_1D_x(vz_h,if_act,PDMB);        
    [p_left, p_right] = reconstruct_1D_x(p_h,if_act,PDMB);  
    [bx_left, bx_right] = reconstruct_1D_x(bx_h,if_act,PDMB);         
    [by_left, by_right] = reconstruct_1D_x(by_h,if_act,PDMB);       
    [bz_left, bz_right] = reconstruct_1D_x(bz_h,if_act,PDMB);      
    
    % step 3: calculate the hydro/magnetic flux/stress through cell faces
    [Frho_p,FrhoVx_p,FrhoVy_p,FrhoVz_p,Feng_p, ...
     ~,~,~,~,~] = getHydroFlux(rho_left,vx_left,vy_left,vz_left,p_left,gamma,1);
    [~,~,~,~,~, ...
     Frho_n,FrhoVx_n,FrhoVy_n,FrhoVz_n,Feng_n] ...
                = getHydroFlux(rho_right,vx_right,vy_right,vz_right,p_right,gamma,1); 
  
    rho_flux_x = Frho_p + Frho_n;
    vx_flux_x = FrhoVx_p + FrhoVx_n;
    vy_flux_x = FrhoVy_p + FrhoVy_n;
    vz_flux_x = FrhoVz_p + FrhoVz_n;
    eng_flux_x = Feng_p + Feng_n;

% add hogs diffusion
% compute alfven velocity
    Va_eff = sqrt((bx.^2+by.^2+bz.^2)./(rho));
    dif_factor = 0.5; 
    diff_va = 0*Va_eff;
    diff_va(if_act) = dif_factor*(Va_eff(if_act-1) + Va_eff(if_act))/2;
    rho_flux_x(if_act) = rho_flux_x(if_act) - diff_va(if_act).*(rho_right(if_act) - rho_left(if_act));
    vx_flux_x(if_act)  = vx_flux_x(if_act)  - diff_va(if_act).*(rho_right(if_act).*vx_right(if_act) - rho_left(if_act).*vx_left(if_act));
    vy_flux_x(if_act)  = vy_flux_x(if_act)  - diff_va(if_act).*(rho_right(if_act).*vy_right(if_act) - rho_left(if_act).*vy_left(if_act));
    vz_flux_x(if_act)  = vz_flux_x(if_act)  - diff_va(if_act).*(rho_right(if_act).*vz_right(if_act) - rho_left(if_act).*vz_left(if_act));
    eng_flux_x(if_act) = eng_flux_x(if_act) - diff_va(if_act).*(0.5.*rho_right(if_act).*(vx_right(if_act).^2+vy_right(if_act).^2+vz_right(if_act).^2) + p_right(if_act)./(gamma-1) ...
                                                               - 0.5.*rho_left(if_act).*(vx_left(if_act).^2+vy_left(if_act).^2+vz_left(if_act).^2) - p_left(if_act)./(gamma-1));
    % calculate magnetic stress/Lorentz force
    [Bstress_x_p, Bstress_y_p, Bstress_z_p, ...
     ~, ~, ~] = getMagneticStress(rho_left,vx_left,vy_left,vz_left,p_left,bx_left,by_left,bz_left,1);
    [~, ~, ~, ...
     Bstress_x_n, Bstress_y_n, Bstress_z_n] ...
              = getMagneticStress(rho_right,vx_right,vy_right,vz_right,p_right,bx_right,by_right,bz_right,1);
   
    BstressX_x = Bstress_x_p + Bstress_x_n;
    BstressY_x = Bstress_y_p + Bstress_y_n;
    BstressZ_x = Bstress_z_p + Bstress_z_n;    

    % maxwell's equation
    [Fbx_p,Fby_p,Fbz_p, ...
           ~,~,~] = getMagneticFlux(rho_left,vx_left,vy_left,vz_left,p_left,bx_left,by_left,bz_left,gamma,1);
    [~,~,~, ...
     Fbx_n,Fby_n,Fbz_n] = getMagneticFlux(rho_right,vx_right,vy_right,vz_right,p_right,bx_right,by_right,bz_right,gamma,1);
    
    by_flux_x = Fby_p + Fby_n;
    bz_flux_x = Fbz_p + Fbz_n;
 
    % step 4, update conserved fluid variables using finite-volume method, to get the correct rho and p, note
    % that in 1-D, it's equivalant to finite difference
    rho(ic_act) = rho0(ic_act) - dt.*(rho_flux_x(ic_act+1)-rho_flux_x(ic_act))./dx(ic_act);
    rhovx(ic_act) = rhovx0(ic_act) - dt.*(vx_flux_x(ic_act+1)-vx_flux_x(ic_act))./dx(ic_act);
    rhovy(ic_act) = rhovy0(ic_act) - dt.*(vy_flux_x(ic_act+1)-vy_flux_x(ic_act))./dx(ic_act);
    rhovz(ic_act) = rhovz0(ic_act) - dt.*(vz_flux_x(ic_act+1)-vz_flux_x(ic_act))./dx(ic_act);
    eng(ic_act) = eng0(ic_act) - dt.*(eng_flux_x(ic_act+1)-eng_flux_x(ic_act))./dx(ic_act);

    % step 5 get premitive variables from the updated conserved variables                                                   
    vx = rhovx./rho;
    vy = rhovy./rho;
    vz = rhovz./rho;
    p = (eng - 0.5.*rho.*(vx.^2+vy.^2+vz.^2)).*(gamma-1);
    
    % then apply the magnetic stress to get the correct velocity
    rhovx(ic_act) = rhovx(ic_act) - dt./dx(ic_act).*(BstressX_x(ic_act+1)-BstressX_x(ic_act)); 
    rhovy(ic_act) = rhovy(ic_act) - dt./dx(ic_act).*(BstressY_x(ic_act+1)-BstressY_x(ic_act)); 
    rhovz(ic_act) = rhovz(ic_act) - dt./dx(ic_act).*(BstressZ_x(ic_act+1)-BstressZ_x(ic_act)); 
    vx = rhovx./rho;
    vy = rhovy./rho;
    vz = rhovz./rho;    
    
    % the apply the magnetic flux to evlove the B fields
    by(ic_act) = by0(ic_act) - dt./dx(ic_act).*(by_flux_x(ic_act+1)-by_flux_x(ic_act)); 
    bz(ic_act) = bz0(ic_act) - dt./dx(ic_act).*(bz_flux_x(ic_act+1)-bz_flux_x(ic_act)); 
    
    Mach = sqrt(vx.^2+vy.^2+vz.^2)./sqrt(gamma.*p./rho);
    Max_Mach = max(Mach(:));
        
    % step 6: applying boundary conditions - use symmetric bc on each side
    rho(ic_rb) = rho(2*nx+NO+1-ic_rb);
    p(ic_rb) = p(2*nx+NO+1-ic_rb);
    vx(ic_rb) = vx(2*nx+NO+1-ic_rb);
    vy(ic_rb) = vy(2*nx+NO+1-ic_rb);
    vz(ic_rb) = vz(2*nx+NO+1-ic_rb);
    by(ic_rb) = by(2*nx+NO+1-ic_rb);
    bz(ic_rb) = bz(2*nx+NO+1-ic_rb);
    
    rho(ic_lb) =rho(NO+1-ic_lb);
    vx(ic_lb)  = vx(NO+1-ic_lb);
    vy(ic_lb)  = vy(NO+1-ic_lb);
    vz(ic_lb)  = vz(NO+1-ic_lb);
    p(ic_lb)   =  p(NO+1-ic_lb);
    by(ic_lb)  = by(NO+1-ic_lb);
    bz(ic_lb)  = bz(NO+1-ic_lb);
    
    % step7: save the primitive variables and dt for the next AB time steping
    toc;
    
    % plot something..
    if (mod(n,10)==1)
        figure(1)
        plot(xc,rho,'-+');pause(0.01)
    end
    if(Time>0.25)
        break;
    end
end

% figure
% plot(xc,rho,'-ro');hold on