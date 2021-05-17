clear all;
close all;

%
% Input parameters
%

% nj= the number of particles (must be even).

nj = 10000;

% nid= the number of data grid points.  These are numbered from 
%  2 to nid+1.  The dimension ni=nid+2 holds 2 extra "buffer" points 
%  numbered 1 and ni=nid+2.  These two points wrap around periodically
%  so that rho(1)=rho(ni-1) and rho(ni)=rho(2).

nid = 32;

% nsamplev= the number of sample velocities used for creating a 
%  Maxwellian velocity distribution function.

nsamplev = 500;

% vth= thermal velocity in units of L*omega_e, where L is the length 
%  of the system.

vth = 0.05;
      
% fmax.  vth*fmax is the largest velocity represented by the sample
%  velocity array sampleR.  Since the Maxwellian velocity distribution
%  falls off fairly rapidly, it is not necessary for this to be very 
%  large.  Somewhere between 3 and 4 is fine.

fmax = 4;

% i_quiet_start = 1 for quiet start

i_quiet_start = 0;

% vpert= amplitude of initial perturbation of particle velocities.

vpert = 1.0;
      
% dt= time step.

dt = .05;

% tsim = time of simulation.

tsim = 30;

% tfldout= number of time steps between field outputs

tfldout = .5;

% ipartinc= increment between plotted particles
%   = 1 to plot all points
%   = larger number to speed up code

ipartinc = 5;

% tpause = time to pause (in s) between field outputs

tpause = 0.1;

% i_subtract_thermal_energy = 1 to subtract off thermal part of kinetic
% energy

i_subtract_thermal_energy = 1;

%
% Subsidiary parameters (depending on input paramters above)
%

if ( mod(nj,2)~=0 ) % corrects nj if odd
  nj = nj+1;
end;
ni = nid+2;
nt = round(tsim/dt);
nfldout = round(tfldout/dt);

% px= particle position at integer time step times.
% pv= particle velocity at first for t=0 and later for half time 
%  step times.

px = zeros(nj,1);
pv = zeros(nj,1);

% x= spatial grid position.
% rho= charge density (average=0 because includes ion charge density).
% e= electric field.

x=zeros(ni,1);
rho=zeros(ni,1);
e=zeros(ni,1);

% sampleR= array of sample velocities used for creating a Maxwellian
%  velocity distribution function.

sampleR = zeros(nsamplev,1);

% Define x positions. x=0 and 1 at grid points 1.5 and ni-.5.
% This array is helpful for plotting rho 
%  and e.  Note that x=0 inbetween the first (buffer)
%  grid point and the first (i=2) data point.

dx= 1./nid;
for i= 1:ni
  x(i)= (i-1.5)*dx*0.1;
end

% Set up array sampleR with probabilities associated with certain 
%  velocities.

[dv,sampleR] = setmaxw(fmax,nsamplev);

% Define initial particle positions (uniform distribution in x with
%  sinusoidal perturbation).  Here pi2=2*pi in radians.

pi2 = 2*pi;

njd2= nj/2;
dpx= 1./njd2;

make_two_figures;

for j= 1:njd2
  
  px(j)= (j-.5)*dpx;
  
% Initialize particle velocities.  nrev gets the base 2 bit reversed 
%  number corresponding to j.

  if( i_quiet_start==1 )
    fract = nrev(j,2);
  else
    fract = rand;    % Use random number generator
  end
        
% vnorm is the velocity normalized to the thermal speed.

  vnorm = getmaxw(fract,sampleR,nsamplev,dv);
  vthpart= vnorm*vth;
  vpertpart= 0.
  %vpert*sin( pi2*px(j) );
  px(j) = px(j)*(1.+0.5*sin( pi2*px(j)));
 % vpertpart= vpert*sin( pi*px(j) )^2;
  pv(j)= vthpart + vpertpart;

% Set up particles with negative velocities.

  px(j+njd2)=  px(j);
  pv(j+njd2)= -vthpart + vpertpart;

end;

% Do simulation for nt-1 time steps (nt times). Set initial time t to 0.

rhomax= -1;
emax= -1;
pvmax= -1;
t= 0.
for it= 1:nt

% Evaluate charge density d on the spatial (i) grid.
%  First initialize density to zero.  This is done every timestep

  for i= 1:ni
    rho(i)= 0.;
  end;
  
% The following charge per particle will make the average electron 
%  charge density equal to -1.

  charge= (-1.)*nid/nj;
  
% Accumulate electron charge density to grid.

  for j= 1:nj
          
% ri is where the particle is on the grid.  ri is a real 
%  number, and ii and ip are the integer grid positions to the left and 
%  right of the particle location.

    ri= ( px(j)-x(1) )/dx + 1;
    ii= floor(ri);
    ip= ii + 1;
    
% The particle charge contributes to density at both points ii and ip
%  with weights fii and fip.  This is called linear weighting.  Note 
%  that if ri=ii (particle is at lower ii grid point), fii=1 and 
%  fip=0, while if the particle is at ri=ip, fip=1 and fii=0.

    fip= ri - ii;
    fii= 1. - fip;
    rho(ii)= rho(ii) + fii*charge;
    rho(ip)= rho(ip) + fip*charge;
    
  end;  % end of 'for j='
  
% Fix density overflow into buffer region.

  rho(2)= rho(2) + rho(ni);
  rho(ni-1)= rho(ni-1) + rho(1);

% Fix buffer region assuming periodic boundary conditions.

  rho(1)= rho(ni-1);
  rho(ni)= rho(2);

% Add ion charge density.

  for i= 1:ni
    rho(i)= rho(i) + 1.;
  end;

% Calculate electic field e.  For now, let e(2)= 0; will correct 
%  this later.

  e(2)= 0.;

% Now integrate Gauss's Law to find e elsewhere.  Stepping from 
%  grid point i-1 to i, use the average charge density 
%  .5*( rho(i-1)+rho(i) ).

  for i= 3:ni-1
    e(i)= e(i-1) + dx*.5*( rho(i-1)+rho(i) );
  end;

% Adjust e so that average e = 0.  First calculate average e.

  ave= 0.;
  for i= 2:ni-1
    ave= ave + e(i);
  end;
  ave= ave/nid;
  for i= 2:ni-1
    e(i)= e(i) - ave;
  end;

% Fix buffer region.  e values in the buffer region will make the 
%  calculation of e felt by the particles (for particle acceleration) 
%  easier.

  e(1)= e(ni-1);
  e(ni)= e(2);

% Calculate e field energy een (per grid point).

  een(it)= 0.;
  for i= 2:ni-1
    een(it)= een(it) + e(i)^2;
  end;
  een(it)= .5*een(it)/nid;      % Normalized to get average E energy density.
        
% Accelerate particles.  At the same time, calculate particle kinetic 
%  energy at t.  First initialize particle kinetic energy to zero.

  pen(it)= 0.;
  for j= 1:nj
    ri= ( px(j)-x(1) )/dx + 1;
    ii= floor(ri);
    ip= ii + 1;
    fip= ri - ii;
    fii= 1. - fip;

% Now get e field felt by particle, pe.

    pe= fii*e(ii) + fip*e(ip);

% At the start of the run, pv is defined for 
%  t=0, but afterwards, pv will be defined at this point for t-dt/2.

    if( it==1 )
      pen(it)= pen(it) + pv(j)^2;
      pv(j)= pv(j) - .5*dt*pe;
    else
      pvit=  pv(j) - .5*dt*pe;
      pen(it)= pen(it) + pvit^2;
      pv(j)= pv(j) -    dt*pe;
    end;
  end;  % end of 'for j=' 

% rmass is the mass associated with the electron

  rmass= (+1.) * nid / nj;

% Divide pen by nid just like we did for een.

  pen(it)= .5*rmass*pen(it)/nid;
  if( i_subtract_thermal_energy==1 )
    if( it==1 )
      pen0= pen(1) - .5*.5*vpert^2;
    end
    pen(it)= pen(it) - pen0;
  end
  ten(it)= een(it) + pen(it);

% Step particles one time step.

  for j= 1:nj
    px(j)= px(j) + dt*pv(j);

% If particles cross periodic boundary, correct their position.

    if( px(j)<0. ) 
      px(j)= px(j) + 1.;
    end;
    if( px(j)>=1. ) 
      px(j)= px(j) - 1.;
    end;
  end;
  
  t= t + dt
  tarray(it)= t;
  
  % Make plots.

  if( mod(it,nfldout)==0 )
    
    figure(1)
    subplot(3,1,1)
    plot(x,rho);
    rhomax= max( [abs(rho) ; rhomax] );
    ylim([-rhomax rhomax])
    xlim([ x(1) x(ni)])
    ylabel('\rho_q');
    xlabel('x');
    title(strcat('time=',num2str(t)));
    
    subplot(3,1,2)
    plot(x,e);
    emax= max([abs(e) ; emax]);
    ylim([-emax emax])
    xlim([ x(1) x(ni)])
    ylabel('E');
    xlabel('x');
    
    subplot(3,1,3)
    plot(px(1:ipartinc:nj),pv(1:ipartinc:nj),'.');
    pvmax= max([abs(pv) ; pvmax]);
    ylim([-pvmax pvmax])
    xlim([ x(1) x(ni)])
    ylabel('pv');
    xlabel('px');
    pause(tpause);
    
    figure(2)
    subplot(2,1,2)
    plot(tarray,pen,'r',...
      tarray,een,'b',...
      tarray,ten,'k')
    xlim([0 tsim])
  end;

end;  % end of 'for it='
