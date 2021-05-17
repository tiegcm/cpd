clear all;
%close all;

%
% Input parameters
%

% Load option, 1 for manual with fmax, 2 using inverse error function
% without fmax, 3 using inverse error function with fmax
i_load_option = 3;

% nj= the number of particles (must be even).

nj = 10000;

% nSample= the number of sample fraction values used for creating a 
%  Maxwellian velocity distribution function.

nSample = 500;

% vth= thermal velocity in units of L*omega_e, where L is the length 
%  of the system.

vth = .01;
      
% fvMax.  vth*fvMax is the largest velocity represented for the sample
%  array fSample.  Since the Maxwellian velocity distribution
%  falls off fairly rapidly, it is not necessary for this to be very 
%  large.  Somewhere between 3 and 4 is fine.

fvMax = 4;

% i_quiet_start = 1 for quiet start

i_quiet_start = 0;

% vpert= amplitude of initial perturbation of particle velocities.

vpert = .0;

% px= particle position at integer time step times.
% pv= particle velocity at first for t=0 and later for half time 
%  step times.

px = zeros(nj,1);
pv = zeros(nj,1);

% fSample= array of sample fractions used for creating a Maxwellian
%  velocity distribution function.

rSample = zeros(nSample,1);

% Set up array fSample with probabilities associated with certain 
%  velocities.

[dv,fSample] = setmaxw(fvMax,nSample);
vSample = (0:dv:fvMax)';
figure(1)
plot(vSample,fSample)

% Define initial particle positions (uniform distribution in x with
%  sinusoidal perturbation).  Here pi2=2*pi in radians.

pi2 = 2*pi;

njd2= nj/2;
dpx= 1./njd2;

% fractMax = largest fraction allowed for simple thermal loader
fractMax = 1 - 1/nj;

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

  if( i_load_option==1 )
    vnorm = getmaxw(fract,fSample,nSample,dv);
  elseif( i_load_option==2 )
    vnorm = erfinv( fract );
  else
    vnorm = erfinv( fract*fractMax );   % largest arg is now fractMax 
  end
  vthpart= vnorm*vth;
  vpertpart= vpert*sin( pi2*px(j) );
  pv(j)= vthpart + vpertpart;

% Set up particles with negative velocities.

  pv(j+njd2)= -vthpart + vpertpart;

end;

figure(2)
hist(pv,31)


