%**********************************************************************
      function [dv,sampleR] = setmaxw (fmax,nsamplev)
%**********************************************************************
%     This routine, modeled after A. B. Langdon's ES1 loader, 
%     stores an integrated maxwellian in the array sampleR,
%     and should be called prior to calls to GETMAXW.
%     The array sampleR should be unchanged between the call to this
%     routine and calls to GETMAXW.
%     vth is assumed to be unity here (this gives velocity normalized
%        to thermal speed
%     sampleR  := contains the integrated maxwellian to be used by 
%                 GETMAXW
%     nsamplev :  dimension of array sampleR; the larger nsamplev is
%                 the more accurate (and slower) GETMAXW will be
%      integer nsamplev
%      real vth, sampleR(nsamplev)

sampleR =zeros(nsamplev,1);

dv = fmax/(nsamplev-1);
sampleR(1) = 0.;

for i=2:nsamplev
vi = (i-1.5)*dv;
fv = exp ( -.5*vi^2 );
sampleR(i) = sampleR(i-1) + fv;
end;

% Normalize the integrated function to 1.

for i = 1:nsamplev
sampleR(i) = sampleR(i)/sampleR(nsamplev);
end;

