function data = MPWENO5(fm2,fm1,f,fp1,fp2)

f0 = 1/3*fm2-7/6*fm1+11/6*f;
f1 = -1/6*fm1+5/6*f+1/3*fp1;
f2 = 1/3*f+5/6*fp1-1/6*fp2;

a = fm1 - fm2;
b = f - fm1;
c = fp1 - f;
d = fp2 - fp1;

IS0 = 13*(a-b).^2+3*(a-3*b).^2;
IS1 = 13*(b-c).^2+3*(b+c).^2;
IS2 = 13*(c-d).^2+3*(3*c-d).^2;

epsilon = 1e-10;

a0 = 1./(epsilon+IS0).^2;
a1 = 6./(epsilon+IS1).^2;
a2 = 3./(epsilon+IS2).^2;

w0 = a0./(a0+a1+a2);
w1 = a1./(a0+a1+a2);
w2 = a2./(a0+a1+a2);

data = w0.*f0+w1.*f1+w2.*f2;

% dj = hp1 - 2*h + hm1;
% djp1 = hp2 - 2*hp1 + h;
% dMM = minmod(dj,djp1);
% 
% alp=2;
% beta = 4;
% 
% hUL = h + alp*(h - hm1);
% hMD = 0.5*(h + hp1) - 0.5*dMM;
% hLC = h +0.5*(h-hm1) + beta/3*dMM;
% 
% UL_min = max(min(h,min(hp1,hMD)),min(h,min(hUL,hLC)));
% UL_max = min(max(h,max(hp1,hMD)),max(h,max(hUL,hLC)));
% 
% data1 = median(data,UL_min,UL_max);

alpha = 4.0;
epsm = 1.e-12;

fMP = f + minmod(fp1-f,alpha*(f-fm1)); %%%

  if ((data - f).*(data - fMP) <= epsm) 
      return;
  else

  d2m = fm2 + f - 2.0*fm1;   % /* -- Eq. (2.19) -- */
  d2  = fm1 + fp1 - 2.0*f;
  d2p = f + fp2 - 2.0*fp1;   % /* -- Eq. (2.19) -- */

  scrh1 = minmod(4.0*d2 - d2p, 4.0*d2p - d2);
  scrh2 = minmod(d2, d2p);
  dMMp  = minmod(scrh1,scrh2);  % /* -- Eq. (2.27) -- */

  scrh1 = minmod(4.0*d2m - d2, 4.0*d2 - d2m);
  scrh2 = minmod(d2, d2m);
  dMMm  = minmod(scrh1,scrh2); %   /* -- Eq. (2.27) -- */

  fUL = f + alpha*(f - fm1);  % /* -- Eq. (2.8) -- */
  fAV = 0.5*(f + fp1);        %/* -- Eq. (2.16) -- */
  fMD = fAV - 0.5*dMMp; %/* -- Eq. (2.28) -- */
  fLC = 0.5*(3.0*f - fm1) + 4.0/3.0*dMMm; % /* -- Eq. (2.29) -- */

  scrh1 = min(f, fp1); scrh1 = min(scrh1, fMD);
  scrh2 = min(f, fUL);    scrh2 = min(scrh2, fLC);
  fmin  = max(scrh1, scrh2);  %/* -- Eq. (2.24a) -- */

  scrh1 = max(f, fp1); scrh1 = max(scrh1, fMD);
  scrh2 = max(f, fUL);    scrh2 = max(scrh2, fLC);
  fmax  = min(scrh1, scrh2);  %/* -- Eq. (2.24b) -- */

  data = median(data, fmin, fmax); %/* -- Eq. (2.26) -- */
  end
end