%**********************************************************************
      function v = getmaxw (fract,sampleR,nsamplev,dv)
%**********************************************************************
%     Returns the value of v = vth*erf**(-1)(fract).
%     v         := value of inverted, integrated maxwellian of thermal
%                  velocity vth, evaluated at value fract.
%     fract     :  value at which the inverted, integrated maxwellian
%                  is to be evaluated -must- be between 0 and 1 
%                  inclusive.
%     sampleR   :  sampleR array initialized by SETMAXW.
%               := unchanged on exit.
%     nsamplev  :  dimension of array sampleR supplied to SETMAXW.
%               := unchanged on exit.

if (fract<0.)|(fract>1.)
  printf('\nGETMAXW : fract not between 0 and 1\n');
  return; 
end;
nl = 1;
nu = nsamplev;

ntrial = floor((nl + nu)/2);

while (ntrial~=nl)
  if (sampleR((ntrial))<fract) 
    nl = ntrial;
  end;

  if (sampleR((ntrial))>=fract)
    nu = ntrial;
  end;

ntrial = floor((nl + nu)/2);

end;

v = dv*( nl-1 + (fract-sampleR(floor(nl))) / (sampleR(floor(nu))-sampleR(floor(nl))) );




