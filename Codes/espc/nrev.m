%**********************************************************************
   function  fract =  nrev (i,n)
%**********************************************************************
%
%     Returns a fraction between 0 and 1 representing the
%     n-reversed number corresponding to i.
%     i     :  index to be n-reversed
%     fract := number between 0 and 1 representing the
%              integer i n-reversed.
%     n     :  the integer i will be reversed in base n.

i;
n;
j = i;
powern = 1.;
fract = 0.;

jnext = floor(j/n);
powern = powern/n;
fract = fract + (j-jnext*n)*powern;

while (jnext~=0.0)
  j = jnext;
  jnext = floor(j/n);
  powern = powern/n;
  fract = fract + (j-jnext*n)*powern;
end
