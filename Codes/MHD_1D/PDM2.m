function [f_left,f_right]=PDM2(f0,f1,f2,f3,f,PDMB)

maxf = max(f1,f2);
minf = min(f1,f2);

f = max(minf,min(f,maxf));

df0 = PDMB.*(f1-f0);
df1 = PDMB.*(f2-f1);
df2 = PDMB.*(f3-f2);

s0 = sign(df0);
s1 = sign(df1);
s2 = sign(df2);

df0 = abs(df0);
df1 = abs(df1);
df2 = abs(df2);

q0 = abs(s0+s1);
q1 = abs(s1+s2);

df_left = f - f1;
df_righ = f2 - f;

f_left = f - s1.*max(0,abs(df_left) - q0.*df0);
f_right= f + s1.*max(0,abs(df_righ) - q1.*df2); 

end
