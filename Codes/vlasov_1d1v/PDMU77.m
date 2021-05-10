function [f_left]=PDMU77(f0,f1,f2,f3,f4,f5,f6,f_itp,PDMB)

maxf = max(f3,f4);
minf = min(f3,f4);

f = max(minf,min(f_itp,maxf));

df0 = PDMB.*(f3-f2);
df1 = PDMB.*(f4-f3);
% df2 = PDMB.*(f3-f2);

s0 = sign(df0);
s1 = sign(df1);
% s2 = sign(df2);

df0 = abs(df0);
% df1 = abs(df1);
% df2 = abs(df2);

q0 = abs(s0+s1);
% q1 = abs(s1+s2);

df_left = f - f3;
% df_righ = f2 - f;

f_pdm = f - s1.*max(0,abs(df_left) - q0.*df0);
% f_right= f + s1.*max(0,abs(df_righ) - q1.*df2); 

    % the compute the deltas
    D1 = f1 - f0;
    D2 = f2 - f1;
    D3 = f3 - f2;
    D4 = f4 - f3;
    D5 = f5 - f4;
    D6 = f6 - f5;
% 
    SmL=(D1>0)&(D2>0)&(D3<0)&(D4<0)&(abs(D2)<abs(D1))&(abs(D3)<abs(D4));
    SmC=(D2>0)&(D3>0)&(D4<0)&(D5<0)&(abs(D3)<abs(D2))&(abs(D4)<abs(D5));
    SmR=(D3>0)&(D4>0)&(D5<0)&(D6<0)&(abs(D4)<abs(D3))&(abs(D5)<abs(D6));
    
    isSmooth = sign(SmL+SmC+SmR);
    f_left = f_pdm.*(1-isSmooth) + isSmooth.*f_itp;
    
end
