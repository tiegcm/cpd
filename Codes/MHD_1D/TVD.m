function [f_left,f_right]=TVD(fm1,f0,f1,f2)

%%

% fm1=1;
% f0=1;
% f1=10;
% f2=10;

% left value
eps=1e-10;
r_i = (f0-fm1)./(f1-f0+eps);
phi_i = (r_i + abs(r_i))./(1+abs(r_i)); % Van Leer slope limiter
% phi_i = max(0,max(min(2.*r_i,1),min(r_i,2))); % Superbee slope limiter

f = 0.5*(f0+f1);
f_left = f0 - phi_i.*(f0-f);

% right value
r_i = (f1-f2)./(f0-f1+eps);
phi_i = (r_i + abs(r_i))./(1+abs(r_i));
% phi_i = max(0,max(min(2.*r_i,1),min(r_i,2)));

f = 0.5*(f1+f0);
f_right = f1 - phi_i.*(f1-f);

end
