function data = WENO5(hm2,hm1,h,hp1,hp2)
%% WENO-5 reconstruction
h0 = 1/3*hm2-7/6*hm1+11/6*h;
h1 = -1/6*hm1+5/6*h+1/3*hp1;
h2 = 1/3*h+5/6*hp1-1/6*hp2;

IS0 = 13/12*(hm2 - 2*hm1 + h).^2 + 1/4*(hm2 - 4*hm1 + 3*h).^2;
IS1 = 13/12*(hm1 - 2*h + hp1).^2 + 1/4*(hm1 - hp1).^2;
IS2 = 13/12*(h - 2*hp1 + hp2).^2 + 1/4*(3*h - 4*hp1 + hp2).^2;
epsilon = 1e-8;

a0 = 1/10./(epsilon+IS0).^2;
a1 = 6/10./(epsilon+IS1).^2;
a2 = 3/10./(epsilon+IS2).^2;


w0 = a0./(a0+a1+a2);
w1 = a1./(a0+a1+a2);
w2 = a2./(a0+a1+a2);

data = w0.*h0 + w1.*h1 + w2.*h2;
% 

%% WENO-Z resonstruction
% 
% a0 = hm2 - 2*hm1 +   h;
% a1 = hm2 - 4*hm1 + 3*h;
% b0 = 13/12.*a0.*a0 + .25.*a1.*a1;
% 
% a0 = hm1 - 2*h + hp1;
% a1 = hm1 - hp1;
% b1 = 13/12.*a0.*a0 + .25.*a1.*a1;
% 
% a0 = h - 2*hp1 + hp2;
% a1 = 3*h - 4*hp1 + hp2;
% b2 = 13/12.*a0.*a0 + .25.*a1.*a1;
% 
% t5 = abs(b0-b2);
% 
% a0 = 1.0.*(1.0 + t5./(b0 + 1.e-40));
% a1 = 6.0.*(1.0 + t5./(b1 + 1.e-40));
% a2 = 3.0.*(1.0 + t5./(b2 + 1.e-40));
% 
% sum_a = 1.0./(a0 + a1 + a2);
% w0 = a0.*sum_a;
% w1 = a1.*sum_a;
% w2 = a2.*sum_a;
% 
% h0 = 1/3*hm2-7/6*hm1+11/6*h;
% h1 = -1/6*hm1+5/6*h+1/3*hp1;
% h2 = 1/3*h+5/6*hp1-1/6*hp2;
% 
% data = w0.*h0 + w1.*h1 + w2.*h2;
end