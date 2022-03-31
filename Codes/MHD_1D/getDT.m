function dt = getDT(rho,vx,vy,vz,bx,by,bz,p,gamma,dx,dy,dz,CFL)
%%
Vfluid = sqrt(vx.^2+vy.^2+vz.^2);
    Btotal = sqrt(bx.^2+by.^2+bz.^2);
    Valfvn = Btotal./sqrt(rho);
    Vsound = sqrt(gamma.*p./rho);
    %
    VCFL = Vfluid+sqrt(Valfvn.^2+Vsound.^2);
    dtCFL = CFL ./(VCFL./dx+VCFL./dy+VCFL./dz);
% dx = min(dx,dy);
%     dtCFL = CFL./VCFL./min(dx,dz);
    dt = min(dtCFL(:));
    
end