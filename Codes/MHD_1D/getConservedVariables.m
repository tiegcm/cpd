function [rho,rhovx,rhovy,rhovz,eng] =...
    getConservedVariables(rho,vx,vy,vz,p,gamma)

rho = rho;
rhovx = rho.*vx;
rhovy = rho.*vy;
rhovz = rho.*vz;
eng= 0.5.*rho.*(vx.^2+vy.^2+vz.^2)+p./(gamma-1);

end
