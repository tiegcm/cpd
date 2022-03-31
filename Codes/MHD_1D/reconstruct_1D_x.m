function [rho_left, rho_right] = reconstruct_1D_x(rho_h,if_act,PDMB)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates the two states on each side of a cell interface,
% using varies types and orders of reconstruction methods
% INPUT: primitive variable, index array of active cell faces, PDMB number
% OUTPUT: left and right states on active cell faces
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%~~~~~ using low-order method without any limiters ~~~~~~~~~~~~~~~~~~~~~~~~
% 1st-order "reconstruction" - there's actually no reconstruction if using
% 1st order since left and right states are just the values on the left and
% right cells of the interface, this would garrantee TVD but over diffusive
%     rho_left(if_act) = rho_h(if_act-1);
%     rho_right(if_act) = rho_h(if_act);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~ using high-order reconstruction with the PDM limiter~~~~~~~~~~~~~~~~
% STEP 1: reconstruction
% 2nd-order reconstruction     
%     rho_interp(if_act) = (rho_h(if_act-1)+rho_h(if_act))/2;

% 4th-order reconstruction
%     rho_interp(if_act) = (-1*rho_h(if_act-2)+7*rho_h(if_act-1)+...
%                                          7*rho_h(if_act)-rho_h(if_act+1))/12;

% 6th-order reconstruction
%     PUT IT IN HERE!

% % %8th-order reconstruction - the default choice which is used in LFM
%      rho_interp(if_act) = (-3*rho_h(if_act-4)+29*rho_h(if_act-3)-139*rho_h(if_act-2)+533*rho_h(if_act-1)+...
%                           533*rho_h(if_act)-139*rho_h(if_act+1)+29*rho_h(if_act+2)-3*rho_h(if_act+3))/840;
% % 
% % 
% % % STEP 2: limiting
% % % Splitting the reconstructed value into left and right states at the
% % % interface
% % 
%     [rho_left(if_act),rho_right(if_act)]= PDM2(rho_h(if_act-2),rho_h(if_act-1),...
%                                                rho_h(if_act),rho_h(if_act+1),rho_interp(if_act),PDMB);
 
%~~~~ Using upwind reconstruction with the PDM limiter~~~~~~~~~~~~~~~~~~~~~

% % 7th order upwind reconstruction
% rho_interp(if_act) = -rho_h(if_act-4)*1/140+rho_h(if_act-3)*5/84-rho_h(if_act-2)*101/420+rho_h(if_act-1)*319/420+...
%                       rho_h(if_act)*107/210-rho_h(if_act+1)*19/210+rho_h(if_act+2)/105;
% % limiting - get the left state
% [rho_left(if_act),~]= PDM2(rho_h(if_act-2),rho_h(if_act-1),rho_h(if_act),rho_h(if_act+1),rho_interp(if_act),PDMB);
% 
% rho_interp(if_act) = -rho_h(if_act+3)*1/140+rho_h(if_act+2)*5/84-rho_h(if_act+1)*101/420+rho_h(if_act)*319/420+...
%                       rho_h(if_act-1)*107/210-rho_h(if_act-2)*19/210+rho_h(if_act-3)/105;
% % limiting - get the right state
% [rho_right(if_act),~]= PDM2(rho_h(if_act+1),rho_h(if_act),rho_h(if_act-1),rho_h(if_act-2),rho_interp(if_act),PDMB);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~ using high-order reconestrucion with the WENO type method~~~~~~~~~~~
%
% 5th-order Monotonicity preserving (MP) WENO
% rho_left(if_act) = MPWENO5(rho_h(if_act-3),rho_h(if_act-2),rho_h(if_act-1),rho_h(if_act),rho_h(if_act+1));
% rho_right(if_act)= MPWENO5(rho_h(if_act+2),rho_h(if_act+1),rho_h(if_act),rho_h(if_act-1),rho_h(if_act-2));

% 5th-order original WENO
% rho_left(if_act) = WENO5(rho_h(if_act-3),rho_h(if_act-2),rho_h(if_act-1),rho_h(if_act),rho_h(if_act+1));
% rho_right(if_act)= WENO5(rho_h(if_act+2),rho_h(if_act+1),rho_h(if_act),rho_h(if_act-1),rho_h(if_act-2));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% 2nd order TVD
[rho_left(if_act),rho_right(if_act)]= TVD(rho_h(if_act-2),rho_h(if_act-1),...
                                                rho_h(if_act),rho_h(if_act+1));

% 3rd order PPM - doesn't work yet... no idea why...
% [rho_left(if_act),rho_right(if_act)]=PPM2(rho_h(if_act-3),rho_h(if_act-2),rho_h(if_act-1),rho_h(if_act),rho_h(if_act+1),rho_h(if_act+2));
% 
end