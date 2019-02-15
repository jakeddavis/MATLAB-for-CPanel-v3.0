
% Jake Davis
% 12/17/2018

% function [ beta2,Mn1,Mn2,M2,P2_P1,T2_T1,Cp ] = Oblique_Solver( M1,theta2,gam )
function [ M2,P2_P1] = Oblique_Solver( M1,theta2,gam )

beta2 = beta_from_theta_M( theta2,M1,1,gam );
Mn1   = M1*sind( beta2 );
Mn2   = M2_from_M1( Mn1,gam );
M2    = Mn2/sind(beta2-theta2);
[P2_P1,T2_T1] = P_T_oblique( Mn1,gam );

Cp = (2 / (gam*M1^2)) * (P2_P1 - 1);

end

