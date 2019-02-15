
% Jake Davis
% 01/04/2019

function [ M2,Cp,ds ] = wedgeSolver( M1,theta2 )

gam = 1.4;

beta2 = beta_from_theta_M( theta2,M1,1,gam );
Mn1   = M1*sind( beta2 );
Mn2   = M2_from_M1( Mn1,gam );
M2    = Mn2/sind(beta2-theta2);
[P2_P1,T2_T1] = P_T_oblique( Mn1,gam );

R = 287; % J/kg*K
SpH_Cp = 1004.5; % J/kg*K

Cp = (2 / (gam*M1^2)) * (P2_P1 - 1);
ds = SpH_Cp*log(T2_T1) - R*log(P2_P1);

end

