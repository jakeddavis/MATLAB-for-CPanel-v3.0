
% Jake Davis
% 12/17/2018

% function [ nu1,nu2,M2,P2_P1,T2_T1 ] = Expansion_Solver( M1,theta2,gam )
function [ M2,P2_P1 ] = Expansion_Solver( M1,theta2,gam )

nu1 = nu_from_M( M1,gam );
nu2 = theta2 + nu1;
M2  = Mach_from_nu_CPG( nu2 );
P2_P1 = (1/P_T_expansion( M2,gam ))*P_T_expansion( M1,gam );

end

