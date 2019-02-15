function [ P0_P, T0_T ] = P_T_expansion( M,gam )

P0_P = (1+M^2*((gam-1)/2))^(gam/(gam-1));
T0_T = 1 + M^2*(gam-1)/2;

end

