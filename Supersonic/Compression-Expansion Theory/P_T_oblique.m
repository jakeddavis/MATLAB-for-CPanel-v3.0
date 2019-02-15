function [ P2_P1,T2_T1 ] = P_T_oblique( M,gam )

T2_T1 = (1 + ((2*gam)/(gam+1))*(M^2-1))*(2 + (gam-1)*M^2)/((gam+1)*M^2);
P2_P1 = 1 + (M^2-1)*(2*gam)/(gam+1);

end

