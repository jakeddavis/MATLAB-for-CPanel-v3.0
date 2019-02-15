
% Jake Davis
% 12/17/2018

function Cp = computeCp(gam,M1,P2_P1,Pb)

Cp = (2 / (gam*M1^2)) * (P2_P1*Pb - 1);

end