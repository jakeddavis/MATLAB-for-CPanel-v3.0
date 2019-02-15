
% Jake Davis
% 10/10/2018

function [Q1, w0] = supEdgeInfSon_preDec13(ym1,ym2,xmc,ym1c,ym2c,R1,R2,lam,z)

ym1 = -ym1;
ym2 = -ym2;
ym1c = -ym1c;
ym2c = -ym2c;

zr = (ym1*R2 - ym2*R1) / (ym1*ym2 + (1-lam^2)*R1*R2);
w0 = zr * (1 - ((1-lam^2)*zr^2)/3 + ((1-lam^2)^2*zr^4)/5 ...
    - ((1-lam^2)^3*zr^6)/7); % + . . . (atan2 expansion)
Q1 = atan2(z*xmc * (ym1c*R2-ym2c*R1) , z^2*ym1c*ym2c+xmc^2*R1*R2);

end