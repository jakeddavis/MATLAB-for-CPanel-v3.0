
% Jake Davis
% 05/24/2018

function [u,v,w] = subConstantDoubletFarFieldVel(mu,A,Q,P)

x0 = Q(1); y0 = Q(2);
x = P(1); y = P(2); z = P(3);

u = ((3*mu*A)/(4*pi)) * (((x-x0)*z) / ((x-x0)^2+(y-y0)^2+z^2)^(5/2));
v = ((3*mu*A)/(4*pi)) * (((y-y0)*z) / ((x-x0)^2+(y-y0)^2+z^2)^(5/2));
w = -(((mu*A)/(4*pi)) * ((x-x0)^2+(y-y0)^2-(2*z^2))) / ((x-x0)^2+(y-y0)^2+z^2)^(5/2);

end