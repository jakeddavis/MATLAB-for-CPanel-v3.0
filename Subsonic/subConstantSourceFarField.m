
% Jake Davis
% 02/03/2017

function Phi = subConstantSourceFarField(sigma,A,Q,P)

x0 = Q(1); y0 = Q(2);
x = P(1);  y = P(2);  z = P(3);

% r = sqrt((x-x0)^2 + (y-y0)^2 + z^2);
r = norm(P - Q);

Phi = (sigma*A) / (4*pi*r);

end