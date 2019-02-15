
% Jake Davis
% 02/03/2017

function Phi = constantDoubletFarField(mu,verts,P)

[A,Q,~] = triParams(verts);
r = ((P(1)-Q(1))^2 + (P(2)-Q(2))^2 + P(3)^2)^(-3/2);

Phi = -((mu*A)/(4*pi)) * (P(3)*r);

end