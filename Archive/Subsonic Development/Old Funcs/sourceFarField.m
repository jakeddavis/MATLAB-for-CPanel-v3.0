
% Jake Davis
% 02/03/2017

function Phi = sourceFarField(sigma,verts,P)

[A,Q,~] = triParams(verts);
r = sqrt((P(1)-Q(1))^2 + (P(2)-Q(2))^2 + P(3)^2);

Phi = -(sigma*A) / (4*pi*r);

end