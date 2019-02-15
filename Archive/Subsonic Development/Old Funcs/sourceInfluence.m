
% Jake Davis
% 02/03/2017

% Velocity Potential from Point Source Influence (Far Field Influence)

function Phi = sourceInfluence(sigma,verts,P,ffCond)

% compute geometric parameters of the triangle
[A,Q,diam] = triParams(verts);

% P: influenced point (in flow field) %
x = P(1); y = P(2); z = P(3);
% Q: influencing point (centroid of triangle) %
x0 = Q(1); y0 = Q(2);

% compute hyperbolic distance
r = sqrt((x-x0)^2 + (y-y0)^2 + z^2);

% check for far field condition
if r >= diam*ffCond
    % far field condition met
    Phi = -(sigma*A) / (4*pi*r);
else
% far field condition not met
%     Phi = 'poop';
    Phi = sourceNearField(sigma,verts,P);
end

disp(Phi)

end