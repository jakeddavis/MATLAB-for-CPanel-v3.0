
% Jake Davis
% 02/03/2018

% Near field calculations for source influence

function Phi = constantDoubletNearField(mu,verts,P)

% influenced point
x = P(1); y = P(2); z = P(3);

[m12,m23,m31] = mFunc(verts);
[r1,r2,r3]    = rFunc(verts,x,y,z);
[e1,e2,e3]    = eFunc(verts,x,z);
[h1,h2,h3]    = hFunc(verts,x,y);

Phi = (mu/(4*pi)) * ...
    (atan((m12*e1-h1)/(z*r1)) - atan((m12*e2-h2)/(z*r2)) +...
    atan((m23*e2-h2)/(z*r2)) - atan((m23*e3-h3)/(z*r3)) +...
    atan((m31*e3-h3)/(z*r3)) - atan((m31*e1-h1)/(z*r1)));

end