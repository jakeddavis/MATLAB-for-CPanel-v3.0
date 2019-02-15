
% Jake Davis
% 02/03/2018

% Near field calculations for source influence

function Phi = sourceNearField(sigma,verts,P)

% influenced point
x = P(1); y = P(2); z = P(3);
% tri vertex 1
x1 = verts(1,1); y1 = verts(1,2);
% tri vertex 2
x2 = verts(2,1); y2 = verts(2,2);
% tri vertex 3
x3 = verts(3,1); y3 = verts(3,2);

[d12,d23,d31] = dFunc(verts);
[m12,m23,m31] = mFunc(verts);
[r1,r2,r3]    = rFunc(verts,x,y,z);
[e1,e2,e3]    = eFunc(verts,x,z);
[h1,h2,h3]    = hFunc(verts,x,y);

Phi = -(sigma/(4*pi)) * ...
    (((((x-x1)*(y2-y1)-(y-y1)*(x2-x1))/d12) * log((r1+r2+d12)/(r1+r2-d12)) +...
    (((x-x2)*(y3-y2)-(y-y2)*(x3-x2))/d23) * log((r2+r3+d23)/(r2+r3-d23)) +...
    (((x-x3)*(y1-y3)-(y-y3)*(x1-x3))/d31) * log((r3+r1+d31)/(r3+r1-d31))) - abs(z)*...
    (atan((m12*e1-h1)/(z*r1)) - atan((m12*e2-h2)/(z*r2)) +...
    atan((m23*e2-h2)/(z*r2)) - atan((m23*e3-h3)/(z*r3)) +...
    atan((m31*e3-h3)/(z*r3)) - atan((m31*e1-h1)/(z*r1))));

end