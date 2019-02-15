
% Jake Davis
% 02/28/2018

function [R,Q,J] = triGeom(pnt1,pnt2,P)

x = P(1); y = P(2); z = P(3);
x1 = pnt1(1); y1 = pnt1(2);
x2 = pnt2(1); y2 = pnt2(2);

d = sqrt((x2-x1)^2 + (y2-y1)^2);
m = (y2-y1) / (x2-x1);
r1 = sqrt((x-x1)^2 + (y-y1)^2 + z^2);
r2 = sqrt((x-x2)^2 + (y-y2)^2 + z^2);
e1 = (x-x1)^2 + z^2;
e2 = (x-x2)^2 + z^2;
h1 = (x-x1) * (y-y1);
h2 = (x-x2) * (y-y2);

R = ((x-x1)*(y2-y1)-(y-y1)*(x2-x1))/d;
Q = log((r1+r2+d)/(r1+r2-d));
J1 = atan((m*e1-h1)/(z*r1)) - atan((m*e2-h2)/(z*r2));

triGeometry = triGeom_Johnson(pnt1,pnt2,P);
a = triGeometry.a; g = triGeometry.g; l1 = triGeometry.l1; l2 = triGeometry.l2;
c1 = triGeometry.c1; c2 = triGeometry.c2;

J = atan2(a*(l2*c1-l1*c2) , c1*c2+a^2*l1*l2);

end