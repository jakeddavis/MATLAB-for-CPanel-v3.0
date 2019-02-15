
% Jake Davis
% 02/24/2018

function [int,a,l1,l2,c1,c2,g,nu_xi,nu_eta] = triEdge_fieldPoint_geom(V1,V2,P)
% Finds the point on a line defined by V1 and V2 where the normal to this
% lines foes through the line Q

x0 = P(1); y0 = P(2); h = P(3);
x1 = V1(1); y1 = V1(2);
x2 = V2(1); y2 = V2(2);

m = (y2-y1)/(x2-x1);
k = y1 - m*x1;
x = (x0 + m*y0 - m*k) / (m^2 + 1);
y = m*x + k;
int = [x y];

a = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)/sqrt((y2-y1)^2 + (x2-x1)^2);

l1 = sqrt((x1-x)^2 + (y1-y)^2);
l2 = sqrt((x2-x)^2 + (y2-y)^2);

g = sqrt(a^2 + h^2);
c1 = g^2 + abs(h)*sqrt(l1^2 + g^2);
c2 = g^2 + abs(h)*sqrt(l2^2 + g^2);

a_nu_xi  = x - x0;
a_nu_eta = y - y0;
% a_nu_xi  = abs(x - x0);
% a_nu_eta = abs(y - y0);

nu_xi  = a_nu_xi / a;
nu_eta = a_nu_eta / a;

end