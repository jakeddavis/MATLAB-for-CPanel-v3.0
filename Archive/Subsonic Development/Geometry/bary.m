
% Jake Davis
% 02/04/2018

function mu_out = bary(mu,verts,P)
% function [mu1,mu2,mu3] = bary(verts,mu,P)

x = P(1); y = P(2);
mu_0 = mu(1); mu_x = mu(2); mu_y = mu(3);
x1 = verts(1,1); x2 = verts(2,1); x3 = verts(3,1);
y1 = verts(1,2); y2 = verts(2,2); y3 = verts(3,2);
detT = (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3);

A0 = (x2*y3 - y2*x3)/detT;
Ax = (y2 - y3)/detT;
Ay = (x3 - x2)/detT;
mu1 = A0*mu_0 + Ax*mu_x + Ay*mu_y;
% mu1 = A0*mu_0 + Ax*mu_x*x + Ay*mu_y*y;

B0 = (y1*x3 - x1*y3)/detT;
Bx = (y3 - y1)/detT;
By = (x1 - x3)/detT;
mu2 = B0*mu_0 + Bx*mu_x + By*mu_y;
% mu2 = B0*mu_0 + Bx*mu_x*x + By*mu_y*y;

C0 = 1 - A0 - B0;
Cx = 1 - Ax - Bx;
Cy = 1 - Ay - By;
mu3 = C0*mu_0 + Cx*mu_x + Cy*mu_y;
% mu3 = C0*mu_0 + Cx*mu_x*x + Cy*mu_y*y;

mu_out = mu1 + mu2 + mu3;
% [A,~,~] = triParams(verts);
% mu_out = A*(mu1 + mu2 + mu3)/3;

end