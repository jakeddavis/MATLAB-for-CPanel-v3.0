
% Jake Davis
% 10/25/2018

function [phiD,phiS] = supLinearFarField(muVec,sigma,verts,P,Q,M,nDir)

trans = supCoordTrans(verts,M,0,0,nDir);
P = trans * (P - Q)';
verts(1,:) = trans * (verts(1,:) - Q)';
verts(2,:) = trans * (verts(2,:) - Q)';
verts(3,:) = trans * (verts(3,:) - Q)';
[A,Q,~] = triParams(verts,P,0);

mu_0 = muVec(1); mu_x = muVec(2); mu_y = muVec(3);
x0 = Q(1); y0 = Q(2); z0 = Q(3);
x = P(1);  y = P(2);  z = P(3);
x1 = verts(1,1); y1 = verts(1,2);
x2 = verts(2,1); y2 = verts(2,2);
x3 = verts(3,1); y3 = verts(3,2);

mu1 = mu_0 + mu_x*x1 + mu_y*y1;
mu2 = mu_0 + mu_x*x2 + mu_y*y2;
mu3 = mu_0 + mu_x*x3 + mu_y*y3;
mu = (A/3)*(mu1 + mu2 + mu3);

RBd = ((x-x0)^2 - (y-y0)^2 - (z-z0)^2)^(-3/2);
RBs = ((x-x0)^2 - (y-y0)^2 - (z-z0)^2)^(-1/2);

phiD = (mu/2/pi) * (z-z0)*RBd;
phiS = (sigma/2/pi) *A*RBs;
% phiS = (sigma/2/pi) * A*RBs;

end

%%

% function Phi = subLinearDoubletFarField(muVec,verts,A,Q,P)
% 
% mu_0 = muVec(1); mu_x = muVec(2); mu_y = muVec(3);
% x0 = Q(1); y0 = Q(2);
% x = P(1);  y = P(2);  z = P(3);
% x1 = verts(1,1); y1 = verts(1,2);
% x2 = verts(2,1); y2 = verts(2,2);
% x3 = verts(3,1); y3 = verts(3,2);
% 
% mu1 = mu_0 + mu_x*x1 + mu_y*y1;
% mu2 = mu_0 + mu_x*x2 + mu_y*y2;
% mu3 = mu_0 + mu_x*x3 + mu_y*y3;
% mu = (A/3)*(mu1 + mu2 + mu3);
% 
% r = ((x-x0)^2 + (y-y0)^2 + z^2)^(-3/2);
% 
% Phi = -(mu/4/pi) * (z*r);
% 
% end
