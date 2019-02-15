
function mu = muBary(verts,P,muVec)

syms mu1 mu2 mu3
% verts = [0 0 ; 1 0 ; cosd(60) sind(60)];
mu_0 = muVec(1); mu_x = muVec(2); mu_y = muVec(3);
x = P(1); y = P(2);
x1 = verts(1,1); x2 = verts(2,1); x3 = verts(3,1);
y1 = verts(1,2); y2 = verts(2,2); y3 = verts(3,2);
detT = (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3);
% x = .1; y = 2;

A0 = (x2*y3 - y2*x3)/detT;
Ax = (y2 - y3)/detT;
Ay = (x3 - x2)/detT;
lam1 = A0 + Ax*x + Ay*y;

B0 = (y1*x3 - x1*y3)/detT;
Bx = (y3 - y1)/detT;
By = (x1 - x3)/detT;
lam2 = B0 + Bx*x + By*y;

C0 = 1 - A0 - B0;
Cx = 1 - Ax - Bx;
Cy = 1 - Ay - By;
lam3 = C0 + Cx*x + Cy*y;

% mu_0 = 1; mu_x = 0; mu_y = 0;
mu_xy = mu_0 + mu_x*x + mu_y*y;

f1 = (mu_xy - mu2*lam2 - mu3*lam3)/lam1 == mu1;
f2 = (mu_xy - mu1*lam1 - mu3*lam3)/lam2 == mu2;
f3 = (mu_xy - mu1*lam1 - mu2*lam2)/lam3 == mu3;

[mu_1,mu_2,mu_3] = solve(f1,f2,f3);
A = abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))/2;

% mu = A*(mu_1 + mu_2 + mu_3)/3;
mu = mu_1 + mu_2 + mu_3;

