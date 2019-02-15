
verts = [0 0 ; 1 0 ; cosd(60) sind(60)];
[A,Q,diam] = triParams(verts);
P = [0 0 1];

mu0 = 1;
x = P(1); y = P(2); z = P(3);
xi2 = verts(3,1); eta2 = verts(3,2);

% syms mu0 x y z xi xi2 eta eta2
eta_a = @(xi) xi*eta2/xi2;

f = @(eta,xi)(-mu0*z/4/pi) * ((x-xi)^2+(y-eta)^2+z^2)^(-3/2);

integral2(f,0,1,0,xi2)

% int(int(f,eta,0,eta_a),xi,0,xi2)