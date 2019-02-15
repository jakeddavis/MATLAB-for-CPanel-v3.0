
% Jake Davis
% 09/08/2018

function triGeom = supTriGeomTest(pnt1,pnt2,P)

% Redfine Stuff
x = P(1); y = P(2); z = P(3);
h = z;
x1 = pnt1(1); y1 = pnt1(2);
x2 = pnt2(1); y2 = pnt2(2);

% Original Stuff
m    = (y2-y1) / (x2-x1);
lam  = 1 / m;
s1   = y - y1;
s2   = y - y2;
r1   = sqrt(s1^2 + z^2);
r2   = sqrt(s2^2 + z^2);
xm   = (x-x1) - (y-y1)/m;
xm2 = x - x1 - (y-y1)/m;
sm1  = xm + s1/m;
sm2  = xm + s2/m;
ym1  = lam*sm1 - s1;
ym2  = lam*sm2 - s2;
xmc  = m*xm;
ym1c = sm1 - m*s1;
ym2c = sm2 - m*s2;
R1   = sqrt(sm1^2 - r1^2);
R2   = sqrt(sm2^2 - r2^2);

% New Stuff
d = sqrt((x2-x1)^2 + (y2-y1)^2);
nu_eta = (x2-x1) / d;    % C
nu_xi = -(y2-y1) / d;    % S
% m2 = -nu_xi/nu_eta;
a = (x-x1)*nu_xi - (y-y1)*nu_eta;
% a2 = (x-x2)*nu_xi - (y-y2)*nu_eta;
l1 = (x1-x)*nu_eta + (y1-y)*nu_xi;
l2 = (x2-x)*nu_eta + (y2-y)*nu_xi;

b = nu_xi^2 - nu_eta^2;
g = sqrt(a^2 - b*h^2);
R1 = sqrt((g^2 - l1^2) / b);
R2 = sqrt((g^2 - l2^2) / b);
% r1 = sqrt((x-x1)^2 + (y-y1)^2 + z^2);
% r2 = sqrt((x-x2)^2 + (y-y2)^2 + z^2);

% Testing

% l1Test1 = -nu_eta*sm1 - nu_xi*s1;   % check
% l1Test2 = nu_xi*(sm1/m - s1);       % check
% l2Test1 = -nu_eta*sm2 - nu_xi*s2;   % check
% l2Test2 = nu_xi*(sm2/m - s2);       % check
% lamTest = -nu_eta/nu_xi;            % check
% mTest = -nu_xi/nu_eta;              % check
% xmTest1 = -a/nu_eta;
% xmTest2 = (a*nu_eta-d*nu_xi)/b;
% xmTest3 = (a*nu_eta-l1*nu_xi)/b;
xmcTest1 = -a/nu_eta;


if b > 0    % super
    F1 = (l1*R2 - l2*R1)/g^2;
    F2 = (b*R1*R2 + l1*l2) / g^2;
elseif b < 0    % sub
    F1 = (R2^2 - R1^2) / (l1*R2 + l2*R1);
    F2 = (g^2 - l1^2 - l2^2) / (b*R1*R2 - l1*l2);
else
    error('Uh oh...')
end

% Define Output Structure

% old stuff
triGeom.m  = m;
triGeom.a  = a;  triGeom.g  = g; 
triGeom.l1 = l1; triGeom.l2 = l2;
triGeom.nu_xi = nu_xi; triGeom.nu_eta = nu_eta;

% new stuff
triGeom.b = b; triGeom.g = g;
triGeom.R1 = R1; triGeom.R2 = R2;
triGeom.F1 = F1; triGeom.F2 = F2;

end