
% Jake Davis
% 02/28/2018

function triGeom = triGeom_func(pnt1,pnt2,P)

% Redfine Stuff
x = P(1); y = P(2); z = P(3);
x1 = pnt1(1); y1 = pnt1(2);
x2 = pnt2(1); y2 = pnt2(2);

% K & P
m  = (y2-y1) / (x2-x1);
e1 = (x-x1)^2 + z^2;
e2 = (x-x2)^2 + z^2;
h1 = (x-x1) * (y-y1);
h2 = (x-x2) * (y-y2);
r1 = sqrt((x-x1)^2 + (y-y1)^2 + z^2);
r2 = sqrt((x-x2)^2 + (y-y2)^2 + z^2);

% Johnson
d = sqrt((x2-x1)^2 + (y2-y1)^2);
C = (x2-x1) / d;    % nu_eta
S = (y2-y1) / d;    % nu_xi

a = (x-x1)*S - (y-y1)*C;
h = z;
g = sqrt(a^2 + h^2);
% thing1 = (x1-x)*C
% thing2 = (y1-y)*S
l1 = (x1-x)*C + (y1-y)*S;
l2 = (x2-x)*C + (y2-y)*S;

s1 = sqrt(l1^2 + g^2);
s2 = sqrt(l2^2 + g^2);
c1 = g^2 + abs(h)*s1;
c2 = g^2 + abs(h)*s2;

% nu_xi  = (y2-y1) / (l2-l1);
nu_xi  = (y2-y1) / (l2-l1);
nu_eta = (x2-x1) / (l2-l1);

nuTest = sqrt(nu_xi^2 + nu_eta^2)
lTest = l2 - l1
disp(d)

% Supersonic Stuff
b = nu_xi^2 - nu_eta^2;
gSup = sqrt(a^2 - b*h^2);
R1 = sqrt((g^2 - l1^2) / b);
R2 = sqrt((g^2 - l2^2) / b);

% Define Output Structure
triGeom.b = b; triGeom.gSup = gSup;
triGeom.R1 = R1; triGeom.R2 = R2;

triGeom.m  = m;  triGeom.d  = d;
triGeom.e1 = e1; triGeom.e2 = e2; 
triGeom.h1 = h1; triGeom.h2 = h2;
triGeom.r1 = r1; triGeom.r2 = r2;
triGeom.a  = a;  triGeom.g  = g; 
triGeom.l1 = l1; triGeom.l2 = l2;
triGeom.c1 = c1; triGeom.c2 = c2; 
triGeom.nu_xi = nu_xi; triGeom.nu_eta = nu_eta;

end

%%

% r1 = sqrt((x-x1)^2 + (y-y1)^2 + z^2);
% r2 = sqrt((x-x2)^2 + (y-y2)^2 + z^2);

%%

% R = ((x-x1)*(y2-y1)-(y-y1)*(x2-x1))/d;
% Q = log((r1+r2+d)/(r1+r2-d));
% J1 = atan((m*e1-h1)/(z*r1)) - atan((m*e2-h2)/(z*r2));
% 
% triGeometry = triGeom_Johnson(pnt1,pnt2,P);
% a = triGeometry.a; g = triGeometry.g; l1 = triGeometry.l1; l2 = triGeometry.l2;
% c1 = triGeometry.c1; c2 = triGeometry.c2;
% 
% J = atan2(a*(l2*c1-l1*c2) , c1*c2+a^2*l1*l2);
