
% Jake Davis
% 04/07/2018

function triGeom = supTriGeom_func(pnt1,pnt2,P,delMag)

% Redfine Stuff
x = P(1); y = P(2); z = P(3);
x1 = pnt1(1); y1 = pnt1(2);
x2 = pnt2(1); y2 = pnt2(2);

% Initialize
xm   = 0; xmc  = 0;
ym1  = 0; ym2  = 0;
ym1c = 0; ym2c = 0;
lam  = 0;

% Get influence coefficient variables
m  = (y2-y1) / (x2-x1);
s1 = y - y1;
s2 = y - y2;
r1 = sqrt(s1^2 + z^2);
r2 = sqrt(s2^2 + z^2);

if abs(m) < 10^-delMag
% panel parallel to freestream
    m    = 0;
    xmc  = -s1;
    ym1c = -(x - x1);
    ym2c = -(x - x2);
    sm1  = x - x1;
    sm2  = x - x2;
elseif abs(m) > 10^delMag
% panel perpendicular to freestream
    lam = 0;
    xm  = x - x1 - (y-y1)*lam;
    sm1 = xm + s1*lam;
    sm2 = xm + s2*lam;
    ym1 = s1 - sm1*lam;
    ym2 = s2 - sm2*lam;
else
    lam  = 1 / m;
    xm   = x - x1 - (y-y1)/m;
    sm1  = xm + s1/m;
    sm2  = xm + s2/m;
    ym1  = s1 - sm1/m;
    ym2  = s2 - sm2/m;
    xmc  = m*xm;
    ym1c = m*ym1;
    ym2c = m*ym2;
end

R1 = sqrt(sm1^2 - r1^2);
R2 = sqrt(sm2^2 - r2^2);

% Define Output Structure
triGeom.xm   = xm;   triGeom.xmc  = xmc;
triGeom.ym1  = ym1;  triGeom.ym2  = ym2;
triGeom.ym1c = ym1c; triGeom.ym2c = ym2c;
triGeom.R1   = R1;   triGeom.R2   = R2;
triGeom.m    = m;    triGeom.lam  = lam;

end

%% Subsonic

% % Redfine Stuff
% x = P(1); y = P(2); z = P(3);
% x1 = pnt1(1); y1 = pnt1(2);
% x2 = pnt2(1); y2 = pnt2(2);
% 
% % K & P
% m  = (y2-y1) / (x2-x1);
% e1 = (x-x1)^2 + z^2;
% e2 = (x-x2)^2 + z^2;
% h1 = (x-x1) * (y-y1);
% h2 = (x-x2) * (y-y2);
% r1 = sqrt((x-x1)^2 + (y-y1)^2 + z^2);
% r2 = sqrt((x-x2)^2 + (y-y2)^2 + z^2);
% 
% % Johnson
% d = sqrt((x2-x1)^2 + (y2-y1)^2);
% C = (x2-x1) / d;
% S = (y2-y1) / d;
% 
% a = (x-x1)*S - (y-y1)*C;
% h = z;
% g = sqrt(a^2 + h^2);
% 
% l1 = (x1-x)*C + (y1-y)*S;
% l2 = (x2-x)*C + (y2-y)*S;
% 
% s1 = sqrt(l1^2 + g^2);
% s2 = sqrt(l2^2 + g^2);
% c1 = g^2 + abs(h)*s1;
% c2 = g^2 + abs(h)*s2;
% 
% nu_xi  = -(y2-y1) / (l2-l1);
% nu_eta = (x2-x1) / (l2-l1);
% 
% % Define Output Structure
% triGeom.m  = m;  triGeom.d  = d;
% triGeom.e1 = e1; triGeom.e2 = e2; 
% triGeom.h1 = h1; triGeom.h2 = h2;
% triGeom.r1 = r1; triGeom.r2 = r2;
% triGeom.a  = a;  triGeom.g  = g; 
% triGeom.l1 = l1; triGeom.l2 = l2;
% triGeom.c1 = c1; triGeom.c2 = c2; 
% triGeom.nu_xi = nu_xi; triGeom.nu_eta = nu_eta;

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
