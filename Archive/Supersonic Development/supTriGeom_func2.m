
% Jake Davis
% 09/08/2018

function triGeom = supTriGeom_func2(pnt1,pnt2,P)

% Redfine Stuff
x = P(1); y = P(2); z = P(3);
h = z;
x1 = pnt1(1); y1 = pnt1(2);
x2 = pnt2(1); y2 = pnt2(2);

% old stuff
m  = (y2-y1) / (x2-x1);

d = sqrt((x2-x1)^2 + (y2-y1)^2);
nu_eta = (x2-x1) / d;    % C
nu_xi = -(y2-y1) / d;    % S

m2 = -nu_xi/nu_eta;

a = (x-x1)*nu_xi - (y-y1)*nu_eta;
l1 = (x1-x)*nu_eta + (y1-y)*nu_xi;
l2 = (x2-x)*nu_eta + (y2-y)*nu_xi;

% new stuff
b = nu_xi^2 - nu_eta^2;
g = sqrt(a^2 - b*h^2);
R1 = sqrt((g^2 - l1^2) / b);
R2 = sqrt((g^2 - l2^2) / b);

r1 = sqrt((x-x1)^2 + (y-y1)^2 + z^2);
r2 = sqrt((x-x2)^2 + (y-y2)^2 + z^2);

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

%% 09/08/2018 - Sub Stuff that's not needed

% % K & P
% m  = (y2-y1) / (x2-x1);
% e1 = (x-x1)^2 + z^2;
% e2 = (x-x2)^2 + z^2;
% h1 = (x-x1) * (y-y1);
% h2 = (x-x2) * (y-y2);
% r1 = sqrt((x-x1)^2 + (y-y1)^2 + z^2);
% r2 = sqrt((x-x2)^2 + (y-y2)^2 + z^2);

% % Johnson
% d = sqrt((x2-x1)^2 + (y2-y1)^2);
% C = (x2-x1) / d;    % nu_eta
% S = (y2-y1) / d;    % nu_xi

% a = (x-x1)*S - (y-y1)*C;
% h = z;
% g = sqrt(a^2 + h^2);
% % thing1 = (x1-x)*C
% % thing2 = (y1-y)*S
% l1 = (x1-x)*C + (y1-y)*S;
% l2 = (x2-x)*C + (y2-y)*S;

% s1 = sqrt(l1^2 + g^2);
% s2 = sqrt(l2^2 + g^2);
% c1 = g^2 + abs(h)*s1;
% c2 = g^2 + abs(h)*s2;

% nu_xi  = (y2-y1) / (l2-l1);
% % nu_xi  = -(y2-y1) / (l2-l1);
% nu_eta = (x2-x1) / (l2-l1);

% triGeom.m  = m;  triGeom.d  = d;
% triGeom.e1 = e1; triGeom.e2 = e2; 
% triGeom.h1 = h1; triGeom.h2 = h2;
% triGeom.r1 = r1; triGeom.r2 = r2;
% triGeom.a  = a;  triGeom.g  = g; 
% triGeom.l1 = l1; triGeom.l2 = l2;
% triGeom.c1 = c1; triGeom.c2 = c2; 
% triGeom.nu_xi = nu_xi; triGeom.nu_eta = nu_eta;

