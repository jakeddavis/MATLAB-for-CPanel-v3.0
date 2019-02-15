
% Jake Davis
% 10/22/2018

%%%%%% Need to control direction of panel normal
%%%%%% Need to take negative of nuXi

function [Q1, w0] = infCoeffComp(pnt1,pnt2,P,delMag)

x = P(1); y = P(2); h = P(3);
x1 = pnt1(1); y1 = pnt1(2);
x2 = pnt2(1); y2 = pnt2(2);
delSmall = 10^-delMag;
delLarge = 10^delMag;

% Generic
m  = (y2-y1) / (x2-x1);

%%%%%%%%%%%%%%% Might be a problem doing out of plane calcs

% Subsonic
d = sqrt((x2-x1)^2 + (y2-y1)^2);
nuEta = -(x2-x1) / d;    % nu_eta
nuXi = (y2-y1) / d;    % nu_xi
% nuEta = (x2-x1) / d;    % nu_eta
% nuXi = -(y2-y1) / d;    % nu_xi
l1 = (x1-x)*nuEta + (y1-y)*nuXi;
l2 = (x2-x)*nuEta + (y2-y)*nuXi;
a = (x1-x)*nuXi + (y1-y)*nuEta;
% a = (x-x1)*nuXi - (y-y1)*nuEta;

% Supersonic
s1 = y - y1;
s2 = y - y2;
if abs(m) < delSmall  % panel parallel to freestream
    m    = 0;
    xmc  = -s1;
    ym1c = -(x - x1);
    ym2c = -(x - x2);
    sm1  = x - x1;
    sm2  = x - x2;
elseif abs(m) > delLarge   % panel perpendicular to freestream
    lam = 0;
    xm  = (x-x1) - (y-y1)*lam;
    sm1 = xm + s1*lam;
    sm2 = xm + s2*lam;
    %%%%%%%%%%%%%%%%%%%%%%%% be sure to check this
    ym1 = s1 - sm1*lam;
    ym2 = s2 - sm2*lam;
else
    lam  = 1 / m;
    xm   = (x-x1) - (y-y1)/m;
    sm1  = xm + s1/m;
    sm2  = xm + s2/m;
    ym1  = s1 - sm1/m;
    ym2  = s2 - sm2/m;
%     ym1 = lam*sm1 - s1;
%     ym2 = lam*sm2 - s2;
    xmc  = m*xm;
    ym1c = m*ym1;
    ym2c = m*ym2;
end

% Comparisons
c_m_SB = -nuXi/nuEta
c_m_SS = m

c_l1_SB = l1
c_l1_SS1 = -nuEta*sm1-nuXi*s1
c_l1_SS2 = nuXi*(sm1/m-s1)
c_l2_SB = l2
c_l2_SS1 = -nuEta*sm2-nuXi*s2
c_l2_SS2 = nuXi*(sm2/m-s2)

c_xm_SB = -a/nuXi
c_xm_SS = xm

b = nuXi^2 - nuEta^2
gSup = sqrt(a^2 - b*h^2);
% R1 = sqrt((g^2 - l1^2) / b);
% R2 = sqrt((g^2 - l2^2) / b);

end