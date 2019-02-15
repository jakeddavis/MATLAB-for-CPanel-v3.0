
% Jake Davis
% 11/21/2018

%%%%%% Need to control direction of panel normal
%%%%%% Need to take negative of nuXi

function [Q1, w0] = supInfCoeffsJohn(pnt1,pnt2,P,interPnt,delMag,B)

x = P(1); y = P(2); z = P(3);
h = z;
x1 = pnt1(1); y1 = pnt1(2);
x2 = pnt2(1); y2 = pnt2(2);
delSmall = 10^-delMag;
delLarge = 10^delMag;

m  = (y2-y1) / (x2-x1);
s1 = y - y1;
s2 = y - y2;
r1 = sqrt(s1^2 + z^2);
r2 = sqrt(s2^2 + z^2);

if abs(m) < delSmall  % panel parallel to freestream
    mFlag = 1;
    m    = 0;
    xmc  = -s1;
    ym1c = -(x - x1);
    ym2c = -(x - x2);
    sm1  = x - x1;
    sm2  = x - x2;
    xm = xmc;
elseif abs(m) > delLarge   % panel perpendicular to freestream
    lam = 0;
    xm  = (x-x1) - (y-y1)*lam;
    sm1 = xm + s1*lam;
    sm2 = xm + s2*lam;
    ym1 = s1 - sm1*lam;
    ym2 = s2 - sm2*lam;
else
    lam  = 1 / m;
    xm   = (x-x1) - (y-y1)/m;
    sm1  = xm + s1/m;
    sm2  = xm + s2/m;
    ym1  = s1 - sm1/m;
    ym2  = s2 - sm2/m;
    xmc  = m*xm;
    ym1c = m*ym1;
    ym2c = m*ym2;
end

% % Check if vertices are in DoD
% R1 = 0;
% R2 = 0;
% if dot(P' - pnt1, [1 0 0]) >= 0
%     R1 = real(sqrt((x-x1)^2 - (y-y1)^2 - z^2));
% end
% if dot(P' - pnt2, [1 0 0]) >= 0
%     R2 = real(sqrt((x-x2)^2 - (y-y2)^2 - z^2));
% end
% 
% if R1 > delSmall || R2 > delSmall
%     if R1 < delSmall
%         edgeInterPnt = supEdgeCheckTesting_1126(pnt1,pnt2,P,B);
%         pnt1 = edgeInterPnt;
%         %     x1 = edgeInterPnt(1); y1 = edgeInterPnt(2);
%         %             plot3(x1,y1,0,'r*')
%     end
%     if R2 < delSmall
%         edgeInterPnt = supEdgeCheckTesting_1126(pnt1,pnt2,P,B);
%         pnt2 = edgeInterPnt;
%         %     x2 = edgeInterPnt(1); y2 = edgeInterPnt(2);
%         %             plot3(x2,y2,0,'r*')
%     end
% end

edge12 = pnt2 - pnt1;
edge12(3) = [];
nuVecO(1) = edge12(2);
nuVecO(2) = -edge12(1);
if dot(nuVecO, edge12) > 0
    nu = -nuVecO / norm(nuVecO);
else
    nu = nuVecO / norm(nuVecO);
end
nuXi = nu(1);
nuEta = nu(2);
a = dot(nu, [x1-x y1-y]);
l1 = dot(nu, [y1-y x1-x]);
l2 = dot(nu, [y2-y x2-x]);
b = nuXi^2 - nuEta^2;
g = sqrt(a^2 - b*h^2);

R1 = real(sqrt((g^2-l1^2)/b));
R2 = real(sqrt((g^2-l2^2)/b));

hH113 = 0;
F111 = 0;

if R1 > 0 || R2 > 0
    if b > 0
        F1 = (l1*R2 - l2*R1) / g^2;
        F2 = (b*R1*R2 + l1*l2) / g^2;
        F111 = -atan2(sqrt(b)*F1, F2) / sqrt(b);
%         F111 = -atan((sqrt(b)*F1) / F2) / sqrt(b);
    elseif b < 0
        F1 = (R2^2 - R1^2) / (l1*R2 + l2*R1);
        F2 = (g^2-l1^2-l2^2) / (b*R1*R2 - l1*l2);
        F111 = -(sign(nuEta)/sqrt(b)) * log((sqrt(b)*R1 + abs(l1)) / (sqrt(b)*R2 + abs(l2)));
    else
        error('uh oh...')
    end
    hH113 = atan2(h*a*F1, R1*R2 + h^2*F2);
%     hH113 = atan((h*a*F1) / (R1*R2 + h^2*F2));

    testFs = F2^2 + b*F1^2
elseif (abs(m) > 1)
    xmInter = (interPnt(1)-x1) - (interPnt(2)-y1)/m;
    if (xmInter > 0) && (sign(ym1) ~= sign(ym2)) && (abs(ym1) > 0.0001 && abs(ym2) > 0.0001)
        hH113 = pi * sign(h*nuXi);
        F111 = pi / sqrt(b);
    end
end

Q1 = -hH113;
w0 = -nuXi * F111;

end

%%

% Subsonic
% d = sqrt((x2-x1)^2 + (y2-y1)^2);
% nuEta = -(x2-x1) / d;    % nu_eta
% nuXi = (y2-y1) / d;    % nu_xi
% nuEta = (x2-x1) / d    % nu_eta
% nuXi = -(y2-y1) / d    % nu_xi
% l1 = (x1-x)*nuEta + (y1-y)*nuXi;
% l2 = (x2-x)*nuEta + (y2-y)*nuXi;
% a = (x1-x)*nuXi + (y1-y)*nuEta;
% b = nuXi^2 - nuEta^2;
% g = sqrt(a^2 - b*h^2);

%%

% x = P(1); y = P(2); z = P(3);
% x1 = pnt1(1); y1 = pnt1(2);
% x2 = pnt2(1); y2 = pnt2(2);
% delSmall = 10^-delMag;
% delLarge = 10^delMag;
% 
% % Generic
% m  = (y2-y1) / (x2-x1);
% 
% %%%%%%%%%%%%%%% Might be a problem doing out of plane calcs
% 
% % Subsonic
% d = sqrt((x2-x1)^2 + (y2-y1)^2);
% nuEta = (x2-x1) / d;    % nu_eta
% nuXi = -(y2-y1) / d;    % nu_xi
% l1 = (x1-x)*nuEta + (y1-y)*nuXi;
% l2 = (x2-x)*nuEta + (y2-y)*nuXi;
% a = (x-x1)*nuXi - (y-y1)*nuEta;
% 
% % Supersonic
% s1 = y - y1;
% s2 = y - y2;
% if abs(m) < delSmall  % panel parallel to freestream
%     m    = 0;
%     xmc  = -s1;
%     ym1c = -(x - x1);
%     ym2c = -(x - x2);
%     sm1  = x - x1;
%     sm2  = x - x2;
% elseif abs(m) > delLarge   % panel perpendicular to freestream
%     lam = 0;
%     xm  = (x-x1) - (y-y1)*lam;
%     sm1 = xm + s1*lam;
%     sm2 = xm + s2*lam;
%     %%%%%%%%%%%%%%%%%%%%%%%% be sure to check this
%     ym1 = s1 - sm1*lam;
%     ym2 = s2 - sm2*lam;
% else
%     lam  = 1 / m;
%     xm   = (x-x1) - (y-y1)/m;
%     sm1  = xm + s1/m;
%     sm2  = xm + s2/m;
%     ym1  = s1 - sm1/m;
%     ym2  = s2 - sm2/m;
% %     ym1 = lam*sm1 - s1;
% %     ym2 = lam*sm2 - s2;
%     xmc  = m*xm;
%     ym1c = m*ym1;
%     ym2c = m*ym2;
% end
% 
% % Comparisons
% c_m_SB = -nuXi/nuEta
% c_m_SS = m
% 
% c_l1_SB = l1
% c_l1_SS1 = -nuEta*sm1-nuXi*s1
% c_l1_SS2 = nuXi*(sm1/m-s1)
% c_l2_SB = l2
% c_l2_SS1 = -nuEta*sm2-nuXi*s2
% c_l2_SS2 = nuXi*(sm2/m-s2)
% 
% c_xm_SB = -a/nuXi
% c_xm_SS = xm
% 
% b = nuXi^2 - nuEta^2
% gSup = sqrt(a^2 - b*h^2);
% % R1 = sqrt((g^2 - l1^2) / b);
% % R2 = sqrt((g^2 - l2^2) / b);

