
% Jake Davis
% 10/09/2018

function [Q1, w0, xm, m] = ...
    supInfCoeffs_10_09(pnt1,pnt2,P,delMag,depFlag,pntD1,pntD2,edgeD)

% Redfine Stuff
x = P(1); y = P(2); z = P(3);
x1 = pnt1(1); y1 = pnt1(2);
x2 = pnt2(1); y2 = pnt2(2);
delSmall = 10^-delMag;
delLarge = 10^delMag;

% Initialize
xm   = 0; xmc  = 0;
ym1  = 0; ym2  = 0;
ym1c = 0; ym2c = 0;
lam  = 0;

Q1 = 0; w0 = 0;

% Compute influence coefficient variables
m  = (y2-y1) / (x2-x1);
s1 = y - y1;
s2 = y - y2;
r1 = sqrt(s1^2 + z^2);
r2 = sqrt(s2^2 + z^2);

if abs(m) < delSmall  % panel parallel to freestream
    m    = 0;
    xmc  = -s1;
    ym1c = -(x - x1);
    ym2c = -(x - x2);
    sm1  = x - x1;
    sm2  = x - x2;
elseif abs(m) > delLarge   % panel perpendicular to freestream
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

% Compute Influence Coefficients

if depFlag == 1 || (pntD1 == 1 && pntD2 == 1)
% Three situations can result in this case
% 1. Sonic edge completely inside of Mach cone
% 2. Subsonic edge completely inside Mach cone
% 3. Supersonic edge completely inside Mach cone
    if abs(z) < delSmall
    % POI is on panel plane
        w0 = 0; % not true here, but mult'pd by z = 0 later
        Q1 = 0;
        if abs(m) < 1
        % subsonic edge
            if R1 > delSmall || R2 > delSmall
                Q1 = -pi * sign(z);
            end
        end
    else
        if abs(1 - abs(m)) < delSmall
        % sonic edge (only a thing for real R1 AND R2)
            [Q1, w0] = supEdgeInfSonic(ym1,ym2,xmc,ym1c,ym2c,R1,R2,lam,z);
        elseif abs(m) < 1
        % subsonic edge
            Q1 = 0;
            w0 = 0;
            if R1 > delSmall && R2 > delSmall
                Q1 = atan2(z*xmc * (ym1c*R2-ym2c*R1) , z^2*ym1c*ym2c + xmc^2*R1*R2);
                w0 = (m/sqrt(1-m^2)) * log((ym2c+R2*sqrt(1-m^2))/(ym1c+R1*sqrt(1-m^2)));
            end
        elseif abs(m) > 1
        % supersonic edge
            if R1 > delSmall && R2 > delSmall
                Q1 = atan2(z*xm * (ym1*R2-ym2*R1) , z^2*ym2*ym1 + xm^2*R1*R2);
                w0 = (1/sqrt(1-lam^2)) * atan2((sqrt(1-lam^2)*(ym1*R2-ym2*R1)) , ...
                    (ym1*ym2+(1-lam^2)*R1*R2));
            elseif R1 < delSmall
                Q1 = sign(z*(s1-lam*sm1)) * pi/2;
                w0 = sign(s1-lam*sm1) * (pi/2)*sqrt(1-lam^2);
            elseif R2 < delSmall
                Q1 = sign(z*(s2-lam*sm2)) * pi/2;
                w0 = sign(s2-lam*sm2) * (pi/2)*sqrt(1-lam^2);
            end
        end
    end
elseif pntD1 == 1 && pntD2 == 0
% Two situations can result in this case
% 1. Pnt 1 of supersonic edge inside Mach cone, Pnt 2 outside
% 2. Pnt 1 of subsonic edge inside Mach cone, Pnt 2 outside
    if abs(z) < delSmall
    % POI is on panel plane
        [Q1,w0] = supEdgeInfInPlane(m,R1,z,ym1,ym2,delSmall);
    else
        [Q1,w0] = supEdgeInfVertex(m,lam,R1,xm,ym1,s1,sm1,z,delSmall);
    end
elseif pntD1 == 0 && pntD2 == 1
% Two situations can result in this case
% 1. Pnt 2 of supersonic edge inside Mach cone, Pnt 1 outside
% 2. Pnt 2 of subsonic edge inside Mach cone, Pnt 1 outside
    if abs(z) < delSmall
    % POI is on panel plane
        [Q1,w0] = supEdgeInfInPlane(m,R2,z,ym1,ym2,delSmall);
    else
        [Q1,w0] = supEdgeInfVertex(m,lam,R2,xm,ym2,s2,sm2,z,delSmall);
    end
elseif edgeD == 1
    if abs(m) < 1
        Q1 = 0;
        w0 = 0;
    elseif abs(m) > 1
        %%%%%%%%%%% test to see if the sign changes whether using 1 or 2
        Q1 = sign(z*(s2-lam*sm2)) * pi/2;
        w0 = sign(s2-lam*sm2) * (pi/2)*sqrt(1-lam^2);
        Q1test = sign(z*(s1-lam*sm1)) * pi/2;
        w0test = sign(s1-lam*sm1) * (pi/2)*sqrt(1-lam^2);
    end
end

end

%% More Stuff

% if abs(z) < delSmall
% % POI is on panel plane
%     w0 = 0; % not true, but mult'pd by z = 0 later
%     if abs(m) < 1
%     % subsonic edge
%         if (isreal(R1) == 1 || isreal(R2) == 1)
%             Q1 = -pi * sign(z);
%         else
%             Q1 = 0;
%         end
%     elseif abs(m) > 1
%     % supersonic edge
%         if (ym1 < 0 && ym2 > 0)
%             Q1 = -pi * sign(z);
%         elseif (sign(ym1) == sign(ym2)) || R < delSmall
%             Q1 = 0;
%         end
%     end
% else
%     if abs(1 - abs(m)) < delSmall
%     % sonic edge
%         zr = (ym1*R2 - ym2*R1) / (ym1*ym2 + (1-lam^2)*R1*R2);
%         w0 = zr * (1 - ((1-lam^2)*zr^2)/3 + ((1-lam^2)^2*zr^4)/5 ...
%             - ((1-lam^2)^3*zr^6)/7); % + . . . (atan2 expansion)
%         if isreal(R1) == 1 && isreal(R2) == 1
%             Q1 = atan2(z*xmc * (ym1c*R2-ym2c*R1) , z^2*ym1c*ym2c+xmc^2*R1*R2);
%         elseif isreal(R1) == 1
%             Q1 = atan((-z*ym1) / (xm*R1));
%         elseif isreal(R2) == 1
%             Q1 = atan((-z*ym2) / (xm*R2));
%         else
%             Q1 = 0;
%         end
%     elseif abs(m) < 1
%     % subsonic edge
%         if isreal(R1) == 1 && R1 < delSmall
%         elseif isreal(R2) == 1 && R2 < delSmall
%         else
%         end
%     elseif abs(m) > 1
%     % supersonic edge
%     end
% end

%% Stuff

% if abs(z) < delSmall
%     Q1 = 0;
%     w0 = 0; % not true, but mult'pd by z = 0 later
% else
%     if abs(1 - abs(m)) < delSmall
%         % sonic edge (use either Q1)
%         zr = (ym1*R2 - ym2*R1) / (ym1*ym2 + (1-lam^2)*R1*R2);
%         w0 = zr * (1 - ((1-lam^2)*zr^2)/3 + ((1-lam^2)^2*zr^4)/5 ...
%             - ((1-lam^2)^3*zr^6)/7); % + . . . (atan2 expansion)
%         Q1 = atan2((z*xm*(ym1*R2-ym2*R1)) , (xm^2*R1*R2+z^2*ym2*ym1));
%     elseif abs(m) < 1
%         % subsonic edge
%         Q1 = atan2(z*xmc * (ym1c*R2-ym2c*R1) , z^2*ym1c*ym2c+xmc^2*R1*R2);
%         w0 = (m/sqrt(1-m^2)) * log((ym2c+R2*sqrt(1-m^2))/(ym1c+R1*sqrt(1-m^2)));
%     elseif abs(m) > 1
%         % supersonic edge
%         Q1 = atan2((z*xm*(ym1*R2-ym2*R1)) , (xm^2*R1*R2+z^2*ym2*ym1));
%         w0 = (1/sqrt(1-lam^2)) * atan2((sqrt(1-lam^2)*(ym1*R2-ym2*R1)) , ...
%             (ym1*ym2+(1-lam^2)*R1*R2));
%     end
% end

% Pure edge intersection
% w0 = sign(s1-lam*sm1) * (pi/2)*sqrt(1-lam^2);
%                 if ym2 > 0 && ym1 < 0
%                     Q1 = -sign(z);
%                 elseif sign(ym1) == sign(ym2)
%                     Q1 = 0;
%                 end


%% Halfway through supersonic edge fix, still wrong

% std1Flag = 0;
% std2Flag = 0;
% if abs(m) < 1       % subsonic edge
%     if abs(R1) < delSmall || abs(R2) < delSmall || abs(z) < delSmall
%     % Field point is close to panel edgepoint Mach cones
%     % Or field point is on panel plane
%         if abs(R1) < delSmall
%             w01 = 0;
%         else
%             w01 = (m/sqrt(1-m^2)) * atan2(m*s1 - sm1, R1*sqrt(1-m^2));
%         end
%         
%         if abs(R2) < delSmall
%             w02 = 0;
%         else
%             w02 = (m/sqrt(1-m^2)) * atan2(m*s2 - sm2, R2*sqrt(1-m^2));
%         end
%         
%         w0 = w01 - w02;
%     else
%         Q1 = atan2(z*xmc * (ym1c*R2-ym2c*R1) , z^2*ym1c*ym2c+xmc^2*R1*R2);
%         w0 = (m/sqrt(1-m^2)) * log((ym2c+R2*sqrt(1-m^2))/(ym1c+R1*sqrt(1-m^2)));
%     end
% elseif abs(m) > 1   % supersonic edge
%     if abs(R1) < delSmall || abs(R2) < delSmall
%     % Field point is close to panel edgepoint Mach cones
%         std2Flag = 1;
%         if abs(R1) < delSmall
% %             w01 = sign(s1-lam*sm1) * (pi/2)*sqrt(1-lam^2);
%             w0 = sign(s1-lam*sm1) * (pi/2)*sqrt(1-lam^2);
% %         else
% % %             w01 = (m/sqrt(1-m^2)) * atan2(m*s1 - sm1, R1*sqrt(1-m^2));
% %             w01 = (1/sqrt(1-lam^2)) * atan2(-ym1, R1*sqrt(1-lam^2));
%         end
%         
%         if abs(R2) < delSmall
% %             w02 = sign(s2-lam*sm2) * (pi/2)*sqrt(1-lam^2);
%             w0 = sign(s2-lam*sm2) * (pi/2)*sqrt(1-lam^2);
% %         else
% % %             w02 = (m/sqrt(1-m^2)) * atan2(m*s2 - sm2, R2*sqrt(1-m^2));
% %             w02 = (1/sqrt(1-lam^2)) * atan2(-ym2, R2*sqrt(1-lam^2));
%         end
% %         w0 = w01 - w02;
%     end
%         
%     if abs(z) < delSmall
%         std2Flag = 1;
%         if ym2 > 0 && ym1 < 0
%             Q1 = -sign(z);
%         elseif sign(ym1) == sign(ym2)
%             Q1 = 0;
%         end
%         
%     elseif abs(R1) < delSmall || abs(R2) < delSmall
%         std2Flag = 1;
%         if abs(R1) < delSmall
% %             Q11 = sign(z*(s1-lam*sm1)) * pi/2;
%             Q1 = sign(z*(s1-lam*sm1)) * pi/2;
% %         else
% %             Q11 = atan2(-z*ym1, xm*R1);
%         end
%         
%         if abs(R2) < delSmall
% %             Q12 = sign(z*(s2-lam*sm2)) * pi/2;
%             Q1 = sign(z*(s2-lam*sm2)) * pi/2;
% %         else
% %             Q12 = atan2(-z*ym2, xm*R2);
%         end
% %         Q1 = Q11 - Q12;
%     end
%         
%     if std2Flag == 0
%         Q1 = atan2((z*xm*(ym1*R2-ym2*R1)) , (xm^2*R1*R2+z^2*ym2*ym1));
%         w0 = (1/sqrt(1-lam^2)) * atan2((sqrt(1-lam^2)*(ym1*R2-ym2*R1)) , ...
%             (ym1*ym2+(1-lam^2)*R1*R2));
%     end
% end
% 
% if abs(1 - abs(m)) < delSmall   % sonic edge (Q1 is fine)
%     zr = (ym1*R2 - ym2*R1) / (ym1*ym2 + (1-lam^2)*R1*R2);
%     w0 = zr * (1 - ((1-lam^2)*zr^2)/3 + ((1-lam^2)^2*zr^4)/5 ...
%         - ((1-lam^2)^3*zr^6)/7); % + . . . (atan2 expansion)
% end

%% Misunderstanding of z = 0 conditions, pretty much just for Mach wedge case

%     if abs(1 - abs(m)) < delSmall   % sonic edge (Q1 is fine)
%         if abs(z) < delSmall
%             Q1 = -pi * sign(z);
%             w0 = 0; % not true, but mult'pd by z = 0 later
%         else
%             Q1 = atan2((z*xm*(ym1*R2-ym2*R1)) , (xm^2*R1*R2+z^2*ym2*ym1));
%             zr = (ym1*R2 - ym2*R1) / (ym1*ym2 + (1-lam^2)*R1*R2);
%             w0 = zr * (1 - ((1-lam^2)*zr^2)/3 + ((1-lam^2)^2*zr^4)/5 ...
%                 - ((1-lam^2)^3*zr^6)/7); % + . . . (atan2 expansion)
%         end
%     else
%         if abs(m) < 1       % subsonic edge
%             if abs(z) < delSmall
%                 Q1 = -pi * sign(z);
%                 w0 = 0; % not true, but mult'pd by z = 0 later
%             end    
%             
%             if abs(z) < delSmall
%                 Q1 = -pi * sign(z);
%             else
%                 Q1 = atan2(z*xmc * (ym1c*R2-ym2c*R1) , z^2*ym1c*ym2c+xmc^2*R1*R2);
%             end
%             w0 = (m/sqrt(1-m^2)) * log((ym2c+R2*sqrt(1-m^2))/(ym1c+R1*sqrt(1-m^2)));
%         elseif abs(m) > 1   % supersonic edge
%             if abs(z) < delSmall
%                 Q1 = 0;
%             else
%                 Q1 = atan2((z*xm*(ym1*R2-ym2*R1)) , (xm^2*R1*R2+z^2*ym2*ym1));
%             end
%             w0 = (1/sqrt(1-lam^2)) * atan2((sqrt(1-lam^2)*(ym1*R2-ym2*R1)) , ...
%                     (ym1*ym2+(1-lam^2)*R1*R2));
%         end
%     end

%% 09/13/2018 - I think these only apply for z = 0

% % Compute Influence Coefficients
% 
% std1Flag = 0;
% std2Flag = 0;
% if abs(m) < 1       % subsonic edge
%     if abs(R1) < delSmall || abs(R2) < delSmall || abs(z) < delSmall
%     % Field point is close to panel edgepoint Mach cones
%     % Or field point is on panel plane
%         if abs(R1) < delSmall
%             w01 = 0;
%         else
%             w01 = (m/sqrt(1-m^2)) * atan2(m*s1 - sm1, R1*sqrt(1-m^2));
%         end
%         
%         if abs(R2) < delSmall
%             w02 = 0;
%         else
%             w02 = (m/sqrt(1-m^2)) * atan2(m*s2 - sm2, R2*sqrt(1-m^2));
%         end
%         
%         w0 = w01 - w02;
%     else
%         Q1 = atan2(z*xmc * (ym1c*R2-ym2c*R1) , z^2*ym1c*ym2c+xmc^2*R1*R2);
%         w0 = (m/sqrt(1-m^2)) * log((ym2c+R2*sqrt(1-m^2))/(ym1c+R1*sqrt(1-m^2)));
%     end
% elseif abs(m) > 1   % supersonic edge
%     if abs(R1) < delSmall || abs(R2) < delSmall
%     % Field point is close to panel edgepoint Mach cones
%         std2Flag = 1;
%         if abs(R1) < delSmall
% %             w01 = sign(s1-lam*sm1) * (pi/2)*sqrt(1-lam^2);
%             w0 = sign(s1-lam*sm1) * (pi/2)*sqrt(1-lam^2);
% %         else
% % %             w01 = (m/sqrt(1-m^2)) * atan2(m*s1 - sm1, R1*sqrt(1-m^2));
% %             w01 = (1/sqrt(1-lam^2)) * atan2(-ym1, R1*sqrt(1-lam^2));
%         end
%         
%         if abs(R2) < delSmall
% %             w02 = sign(s2-lam*sm2) * (pi/2)*sqrt(1-lam^2);
%             w0 = sign(s2-lam*sm2) * (pi/2)*sqrt(1-lam^2);
% %         else
% % %             w02 = (m/sqrt(1-m^2)) * atan2(m*s2 - sm2, R2*sqrt(1-m^2));
% %             w02 = (1/sqrt(1-lam^2)) * atan2(-ym2, R2*sqrt(1-lam^2));
%         end
% %         w0 = w01 - w02;
%     end
%         
%     if abs(z) < delSmall
%         std2Flag = 1;
%         if ym2 > 0 && ym1 < 0
%             Q1 = -sign(z);
%         elseif sign(ym1) == sign(ym2)
%             Q1 = 0;
%         end
%         
%     elseif abs(R1) < delSmall || abs(R2) < delSmall
%         std2Flag = 1;
%         if abs(R1) < delSmall
% %             Q11 = sign(z*(s1-lam*sm1)) * pi/2;
%             Q1 = sign(z*(s1-lam*sm1)) * pi/2;
% %         else
% %             Q11 = atan2(-z*ym1, xm*R1);
%         end
%         
%         if abs(R2) < delSmall
% %             Q12 = sign(z*(s2-lam*sm2)) * pi/2;
%             Q1 = sign(z*(s2-lam*sm2)) * pi/2;
% %         else
% %             Q12 = atan2(-z*ym2, xm*R2);
%         end
% %         Q1 = Q11 - Q12;
%     end
%         
%     if std2Flag == 0
%         Q1 = atan2((z*xm*(ym1*R2-ym2*R1)) , (xm^2*R1*R2+z^2*ym2*ym1));
%         w0 = (1/sqrt(1-lam^2)) * atan2((sqrt(1-lam^2)*(ym1*R2-ym2*R1)) , ...
%             (ym1*ym2+(1-lam^2)*R1*R2));
%     end
% end

%%

% function triGeom = supInfCoeffs(pnt1,pnt2,P,delMag)

% 
% % Define Output Structure
% triGeom.xm   = xm;   triGeom.xmc  = xmc;
% triGeom.ym1  = ym1;  triGeom.ym2  = ym2;
% triGeom.ym1c = ym1c; triGeom.ym2c = ym2c;
% triGeom.R1   = R1;   triGeom.R2   = R2;
% triGeom.m    = m;    triGeom.lam  = lam;

