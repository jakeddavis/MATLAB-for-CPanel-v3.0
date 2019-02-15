
% Jake Davis
% 10/10/2018

function [Q1,w0] = supEdgeInfVertex(sign,m,lam,R,xm,ym,s,sm,z,delSmall)

if abs(m) < 1
    % subsonic edge
    if R > delSmall
        %%%%%%%%%%%%%% atan or atan2 ?????????
%         Q1 = atan2((-z*ym), (xm*R));
        Q1 = atan((-z*ym) / (xm*R));
        w0 = (m/(2*sqrt(1-m^2))) * log(((sm-m*s+R*sqrt(1-m^2))) / ...
            (sm-m*s-R*sqrt(1-m^2)));
    elseif R < delSmall
        Q1 = 0;
        w0 = 0;
    end
elseif abs(m) > 1
    % supersonic edge
    if R > delSmall
%         Q1 = atan2((-z*ym), (xm*R));
        Q1 = atan((-z*ym) / (xm*R));
        w0 = (1/sqrt(1-lam^2)) * atan2(-ym, (R*sqrt(1-lam^2)));
    elseif R < delSmall
        Q1 = sign(z*(s-lam*sm)) * pi/2;
        w0 = sign(s-lam*sm) * (pi/2)*sqrt(1-lam^2);
    end
end

end
