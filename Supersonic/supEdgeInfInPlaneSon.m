
% Jake Davis
% 10/10/2018

function [Q1,w0] = supEdgeInfInPlaneSon(R1,R2,ym1,ym2,lam,z,delSmall)

zr = (ym1*R2 - ym2*R1) / (ym1*ym2 + (1-lam^2)*R1*R2);
w0 = zr * (1 - ((1-lam^2)*zr^2)/3 + ((1-lam^2)^2*zr^4)/5 ...
    - ((1-lam^2)^3*zr^6)/7); % +... (atan2 expansion)
if abs(m) < 1
    Q11 = 0;
    Q12 = 0;
    %%%%%%%%%%%%%%%%% not sure about this split condition
    if R1 > delSmall
        Q11 = sign(z*xmc) * pi/2;
    end
    if R2 > delSmall
        Q12 = sign(z*xmc) * pi/2;
    end
    Q1 = Q12 - Q11;
    Q1t = -pi * sign(z);
    %                 Q1 = Q11 - Q12;
elseif abs(m) > 1
    Q11 = sign(z*ym1) * pi/2;
    Q12 = atan2(z*ym2, xm*R2);
    Q1 = Q12 - Q11;
    
    if ym1 < 0 && ym2 > 0
        Q1 = -pi * sign(z);
        % else Q1 = 0;
    end
end

end