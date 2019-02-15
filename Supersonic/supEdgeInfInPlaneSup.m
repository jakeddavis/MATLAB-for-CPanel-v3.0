
% Jake Davis
% 10/10/2018

function [Q1,w0] = supEdgeInfInPlaneSup(R1,R2,xm,ym1,ym2,lam,z,delSmall)

ym1neg = -ym1;
ym2neg = -ym2;

Q1 = 0;

if R1 > delSmall && R2 > delSmall
    w0 = (1/sqrt(1-lam^2)) * atan2((sqrt(1-lam^2)*(ym1neg*R2-ym2neg*R1)) , ...
        (ym1neg*ym2neg+(1-lam^2)*R1*R2));
else
    if R1 > delSmall
%         Q11 = atan2(z*ym1, xm*R1);
        w01 = (1/sqrt(1-lam^2)) * atan2(ym1, (R1*sqrt(1-lam^2)));
    else
%         Q11 = sign(z*ym1) * pi/2;
        w01 = sign(ym1) * pi/(2*sqrt(1-lam^2));
    end
    if R2 > delSmall
%         Q12 = atan2(z*ym2, xm*R2);
        w02 = (1/sqrt(1-lam^2)) * atan2(ym2, (R2*sqrt(1-lam^2)));
    else
%         Q12 = sign(z*ym2) * pi/2;
        w02 = sign(ym2) * pi/(2*sqrt(1-lam^2));
    end
%     Q1 = Q12 - Q11;
    w0 = w02 - w01;
end




% if R1 > delSmall
%     w01 = (1/sqrt(1-lam^2)) * atan2(ym1, (R1*sqrt(1-lam^2)));
% else
%     w01 = sign(ym1) * pi/(2*sqrt(1-lam^2));
% end
% if R2 > delSmall
%     w02 = (1/sqrt(1-lam^2)) * atan2(ym2, (R2*sqrt(1-lam^2)));
% else
%     w02 = sign(ym2) * pi/(2*sqrt(1-lam^2));
% end
% w0 = w02 - w01;
% w0t = w01 - w02;
% 
% Q11 = sign(z*ym1) * pi/2;
% Q12 = sign(z*ym2) * pi/2;
% Q1t1 = Q12 - Q11;
% Q1t2 = Q11 - Q12;

end

%% 01/27/2019

% Q1 = 0;
% if R1 > delSmall
%     w01 = (1/sqrt(1-lam^2)) * atan2(ym1, (R1*sqrt(1-lam^2)));
% else
%     w01 = sign(ym1) * pi/(2*sqrt(1-lam^2));
% end
% if R2 > delSmall
%     w02 = (1/sqrt(1-lam^2)) * atan2(ym2, (R2*sqrt(1-lam^2)));
% else
%     w02 = sign(ym2) * pi/(2*sqrt(1-lam^2));
% end
% w0 = w02 - w01;
% w0t = w01 - w02;
% 
% Q11 = sign(z*ym1) * pi/2;
% Q12 = sign(z*ym2) * pi/2;
% Q1t1 = Q12 - Q11;
% Q1t2 = Q11 - Q12;
