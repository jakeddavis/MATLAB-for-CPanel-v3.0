
% Jake Davis
% 11/26/2018

function [Q1,w0] = supEdgeInfSupInters(R1,R2,ym1,ym2,xm,lam,z,delSmall)

ym1neg = -ym1;
ym2neg = -ym2;
% ym1neg = ym1;
% ym2neg = ym2;
% m = 1/lam;

if R1 > delSmall && R2 > delSmall
    Q1 = atan2(z*xm * (ym1neg*R2-ym2neg*R1) , z^2*ym2neg*ym1neg + xm^2*R1*R2);
    w0 = (1/sqrt(1-lam^2)) * atan2((sqrt(1-lam^2)*(ym1neg*R2-ym2neg*R1)) , ...
        (ym1neg*ym2neg+(1-lam^2)*R1*R2));
else
    if R1 > delSmall
        Q11 = atan2(z*ym1, xm*R1);
        w01 = (1/sqrt(1-lam^2)) * atan2(ym1, (R1*sqrt(1-lam^2)));
    else
        Q11 = sign(z*ym1) * pi/2;
        w01 = sign(ym1) * pi/(2*sqrt(1-lam^2));
%         Q11 = 0;
%         w01 = 0;
    end
    if R2 > delSmall
        Q12 = atan2(z*ym2, xm*R2);
        w02 = (1/sqrt(1-lam^2)) * atan2(ym2, (R2*sqrt(1-lam^2)));
    else
        Q12 = sign(z*ym2) * pi/2;
        w02 = sign(ym2) * pi/(2*sqrt(1-lam^2));
%         Q12 = 0;
%         w02 = 0;
    end
%     Q1 = Q11 - Q12;
    Q1 = Q12 - Q11;
%     w0 = w01 - w02;
    w0 = w02 - w01;
end

end

%% Checks

%     Q1t = atan2(z*xm * (ym1neg*R2-ym2neg*R1) , z^2*ym2neg*ym1neg + xm^2*R1*R2);
%     Q1t2 = atan2(z*ym1, xm*R1) - atan2(z*ym2, xm*R2);
%     Q1t3 = atan2(z*ym2, xm*R2) - atan2(z*ym1, xm*R1);

%     Q11t = 0.5 * atan2((2*z*R1*(xm*s1-z^2/m)), (2*z^2*R1^2-(xm^2+z^2/m-z^2)*(s1^2+z^2)))
%     Q12t = 0.5 * atan2((2*z*R2*(xm*s2-z^2/m)), (2*z^2*R2^2-(xm^2+z^2/m-z^2)*(s2^2+z^2)))

%% - 10/16/18

% % m = 1/lam;
% % check = xm^2 + z^2/m^2 - z^2;
% 
% if R1 > delSmall && R2 > delSmall
%     Q1 = atan2(z*xm * (ym1neg*R2-ym2neg*R1) , z^2*ym2neg*ym1neg + xm^2*R1*R2);
%     w0 = (1/sqrt(1-lam^2)) * atan2((sqrt(1-lam^2)*(ym1neg*R2-ym2neg*R1)) , ...
%         (ym1neg*ym2neg+(1-lam^2)*R1*R2));
%     
% %     Q11 = atan2(z*ym1, xm*R1);
% %     Q12 = atan2(z*ym2, xm*R2);
% %     Q1t = Q12 - Q11;
% else
% %     Q11 = 0;
% %     Q12 = 0;
% %     w01 = 0;
% %     w02 = 0;
%     if R1 > delSmall
% %         Q1 = atan2(z*ym1, xm*R1);
%         Q11 = atan2(z*ym1, xm*R1);
% %         w0 = (1/sqrt(1-lam^2)) * atan2(ym1, (R1*sqrt(1-lam^2)));
%         w01 = (1/sqrt(1-lam^2)) * atan2(ym1, (R1*sqrt(1-lam^2)));
%     else
%         Q11 = sign(z*ym1) * pi/2;
%         w01 = sign(ym1) * pi/(2*sqrt(1-lam^2));
%     end
%     if R2 > delSmall
% %         Q1 = atan2(z*ym2, xm*R2);
%         Q12 = atan2(z*ym2, xm*R2);
% %         w0 = (1/sqrt(1-lam^2)) * atan2(ym2, (R2*sqrt(1-lam^2)));
%         w02 = (1/sqrt(1-lam^2)) * atan2(ym2, (R2*sqrt(1-lam^2)));
%     else
%         Q12 = sign(z*ym2) * pi/2;
%         w02 = sign(ym2) * pi/(2*sqrt(1-lam^2));
%     end
% %     Q1 = Q11 - Q12;
%     Q1 = Q12 - Q11;
% %     Q1 = atan2(z*xm * (ym1neg*R2-ym2neg*R1) , z^2*ym2neg*ym1neg + xm^2*R1*R2);
% %     w0 = w01 - w02;
%     w0 = w02 - w01;
% %     w0 = (1/sqrt(1-lam^2)) * atan2((sqrt(1-lam^2)*(ym1neg*R2-ym2neg*R1)) , ...
% %         (ym1neg*ym2neg+(1-lam^2)*R1*R2));
% end

%%

% sign(z*ym1)
% sign(z*ym2)
% Q11t = atan2(z*ym1, xm*R1);
% Q12t = atan2(z*ym2, xm*R2);
% Q1t = Q12t - Q11t;
% 
% w01t = (1/sqrt(1-lam^2)) * atan2(ym1, (R1*sqrt(1-lam^2)));
% w02t = (1/sqrt(1-lam^2)) * atan2(ym2, (R2*sqrt(1-lam^2)));
% w0t = w02t - w01t;



% w0t = w01t - w02t;

% if R1 > delSmall && R2 > delSmall
% elseif R1 < delSmall
%     Q1 = sign(z*(s1-lam*sm1)) * pi/2;
%     w0 = sign(s1-lam*sm1) * (pi/2)*sqrt(1-lam^2);
% elseif R2 < delSmall
%     Q1 = sign(z*(s2-lam*sm2)) * pi/2;
%     w0 = sign(s2-lam*sm2) * (pi/2)*sqrt(1-lam^2);
% end
