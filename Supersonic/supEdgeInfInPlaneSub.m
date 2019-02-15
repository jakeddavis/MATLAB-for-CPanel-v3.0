
% Jake Davis
% 10/10/2018

function [Q1,w0] = supEdgeInfInPlaneSub(R1,R2,xmc,ym1c,ym2c,m,z,delSmall,mFlag)

ym1c_neg = -ym1c;
ym2c_neg = -ym2c;

if R1 > delSmall && R2 > delSmall
    Q1 = 0;
    if mFlag == 1
        w0 = (1/sqrt(1-m^2)) * log((ym2c_neg+R2*sqrt(1-m^2))/(ym1c_neg+R1*sqrt(1-m^2)));
        %         w0t = log((ym2c_neg+R2)/(ym1c_neg+R1));
    else
        w0 = (m/sqrt(1-m^2)) * log((ym2c_neg+R2*sqrt(1-m^2))/(ym1c_neg+R1*sqrt(1-m^2)));
    end
else
    Q11 = 0;
    Q12 = 0;
    w01 = 0;
    w02 = 0;
    
    if R1 > delSmall
        Q11 = sign(z*xmc) * pi/2;
        w01 = (m/(2*sqrt(1-m^2))) * log((-ym1c+R1*sqrt(1-m^2)) / (-ym1c-R1*sqrt(1-m^2)));
    end
    if R2 > delSmall
        Q12 = sign(z*xmc) * pi/2;
        w02 = (m/(2*sqrt(1-m^2))) * log((-ym2c+R2*sqrt(1-m^2)) / (-ym2c-R2*sqrt(1-m^2)));
    end
    
    Q1 = Q12 - Q11;
    w0 = w02 - w01;
end

end

%% 01/27/2019

% Q11 = 0;
% Q12 = 0;
% w01 = 0;
% w02 = 0;
% %%%%%%%%%%%%%%%%% not sure about this split condition
% if R1 > delSmall
%     Q11 = sign(z*xmc) * pi/2;
%     w01 = (m/(2*sqrt(1-m^2))) * log((-ym1c+R1*sqrt(1-m^2)) / (-ym1c-R1*sqrt(1-m^2)));
% end
% if R2 > delSmall
%     Q12 = sign(z*xmc) * pi/2;
%     w02 = (m/(2*sqrt(1-m^2))) * log((-ym2c+R2*sqrt(1-m^2)) / (-ym2c-R2*sqrt(1-m^2)));
% end
% Q1 = Q12 - Q11;
% Q1t1 = -pi * sign(z);
% Q1t2 = Q11 - Q12;
% w0 = w02 - w01;
% w0t = w01 - w02;
