
% Jake Davis
% 10/10/2018

function [Q1,w0] = supEdgeInfSub1(R1,R2,ym1c,ym2c,xmc,m,z,delSmall)

ym1c_neg = -ym1c;
ym2c_neg = -ym2c;

if R1 > delSmall && R2 > delSmall
    Q1 = atan2(z*xmc * (ym1c_neg*R2-ym2c_neg*R1) , z^2*ym1c_neg*ym2c_neg+xmc^2*R1*R2);
    w0 = (m/sqrt(1-m^2)) * log((ym2c_neg+R2*sqrt(1-m^2))/(ym1c_neg+R1*sqrt(1-m^2)));
else
    if R1 > delSmall
        Q11 = sign(z)*atan2(xmc*R1, -abs(z)*ym1c);
        w01 = (m/(2*sqrt(1-m^2))) * log((-ym1c+R1*sqrt(1-m^2)) / (-ym1c-R1*sqrt(1-m^2)));
    else
        Q11 = 0;
        w01 = 0;
    end
    if R2 > delSmall
        Q12 = sign(z)*atan2(xmc*R2, -abs(z)*ym2c);
        w02 = (m/(2*sqrt(1-m^2))) * log((-ym2c+R2*sqrt(1-m^2)) / (-ym2c-R2*sqrt(1-m^2)));
    else
        Q12 = 0;
        w02 = 0;
    end
    Q1 = Q12 - Q11;
%     Q1 = Q11 - Q12;
    w0 = w02 - w01;
%     w0 = w01 - w02;
end

end

%% 10/16/2018

% if R1 > delSmall && R2 > delSmall
%     Q1 = atan2(z*xmc * (ym1c_neg*R2-ym2c_neg*R1) , z^2*ym1c_neg*ym2c_neg+xmc^2*R1*R2);
%     w0 = (m/sqrt(1-m^2)) * log((ym2c_neg+R2*sqrt(1-m^2))/(ym1c_neg+R1*sqrt(1-m^2)));
% 
% %     Q11 = sign(z)*atan2(xmc*R1, -abs(z)*ym1c);
% %     Q12 = sign(z)*atan2(xmc*R2, -abs(z)*ym2c);
% %     Q1t = Q12 - Q11;
% else
% %     Q11 = 0;
% %     Q12 = 0;
% %     w01 = 0;
% %     w02 = 0;
%     if R1 > delSmall
% %         Q1t = sign(z)*atan((xmc*R1) / (-abs(z)*ym1c));
% %         Q1 = sign(z)*atan2(xmc*R1, -abs(z)*ym1c);
%         Q11 = sign(z)*atan2(xmc*R1, -abs(z)*ym1c);
% %         w0 = (m/(2*sqrt(1-m^2))) * log((-ym1c+R1*sqrt(1-m^2)) / (-ym1c-R1*sqrt(1-m^2)));
%         w01 = (m/(2*sqrt(1-m^2))) * log((-ym1c+R1*sqrt(1-m^2)) / (-ym1c-R1*sqrt(1-m^2)));
%     else
%         Q11 = 0;
% %         Q11 = sign(xmc) * pi/2;
%         w01 = 0;
%     end
%     if R2 > delSmall
% %         Q1 = sign(z)*atan2(xmc*R2, -abs(z)*ym2c);
%         Q12 = sign(z)*atan2(xmc*R2, -abs(z)*ym2c);
% %         w0 = (m/(2*sqrt(1-m^2))) * log((-ym2c+R2*sqrt(1-m^2)) / (-ym2c-R2*sqrt(1-m^2)));
%         w02 = (m/(2*sqrt(1-m^2))) * log((-ym2c+R2*sqrt(1-m^2)) / (-ym2c-R2*sqrt(1-m^2)));
%     else
%         Q12 = 0;
% %         Q12 = sign(xmc) * pi/2;
%         w02 = 0;
%     end
%     Q1 = Q12 - Q11;
% %     Q1 = Q11 - Q12;
% %     Q1t2 = atan2(z*xmc * (ym1c_neg*R2-ym2c_neg*R1) , z^2*ym1c_neg*ym2c_neg+xmc^2*R1*R2);
%     w0 = w02 - w01;
% %     w0 = w01 - w02;
% %     w0t2 = (m/sqrt(1-m^2)) * log((ym2c_neg+R2*sqrt(1-m^2))/(ym1c_neg+R1*sqrt(1-m^2)));
% end
%     
%     
% % Q11t = sign(z)*atan2(xmc*R1, abs(z)*ym1c);
% % Q12t = sign(z)*atan2(xmc*R2, abs(z)*ym2c);
% % Q1t = Q12t - Q11t;

