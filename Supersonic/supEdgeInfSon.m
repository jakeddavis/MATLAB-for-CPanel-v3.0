
% Jake Davis
% 12/13/2018

function [Q1, w0] = supEdgeInfSon(xm,ym1,ym2,R1,R2,lam,z,delSmall)

ym1neg = -ym1;
ym2neg = -ym2;

% zr = (ym1*R2 - ym2*R1) / (ym1*ym2 + (1-lam^2)*R1*R2);
zr = (ym1neg*R2 - ym2neg*R1) / (ym1neg*ym2neg + (1-lam^2)*R1*R2);
w0 = zr * (1 - ((1-lam^2)*zr^2)/3 + ((1-lam^2)^2*zr^4)/5 ...
    - ((1-lam^2)^3*zr^6)/7); % + . . . (atan2 expansion)

if R1 > delSmall && R2 > delSmall    
    Q1 = atan2(z*xm * (ym1neg*R2-ym2neg*R1) , z^2*ym2neg*ym1neg + xm^2*R1*R2);
else
    if R1 > delSmall
        Q11 = atan2(z*ym1, xm*R1);
    else
        Q11 = sign(z*ym1) * pi/2;
    end
    if R2 > delSmall
        Q12 = atan2(z*ym2, xm*R2);
    else
        Q12 = sign(z*ym2) * pi/2;
    end
    
    Q1 = Q12 - Q11;
end

end

%% w0 is only part that is special. Either sub or sup Q1 can be used

% Q11 = 0;
% Q12 = 0;
% % ym1neg = -ym1;
% % ym2neg = -ym2;
% ym1neg = -ym1;
% ym2neg = -ym2;
% ym1c_neg = -ym1c;
% ym2c_neg = -ym2c;
% 
% % zr = (ym1*R2 - ym2*R1) / (ym1*ym2 + (1-lam^2)*R1*R2);
% zr = (ym1neg*R2 - ym2neg*R1) / (ym1neg*ym2neg + (1-lam^2)*R1*R2);
% w0 = zr * (1 - ((1-lam^2)*zr^2)/3 + ((1-lam^2)^2*zr^4)/5 ...
%     - ((1-lam^2)^3*zr^6)/7); % + . . . (atan2 expansion)
% 
% if R1 > delSmall && R2 > delSmall    
%     % sub
%     Q1 = atan2(z*xmc * (ym1c_neg*R2-ym2c_neg*R1) , z^2*ym1c_neg*ym2c_neg+xmc^2*R1*R2);
% %     sup
% %     Q1 = atan2(z*xm * (ym1neg*R2-ym2neg*R1) , z^2*ym2neg*ym1neg + xm^2*R1*R2);
% else
% %     % sup
% %     if R1 > delSmall
% %         Q11 = atan2(z*ym1, xm*R1);
% %     else
% %         Q11 = sign(z*ym1) * pi/2;
% %     end
% %     if R2 > delSmall
% %         Q12 = atan2(z*ym2, xm*R2);
% %     else
% %         Q12 = sign(z*ym2) * pi/2;
% %     end
% %     
% %     Q1 = Q12 - Q11;
%     
%     % sub
%     if R1 > delSmall
%         Q11 = sign(z)*atan2(xmc*R1, -abs(z)*ym1c);
%     end
%     if R2 > delSmall
%         Q12 = sign(z)*atan2(xmc*R2, -abs(z)*ym2c);
%     end
%     
%     Q1 = Q12 - Q11;
% end

