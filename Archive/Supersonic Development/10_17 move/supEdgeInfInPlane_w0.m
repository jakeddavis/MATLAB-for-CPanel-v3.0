
% Jake Davis
% 10/14/2018

function w0 = supEdgeInfInPlane_w0()

if abs(1 - abs(m)) < delSmall % sonic edge
    zr = (ym1*R2 - ym2*R1) / (ym1*ym2 + (1-lam^2)*R1*R2);
    w0 = zr * (1 - ((1-lam^2)*zr^2)/3 + ((1-lam^2)^2*zr^4)/5 ...
        - ((1-lam^2)^3*zr^6)/7); % + . . . (atan2 expansion)
elseif abs(m) < 1 % subsonic edge
    
elseif abs(m) > 1 % supersonic edge
    
end

end