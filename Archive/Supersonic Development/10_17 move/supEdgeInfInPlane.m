
% Jake Davis
% 10/10/2018

function [Q1,w0] = supEdgeInfInPlane(m,R,z,ym1,ym2,delSmall)

w0 = 0; % not true, but mult'pd by z = 0 later
Q1 = 0;
if abs(m) < 1
    % subsonic edge
    if R > delSmall
        Q1 = -pi * sign(z);
    end
elseif abs(m) > 1
    if ym1 < 0 && ym2 > 0
        Q1 = -pi * sign(z);
    end
end

end