% 306

% computes the induced velocity at a point by a point vortex

% GAMMA --------- vortex strength
% xv ------------ x-location of vortex
% yv ------------ y-location of vortex
% xpt ----------- x-coordinate of the point we want the induced velocity
% ypt ----------- y-coordinate of the point we want the induced velocity

function [u, v] = point_vortex(GAMMA,xv,yv,xpt,ypt)

    xdif = xpt - xv;
    ydif = ypt - yv;
    
    r = sqrt(xdif^2 + ydif^2);
    
    u = (GAMMA/(2*pi))*(ydif/r^2);
    v = -(GAMMA/(2*pi))*(xdif/r^2);
    
end