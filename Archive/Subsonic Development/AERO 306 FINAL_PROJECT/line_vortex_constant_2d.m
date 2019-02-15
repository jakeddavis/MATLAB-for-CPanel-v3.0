
% 306

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the induced velocity at a point caused by a line vortex %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% u ---------- x-component of the induced velocity
% v ---------- y-component of the induced velocity
% gamma ------ vortex strength
% xpanel ----- vector of x-coordinates of the start and end points of the panel
% ypanel ----- vector of y-coordinates of the start and end points of the panel
% xpt -------- x-coordinate of the point we want the induced velocity at
% ypt -------- y-coordinate of the point we want the induced velocity at

function [u, v] = line_vortex_constant_2d(gamma, xpanel, ypanel, xpt, ypt)

  % calculate theta
  theta = atan2(ypanel(2) - ypanel(1), xpanel(2) - xpanel(1));
  
  % finding '1' info
  xi1 = (xpt - xpanel(1));
  eta1 = (ypt - ypanel(1));
  r1 = sqrt(xi1^2 + eta1^2);
  beta1 = atan2(eta1, xi1) - theta;
  
  % finding '2' info
  xi2 = (xpt - xpanel(2));
  eta2 = (ypt - ypanel(2));
  r2 = sqrt(xi2^2 + eta2^2);
  beta2 = atan2(eta2, xi2) - theta;
  
  % calculate xi- and eta-velocity terms
  I0 = log(r1/r2);
  I1 = beta2 - beta1;
  if (I1 < -pi)
    I1 = I1 + 2*pi;
  elseif (I1 > pi)
    I1 = I1 - 2*pi;
  end
  
  uxi = (gamma/(2*pi))*I1;
  ueta = -(gamma/(2*pi))*I0;
  
  u = uxi*cos(theta) - ueta*sin(theta);
  v = uxi*sin(theta) + ueta*cos(theta);
end

