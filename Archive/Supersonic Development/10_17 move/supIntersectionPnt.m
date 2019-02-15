
% Jake Davis
% 09/15/2018

function pnts = supIntersectionPnt(P,v1,v2,B)

pnts = zeros(1,3);

% Define stuff
i = 1;
x = P(1); y = P(2); z = P(3);
xi1 = v1(1); eta1 = v1(2);
xi2 = v2(1); eta2 = v2(2);
m = (eta2 - eta1) / (xi2 - xi1);

% Intersection 1
x1 = (B*(- B^2*m^2*z^2 + eta1^2 + 2*eta1*m*x - 2*eta1*m*xi1 - 2*eta1*y + ...
    m^2*x^2 - 2*m^2*x*xi1 + m^2*xi1^2 - 2*m*x*y + 2*m*xi1*y + y^2 + z^2)^(1/2) - ...
    x + B^2*m^2*xi1 - B^2*eta1*m + B^2*m*y)/(B^2*m^2 - 1);
y1 = m*(x1-xi1) + eta1;
pnt1 = [x1 y1 0];

% Intersection 2
x2 = -(x + B*(- B^2*m^2*z^2 + eta1^2 + 2*eta1*m*x - 2*eta1*m*xi1 - 2*eta1*y + ...
    m^2*x^2 - 2*m^2*x*xi1 + m^2*xi1^2 - 2*m*x*y + 2*m*xi1*y + y^2 + z^2)^(1/2) - ...
    B^2*m^2*xi1 + B^2*eta1*m - B^2*m*y)/(B^2*m^2 - 1);
y2 = m*(x2-xi1) + eta1;
pnt2 = [x2 y2 0];

pntMin = [x - sqrt(B^2*z^2), y, 0];
c0 = [1 0 0];

% Check if point 1 is on panel edge
if abs(eta1 - y1) < abs(eta2 - eta1) && dot(pntMin - pnt1, c0) >= 0
    pnts(i,:) = pnt1;
    i = i + 1;
end

% Check if point 2 is on panel edge
if abs(eta1 - y2) < abs(eta2 - eta1) && dot(pntMin - pnt2, c0) >= 0
    pnts(i,:) = pnt2;
end

end

%%

% % Compute intersection points
% xPnt1 = (m*y - eta1*m - x + (eta1^2 + 2*eta1*m*x - 2*eta1*m*xi1 - ...
%         2*eta1*y + m^2*x^2 - 2*m^2*x*xi1 + m^2*xi1^2 - m^2*z^2 - ...
%         2*m*x*y + 2*m*xi1*y + y^2 + z^2)^(1/2) + m^2*xi1)/(m^2 - 1);
% yPnt1 = m*(xPnt1-xi1) + eta1;
% 
% xPnt2 = -(x + eta1*m - m*y + (eta1^2 + 2*eta1*m*x - 2*eta1*m*xi1 - ...
%         2*eta1*y + m^2*x^2 - 2*m^2*x*xi1 + m^2*xi1^2 - m^2*z^2 - ...
%         2*m*x*y + 2*m*xi1*y + y^2 + z^2)^(1/2) - m^2*xi1)/(m^2 - 1);
% yPnt2 = m*(xPnt2-xi1) + eta1;
% 
% pnt1 = [xPnt1 yPnt1 0];
% pnt2 = [xPnt2 yPnt2 0];


%%

% if side == 1
%     eta = m*(xi-xi1) + eta1;
% else
%     eta = -(m*(xi-xi1) + eta1);
% end
    
% if side == 1
%     xi = (m*y - eta1*m - x + (eta1^2 + 2*eta1*m*x - 2*eta1*m*xi1 - ...
%         2*eta1*y + m^2*x^2 - 2*m^2*x*xi1 + m^2*xi1^2 - m^2*z^2 - ...
%         2*m*x*y + 2*m*xi1*y + y^2 + z^2)^(1/2) + m^2*xi1)/(m^2 - 1);
% else
%     xi = -(x + eta1*m - m*y + (eta1^2 + 2*eta1*m*x - 2*eta1*m*xi1 - ...
%         2*eta1*y + m^2*x^2 - 2*m^2*x*xi1 + m^2*xi1^2 - m^2*z^2 - ...
%         2*m*x*y + 2*m*xi1*y + y^2 + z^2)^(1/2) - m^2*xi1)/(m^2 - 1);
% end

%%

% syms m xi xi1 eta1 y x z
% 
% eqn1 = m*(xi-xi1) + eta1 == y + sqrt((xi-x)^2 - z^2);
% eqn2 = m*(xi-xi1) + eta1 == y + sqrt((xi-x)^2 - z^2);
% 
% solve(eqn1, xi)
% solve(eqn2, xi)