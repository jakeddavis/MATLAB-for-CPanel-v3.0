
% Jake Davis
% 10/13/2018

function edgeFlag = supEdgeCheck(v1,v2,P,B)

% Define stuff
x = P(1); y = P(2); z = P(3);
xi1 = v1(1); eta1 = v1(2); zeta1 = v1(3);
xi2 = v2(1); eta2 = v2(2); zeta2 = v2(3);
m = (eta2 - eta1) / (xi2 - xi1);
n = (zeta2 - zeta1) / (eta2 - eta1);

% Intersection 1
x1 = (B*(- B^2*m^2*z^2 + eta1^2 + 2*eta1*m*x - 2*eta1*m*xi1 - 2*eta1*y + ...
    m^2*x^2 - 2*m^2*x*xi1 + m^2*xi1^2 - 2*m*x*y + 2*m*xi1*y + y^2 + z^2)^(1/2) - ...
    x + B^2*m^2*xi1 - B^2*eta1*m + B^2*m*y)/(B^2*m^2 - 1);
y1 = m*(x1-xi1) + eta1;
z1 = n*(y1-eta1) + zeta1;
% pnt1 = [x1 y1 0];
pnt1 = [x1 y1 z1];
plot3(x1,y1,0,'g*')

% Intersection 2
x2 = -(x + B*(- B^2*m^2*z^2 + eta1^2 + 2*eta1*m*x - 2*eta1*m*xi1 - 2*eta1*y + ...
    m^2*x^2 - 2*m^2*x*xi1 + m^2*xi1^2 - 2*m*x*y + 2*m*xi1*y + y^2 + z^2)^(1/2) - ...
    B^2*m^2*xi1 + B^2*eta1*m - B^2*m*y)/(B^2*m^2 - 1);
y2 = m*(x2-xi1) + eta1;
pnt2 = [x2 y2 0];
plot3(x2,y2,0,'g*')

pntMin = [x - sqrt(B^2*z^2), y, 0];
c0 = [1 0 0];

% Check if point(s) are on panel edge
edgeFlag = 0;
if abs(eta2 - y1) < abs(eta2 - eta1) && dot(pntMin - pnt1, c0) >= 0
    edgeFlag = 1;
elseif abs(eta2 - y2) < abs(eta2 - eta1) && dot(pntMin - pnt2, c0) >= 0
    edgeFlag = 1;
end

end

%% old post int. pnt calc check scheme

% [pnt1,pnt2] = supIntersectionPnt(P,vert1,vert2);
% c0 = [1 0 0];
% 
% edgeFlag = 0;
% if dot(P - pnt1, c0) >= 0
%     if vert1(2) < vert2(2)
%         if pnt1(2) > vert1(2) && pnt1(2) < vert2(2)
%             edgeFlag = 1;
%         end
%     else
%         if pnt1(2) < vert1(2) && pnt1(2) > vert2(2)
%             edgeFlag = 1;
%         end
%     end
% elseif dot(P - pnt2, c0) >= 0
%     if vert1(2) < vert2(2)
%         if pnt2(2) > vert1(2) && pnt2(2) < vert2(2)
%             edgeFlag = 1;
%         end
%     else
%         if pnt2(2) < vert1(2) && pnt2(2) > vert2(2)
%             edgeFlag = 1;
%         end
%     end
% end
