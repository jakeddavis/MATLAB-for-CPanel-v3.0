
% Jake Davis
% 10/18/2018

function edgeFlag = supEdgeCheck(v1,v2,P,B)

% Define stuff
x = P(1); y = P(2); z = P(3);
xi1 = v1(1); eta1 = v1(2); zeta1 = v1(3);
xi2 = v2(1); eta2 = v2(2); zeta2 = v2(3);

m = (eta2 - eta1) / (xi2 - xi1);
n = (zeta2 - zeta1) / (eta2 - eta1);
l = (zeta2 - zeta1) / (xi2 - xi1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% look into imaginaries
        %%%% --> they happen when there is no intersection
        %%%% --> probably best to project Mach cone onto plane of panel so
        %%%%     that z = 0 (also maybe rotate to align edge with x/y-axis)

%%%%% need to also account for when edge is constant in 2 directions

if xi1 == xi2
    x1 = xi1;
    x2 = xi1;
    if eta1 == eta2
        y1 = eta1;
        z1 = z + sqrt((xi1-x)^2/B^2 - (eta1-y)^2);
        y2 = eta1;
        z2 = z - sqrt((xi1-x)^2/B^2 - (eta1-y)^2);
    elseif zeta1 == zeta2
        z1 = zeta1;
        y1 = y + sqrt((xi1-x)^2/B^2 - (zeta1-z)^2);
        z2 = zeta1;
        y2 = y - sqrt((xi1-x)^2/B^2 - (zeta1-z)^2);
    else
        % Intersection 1
        z1 = (B*zeta1 + n*(- B^2*eta1^2*n^2 + 2*B^2*eta1*n^2*y - 2*B^2*eta1*n*z + 2*B^2*eta1*n*zeta1 - B^2*n^2*y^2 + 2*B^2*n*y*z - 2*B^2*n*y*zeta1 - B^2*z^2 + 2*B^2*z*zeta1 - B^2*zeta1^2 + n^2*x^2 - 2*n^2*x*xi1 + n^2*xi1^2 + x^2 - 2*x*xi1 + xi1^2)^(1/2) - B*eta1*n + B*n*y + B*n^2*z)/(B*n^2 + B);
        y1 = (z1-zeta1)/n + eta1;
        
        % Intersection 2
        z2 = (B*zeta1 - n*(- B^2*eta1^2*n^2 + 2*B^2*eta1*n^2*y - 2*B^2*eta1*n*z + 2*B^2*eta1*n*zeta1 - B^2*n^2*y^2 + 2*B^2*n*y*z - 2*B^2*n*y*zeta1 - B^2*z^2 + 2*B^2*z*zeta1 - B^2*zeta1^2 + n^2*x^2 - 2*n^2*x*xi1 + n^2*xi1^2 + x^2 - 2*x*xi1 + xi1^2)^(1/2) - B*eta1*n + B*n*y + B*n^2*z)/(B*n^2 + B);
        y2 = (z2-zeta1)/n + eta1;
    end
    
    pnt1 = [x1 y1 z1];
    plot3(x1,y1,z1,'b*')
    
    pnt2 = [x2 y2 z2];
    plot3(x2,y2,z2,'r*')
elseif eta1 == eta2 % I think this isn't really needed, in fs direction
    % won't ever intersect without leading point being inside
    % Intersection 1
    y1 = eta1;
    z1 = (l*xi1 - zeta1 + l*(- B^2*eta1^2*l^2 + B^2*eta1^2 + 2*B^2*eta1*l^2*y - 2*B^2*eta1*y - B^2*l^2*y^2 + B^2*y^2 + l^2*xi1^2 + 2*l*xi1*z - 2*l*xi1*zeta1 + z^2 - 2*z*zeta1 + zeta1^2)^(1/2) + l^2*z)/(l^2 - 1);
    x1 = (z1-zeta1)/l + xi1;
    pnt1 = [x1 y1 z1];
    plot3(x1,y1,z1,'b*')
    
    % Intersection 2
    y2 = eta1;
    z2 = -(zeta1 - l*xi1 + l*(- B^2*eta1^2*l^2 + B^2*eta1^2 + 2*B^2*eta1*l^2*y - 2*B^2*eta1*y - B^2*l^2*y^2 + B^2*y^2 + l^2*xi1^2 + 2*l*xi1*z - 2*l*xi1*zeta1 + z^2 - 2*z*zeta1 + zeta1^2)^(1/2) - l^2*z)/(l^2 - 1);
    x2 = (z2-zeta1)/l + xi1;
    pnt2 = [x2 y2 z2];
    plot3(x2,y2,z2,'r*')
elseif zeta1 == zeta2
    % Intersection 1
    z1 = zeta1;
    x1 = (B*(- B^2*m^2*z^2 + 2*B^2*m^2*z*zeta1 - B^2*m^2*zeta1^2 + eta1^2 + 2*eta1*m*x - 2*eta1*m*xi1 - 2*eta1*y + m^2*x^2 - 2*m^2*x*xi1 + m^2*xi1^2 - 2*m*x*y + 2*m*xi1*y + y^2 + z^2 - 2*z*zeta1 + zeta1^2)^(1/2) - x + B^2*m^2*xi1 - B^2*eta1*m + B^2*m*y)/(B^2*m^2 - 1);
    y1 = m*(x1-xi1) + eta1;
    pnt1 = [x1 y1 z1];
    plot3(x1,y1,z1,'b*')
    
    % Intersection 2
    z2 = zeta1;
    x2 = -(x + B*(- B^2*m^2*z^2 + 2*B^2*m^2*z*zeta1 - B^2*m^2*zeta1^2 + eta1^2 + 2*eta1*m*x - 2*eta1*m*xi1 - 2*eta1*y + m^2*x^2 - 2*m^2*x*xi1 + m^2*xi1^2 - 2*m*x*y + 2*m*xi1*y + y^2 + z^2 - 2*z*zeta1 + zeta1^2)^(1/2) - B^2*m^2*xi1 + B^2*eta1*m - B^2*m*y)/(B^2*m^2 - 1);
    y2 = m*(x2-xi1) + eta1;
    pnt2 = [x2 y2 z2];
    plot3(x2,y2,z2,'r*')
else    
    % Intersection 1
    x1 = (y - eta1 + m*xi1 + (y - eta1 - m*x + m*xi1 + B*m*(- B^2*eta1^2*m^2*n^2 + 2*B^2*eta1*m^2*n^2*y - 2*B^2*eta1*m^2*n*z + 2*B^2*eta1*m^2*n*zeta1 - B^2*m^2*n^2*y^2 + 2*B^2*m^2*n*y*z - 2*B^2*m^2*n*y*zeta1 - B^2*m^2*z^2 + 2*B^2*m^2*z*zeta1 - B^2*m^2*zeta1^2 + eta1^2 + 2*eta1*m*x - 2*eta1*m*xi1 - 2*eta1*y + m^2*n^2*x^2 - 2*m^2*n^2*x*xi1 + m^2*n^2*xi1^2 + m^2*x^2 - 2*m^2*x*xi1 + m^2*xi1^2 - 2*m*n*x*z + 2*m*n*x*zeta1 + 2*m*n*xi1*z - 2*m*n*xi1*zeta1 - 2*m*x*y + 2*m*xi1*y + y^2 + z^2 - 2*z*zeta1 + zeta1^2)^(1/2) + B^2*eta1*m^2*n^2 - B^2*m^2*n^2*y + B^2*m^2*n*z - B^2*m^2*n*zeta1)/(B^2*m^2*n^2 + B^2*m^2 - 1))/m;
    z1 = zeta1 - eta1*n + n*y + (n*(y - eta1 - m*x + m*xi1 + B*m*(- B^2*eta1^2*m^2*n^2 + 2*B^2*eta1*m^2*n^2*y - 2*B^2*eta1*m^2*n*z + 2*B^2*eta1*m^2*n*zeta1 - B^2*m^2*n^2*y^2 + 2*B^2*m^2*n*y*z - 2*B^2*m^2*n*y*zeta1 - B^2*m^2*z^2 + 2*B^2*m^2*z*zeta1 - B^2*m^2*zeta1^2 + eta1^2 + 2*eta1*m*x - 2*eta1*m*xi1 - 2*eta1*y + m^2*n^2*x^2 - 2*m^2*n^2*x*xi1 + m^2*n^2*xi1^2 + m^2*x^2 - 2*m^2*x*xi1 + m^2*xi1^2 - 2*m*n*x*z + 2*m*n*x*zeta1 + 2*m*n*xi1*z - 2*m*n*xi1*zeta1 - 2*m*x*y + 2*m*xi1*y + y^2 + z^2 - 2*z*zeta1 + zeta1^2)^(1/2) + B^2*eta1*m^2*n^2 - B^2*m^2*n^2*y + B^2*m^2*n*z - B^2*m^2*n*zeta1))/(B^2*m^2*n^2 + B^2*m^2 - 1);
    y1 = m*(x1-xi1) + eta1;
    pnt1 = [x1 y1 z1];
    plot3(x1,y1,z1,'b*')
    
    % Intersection 2
    x2 = -(eta1 - y - m*xi1 + (eta1 - y + m*x - m*xi1 + B*m*(- B^2*eta1^2*m^2*n^2 + 2*B^2*eta1*m^2*n^2*y - 2*B^2*eta1*m^2*n*z + 2*B^2*eta1*m^2*n*zeta1 - B^2*m^2*n^2*y^2 + 2*B^2*m^2*n*y*z - 2*B^2*m^2*n*y*zeta1 - B^2*m^2*z^2 + 2*B^2*m^2*z*zeta1 - B^2*m^2*zeta1^2 + eta1^2 + 2*eta1*m*x - 2*eta1*m*xi1 - 2*eta1*y + m^2*n^2*x^2 - 2*m^2*n^2*x*xi1 + m^2*n^2*xi1^2 + m^2*x^2 - 2*m^2*x*xi1 + m^2*xi1^2 - 2*m*n*x*z + 2*m*n*x*zeta1 + 2*m*n*xi1*z - 2*m*n*xi1*zeta1 - 2*m*x*y + 2*m*xi1*y + y^2 + z^2 - 2*z*zeta1 + zeta1^2)^(1/2) - B^2*eta1*m^2*n^2 + B^2*m^2*n^2*y - B^2*m^2*n*z + B^2*m^2*n*zeta1)/(B^2*m^2*n^2 + B^2*m^2 - 1))/m;
    z2 = zeta1 - eta1*n + n*y + (n*(y - eta1 - m*x + m*xi1 + B*m*(- B^2*eta1^2*m^2*n^2 + 2*B^2*eta1*m^2*n^2*y - 2*B^2*eta1*m^2*n*z + 2*B^2*eta1*m^2*n*zeta1 - B^2*m^2*n^2*y^2 + 2*B^2*m^2*n*y*z - 2*B^2*m^2*n*y*zeta1 - B^2*m^2*z^2 + 2*B^2*m^2*z*zeta1 - B^2*m^2*zeta1^2 + eta1^2 + 2*eta1*m*x - 2*eta1*m*xi1 - 2*eta1*y + m^2*n^2*x^2 - 2*m^2*n^2*x*xi1 + m^2*n^2*xi1^2 + m^2*x^2 - 2*m^2*x*xi1 + m^2*xi1^2 - 2*m*n*x*z + 2*m*n*x*zeta1 + 2*m*n*xi1*z - 2*m*n*xi1*zeta1 - 2*m*x*y + 2*m*xi1*y + y^2 + z^2 - 2*z*zeta1 + zeta1^2)^(1/2) + B^2*eta1*m^2*n^2 - B^2*m^2*n^2*y + B^2*m^2*n*z - B^2*m^2*n*zeta1))/(B^2*m^2*n^2 + B^2*m^2 - 1);
    y2 = m*(x2-xi1) + eta1;
    pnt2 = [x2 y2 z2];
    plot3(x2,y2,z2,'r*')
end

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
