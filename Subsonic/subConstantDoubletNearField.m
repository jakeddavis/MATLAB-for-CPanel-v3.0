
% Jake Davis
% 02/03/2018

% Near field calculations for source influence

function Phi = subConstantDoubletNearField(mu,verts,P)

z = P(3);
Phi = zeros(1,3);
verts(4,:) = verts(1,:);

for i = 1:3
    pnt1 = [verts(i,1) verts(i,2)];
    pnt2 = [verts(i+1,1) verts(i+1,2)];
    
    triGeom = triGeom_func(pnt1,pnt2,P);
    m = triGeom.m;
    e1 = triGeom.e1; e2 = triGeom.e2;
    h1 = triGeom.h1; h2 = triGeom.h2;
    r1 = triGeom.r1; r2 = triGeom.r2;
    
    J = atan((m*e1-h1)/(z*r1)) - atan((m*e2-h2)/(z*r2));
    Phi(i) = J;  
end

% phiSum = sum(Phi)
Phi = (mu/4/pi)*sum(Phi);

end

%%

% J_sum = sum(Phi)/P(3);

% % influenced point
% x = P(1); y = P(2); z = P(3);
% verts(4,:) = verts(1,:);
% 
% Phi = zeros(1,3);

% for i = 1:3
% %     % tri edge pnt 1
% %     x1 = verts(i,1); y1 = verts(i,2);
% %     % tri edge pnt 2
% %     x2 = verts(i+1,1); y2 = verts(i+1,2);
% % 
% %     m = (y2-y1) / (x2-x1);
% %     r1 = sqrt((x-x1)^2 + (y-y1)^2 + z^2);
% %     r2 = sqrt((x-x2)^2 + (y-y2)^2 + z^2);
% %     e1 = (x-x1)^2 + z^2;
% %     e2 = (x-x2)^2 + z^2;
% %     h1 = (x-x1) * (y-y1);
% %     h2 = (x-x2) * (y-y2);
% 
% %     J = atan((m*e1-h1)/(z*r1)) - atan((m*e2-h2)/(z*r2));
%     
%     Phi(i) = J;
% end
% 
% Phi = (mu/4/pi)*sum(Phi);
% 
% end

% [m12,m23,m31] = mFunc(verts);
% [r1,r2,r3]    = rFunc(verts,x,y,z);
% [e1,e2,e3]    = eFunc(verts,x,z);
% [h1,h2,h3]    = hFunc(verts,x,y);

% Phi = (mu/(4*pi)) * ...
%     (atan((m12*e1-h1)/(z*r1)) - atan((m12*e2-h2)/(z*r2)) +...
%     atan((m23*e2-h2)/(z*r2)) - atan((m23*e3-h3)/(z*r3)) +...
%     atan((m31*e3-h3)/(z*r3)) - atan((m31*e1-h1)/(z*r1)));