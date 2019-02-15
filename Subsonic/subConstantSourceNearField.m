
% Jake Davis
% 02/03/2018

% Near field calculations for source influence

function Phi = subConstantSourceNearField(sigma,verts,P)

z = P(3);
Phi = zeros(1,3);
verts(4,:) = verts(1,:);
for i = 1:3
    pnt1 = [verts(i,1) verts(i,2)];
    pnt2 = [verts(i+1,1) verts(i+1,2)];
    
    triGeom = triGeom_func(pnt1,pnt2,P);
    a  = triGeom.a;
    m  = triGeom.m;  d  = triGeom.d;
    e1 = triGeom.e1; e2 = triGeom.e2;
    h1 = triGeom.h1; h2 = triGeom.h2;
    r1 = triGeom.r1; r2 = triGeom.r2;
    
    Q = log((r1+r2+d)/(r1+r2-d));
    J = atan((m*e1-h1)/(z*r1)) - atan((m*e2-h2)/(z*r2));

    Phi(i) = a*Q - abs(z)*J;   
end

Phi = (-sigma/4/pi)*sum(Phi);

end

%%

%     x1 = pnt1(1); y1 = pnt1(2);
%     x2 = pnt2(1); y2 = pnt2(2);
    
%     d = sqrt((x2-x1)^2 + (y2-y1)^2);
%     m = (y2-y1) / (x2-x1);
%     r1 = sqrt((x-x1)^2 + (y-y1)^2 + z^2);
%     r2 = sqrt((x-x2)^2 + (y-y2)^2 + z^2);
%     e1 = (x-x1)^2 + z^2;
%     e2 = (x-x2)^2 + z^2;
%     h1 = (x-x1) * (y-y1);
%     h2 = (x-x2) * (y-y2);

%     triGeom = triGeom_func(pnt1,pnt2,P);
%     a  = triGeom.a;
%     m  = triGeom.m;  d  = triGeom.d;
%     e1 = triGeom.e1; e2 = triGeom.e2;
%     h1 = triGeom.h1; h2 = triGeom.h2;
%     r1 = triGeom.r1; r2 = triGeom.r2;
%     
% %     if l1 >= 0 && l2 >= 0
% %         Q = log((sqrt(l2^2+g^2)+l2) / (sqrt(l1^2+g^2)+l1));
% %     elseif l1 < 0 && l2 < 0
% %         Q = log((sqrt(l1^2+g^2)-l1) / (sqrt(l2^2+g^2)-l2));
% %     elseif l2 >= 0 && l1 < 0
% %         Q = log(((sqrt(l1^2+g^2)-l1) * (sqrt(l2^2+g^2)+l2))/g^2);
% %     end
%     
% %     J = atan((m*e1-h1)/(z*abs(l1))) - atan((m*e2-h2)/(z*abs(l2)));
%     
% %     x = P(1); y = P(2); z = P(3);
% %     x1 = pnt1(1); y1 = pnt1(2);
% %     x2 = pnt2(1); y2 = pnt2(2);
% %     d = sqrt((x2-x1)^2 + (y2-y1)^2);
% %     r1 = sqrt((x-x1)^2 + (y-y1)^2 + z^2);
% %     r2 = sqrt((x-x2)^2 + (y-y2)^2 + z^2);
% %     R = ((x-x1)*(y2-y1)-(y-y1)*(x2-x1))/d;
%     Q = log((r1+r2+d)/(r1+r2-d));
%     J = atan((m*e1-h1)/(z*r1)) - atan((m*e2-h2)/(z*r2));


%%

% x = P(1); y = P(2); z = P(3);

%     % tri edge pnt 1
%     x1 = verts(i,1); y1 = verts(i,2);
%     % tri edge pnt 2
%     x2 = verts(i+1,1); y2 = verts(i+1,2);

%     d = sqrt((x2-x1)^2 + (y2-y1)^2);
%     m = (y2-y1) / (x2-x1);
%     r1 = sqrt((x-x1)^2 + (y-y1)^2 + z^2);
%     r2 = sqrt((x-x2)^2 + (y-y2)^2 + z^2);
%     e1 = (x-x1)^2 + z^2;
%     e2 = (x-x2)^2 + z^2;
%     h1 = (x-x1) * (y-y1);
%     h2 = (x-x2) * (y-y2);
%     
%     R = ((x-x1)*(y2-y1)-(y-y1)*(x2-x1))/d;
%     Q = log((r1+r2+d)/(r1+r2-d));
%     J = atan((m*e1-h1)/(z*r1)) - atan((m*e2-h2)/(z*r2));

%%
% [d12,d23,d31] = dFunc(verts);
% [m12,m23,m31] = mFunc(verts);
% [r1,r2,r3]    = rFunc(verts,x,y,z);
% [e1,e2,e3]    = eFunc(verts,x,z);
% [h1,h2,h3]    = hFunc(verts,x,y);
% 
% Phi = -(sigma/(4*pi)) * ...
%     (((((x-x1)*(y2-y1)-(y-y1)*(x2-x1))/d12) * log((r1+r2+d12)/(r1+r2-d12)) +...
%     (((x-x2)*(y3-y2)-(y-y2)*(x3-x2))/d23) * log((r2+r3+d23)/(r2+r3-d23)) +...
%     (((x-x3)*(y1-y3)-(y-y3)*(x1-x3))/d31) * log((r3+r1+d31)/(r3+r1-d31))) - abs(z)*...
%     (atan((m12*e1-h1)/(z*r1)) - atan((m12*e2-h2)/(z*r2)) +...
%     atan((m23*e2-h2)/(z*r2)) - atan((m23*e3-h3)/(z*r3)) +...
%     atan((m31*e3-h3)/(z*r3)) - atan((m31*e1-h1)/(z*r1))));