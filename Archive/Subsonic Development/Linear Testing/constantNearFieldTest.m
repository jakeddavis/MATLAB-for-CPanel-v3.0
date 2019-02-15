
function Phi = constantNearFieldTest(sigma,verts,P)

verts(4,:) = verts(1,:);
x = P(1); y = P(2); z = P(3);

Phi = zeros(1,3);
for i = 1:3
    % tri edge pnt 1
    x1 = verts(i,1); y1 = verts(i,2);
    % tri edge pnt 2
    x2 = verts(i+1,1); y2 = verts(i+1,2);
    
    d = sqrt((x2-x1)^2+(y2-y1)^2);
    
    C = (x2-x1)/d;
    S = (y2-y1)/d;
    
    s1 = (x1-x)*C + (y1-y)*S;
    s2 = (x2-x)*C + (y2-y)*S;
    
    R = ((x-x1)*(y2-y1)-(y-y1)*(x2-x1))/d;
    R_hs = (x-x1)*S - (y-y1)*C
    R = (x-x1)*S - (y-y1)*C;
    
    r1 = sqrt((x-x1)^2 + (y-y1)^2 + z^2);
    r2 = sqrt((x-x2)^2 + (y-y2)^2 + z^2);
    
    Q = log((r1+r2+d)/(r1+r2-d))
    Q_hs = log((r1+r2+d)/(r1+r2-d))
%     J = atan((m*e1-h1)/(z*r1)) - atan((m*e2-h2)/(z*r2));
    m = (y2-y1) / (x2-x1);
    e1 = (x-x1)^2 + z^2;
    e2 = (x-x2)^2 + z^2;
    h1 = (x-x1) * (y-y1);
    h2 = (x-x2) * (y-y2);
    J = atan((m*e1-h1)/(z*r1)) - atan((m*e2-h2)/(z*r2));
    
    Y1 = max(y1,y2);
    Y2 = min(y1,y2);
%     gamma = acosd((Y1-Y2)/sqrt((x1-x2)^2+(Y1-Y2)^2))
%     
%     gamma2 = 90 + atand((y2-y1)/(x2-x1))
    J_hs_2 = atan2((R*abs(z)*(r1*s2 - r2*s1)) , (r1*r2*R^2 + z^2*s2*s1))
    J_hs_1 = sign(R) * (atan(abs(z/R)*(s2/r2)) - atan(abs(z/R)*(s1/r1)))
    a = R; h = z;
    l1 = s1; l2 = s2;
    g = sqrt(a^2 + h^2);
    s1_j = sqrt(l1^2 + g^2);
    s2_j = sqrt(l2^2 + g^2);
    c1 = g^2 + abs(h)*s1_j;
    c2 = g^2 + abs(h)*s2_j;
    
    J_john_2 = atan2(a*(l2*c1-l1*c2) , c1*c2+a^2*l1*l2)
    J_john_1 = atan((a*l2)/c2) - atan((a*l1)/c1)
    
    Phi(i) = R*Q - abs(z)*J;
    
%     J12_sum = -J12_sum - atan2((R*abs(z)*(r1*s2 - r2*s1)), ...
%         (r1*r2*R^2 + z^2*s2*s1));
    
%     d = sqrt((x2-x1)^2 + (y2-y1)^2);
%     m = (y2-y1) / (x2-x1);
%     r1 = sqrt((x-x1)^2 + (y-y1)^2 + z^2);
%     r2 = sqrt((x-x2)^2 + (y-y2)^2 + z^2);
%     e1 = (x-x1)^2 + z^2;
%     e2 = (x-x2)^2 + z^2;
%     h1 = (x-x1) * (y-y1);
%     h2 = (x-x2) * (y-y2);
%     
%     R12 = ((x-x1)*(y2-y1)-(y-y1)*(x2-x1))/d;
%     Q12 = log((r1+r2+d)/(r1+r2-d));
%     J12 = atan((m*e1-h1)/(z*r1)) - atan((m*e2-h2)/(z*r2));
%     
%     Phi(i) = R12*Q12 - abs(z)*J12;
    
end

Phi = (-sigma/4/pi)*sum(Phi);

% verts(4,:) = verts(1,:);
% x0 = P(1); y0 = P(2); z0 = P(3);
% J12_sum = 0;
% 
% for i = 1:3
% 
%     x1 = verts(i,1); y1 = verts(i,2);
%     x2 = verts(i+1,1); y2 = verts(i+1,2);
%     
%     d12 = sqrt((x2-x1)^2+(y2-y1^2));
%     
%     C12 = (x2-x1)/d12;
%     S12 = (y2-y1)/d12;
%     
%     s12_1 = (x1-x0)*C12 + (y1-y0)*S12;
%     s12_2 = (x2-x0)*C12 + (y2-y0)*S12;
%     
%     R12 = (x0 - x1)*S12 - (y0-y1)*C12;
%     
%     r1 = sqrt((x0-x1)^2 + (y0-y1)^2 + z0^2);
%     r2 = sqrt((x0-x2)^2 + (y0-y2)^2 + z0^2);
%     
%     J12_sum = -J12_sum - atan2((R12*abs(z0)*(r1*s12_2 - r2*s12_1)), ...
%         (r1*r2*R12^2 + z0^2*s12_2*s12_1));
% 
% end
% 
% Phi = mu*J12_sum/4/pi;

end