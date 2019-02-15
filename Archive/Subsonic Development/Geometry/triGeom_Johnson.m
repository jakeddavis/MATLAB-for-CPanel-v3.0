
% Jake Davis
% 02/28/2018

% function [H_113,F_111,nu_xi,nu_eta] = triGeom_Johnson(pnt1,pnt2,P)
% function [a,g,l1,l2,c1,c2,nu_xi,nu_eta] = triGeom_Johnson(pnt1,pnt2,P)
function triGeom = triGeom_Johnson(pnt1,pnt2,P)

x = P(1); y = P(2); z = P(3);
x1 = pnt1(1); y1 = pnt1(2);
x2 = pnt2(1); y2 = pnt2(2);

d = sqrt((x2-x1)^2 + (y2-y1)^2);
C = (x2-x1) / d;
S = (y2-y1) / d;

a = (x-x1)*S - (y-y1)*C;
h = z;
g = sqrt(a^2 + h^2);

l1 = (x1-x)*C + (y1-y)*S;
l2 = (x2-x)*C + (y2-y)*S;

s1 = sqrt(l1^2 + g^2);
s2 = sqrt(l2^2 + g^2);
c1 = g^2 + abs(h)*s1;
c2 = g^2 + abs(h)*s2;

nu_xi  = -(y2-y1) / (l2-l1);
nu_eta = (x2-x1) / (l2-l1);

%%%%%%%%%%%%%%%%%%%
% x0 = P(1); y0 = P(2); h = P(3);
% x1 = V1(1); y1 = V1(2);
% x2 = V2(1); y2 = V2(2);
% 
% h = z;
% m = (y2-y1)/(x2-x1);
% k = y1 - m*x1;
% x_int = (x + m*y - m*k) / (m^2 + 1);
% y_int = m*x_int + k;
% % int_pnt = [x_int y_int];
% 
% % a = abs((y2-y1)*x - (x2-x1)*y + x2*y1 - y2*x1)/sqrt((y2-y1)^2 + (x2-x1)^2);
% % if x > x_int
% %     a = -a;
% % end
% 
% d = sqrt((x2-x1)^2 + (y2-y1)^2);
% C = (x2-x1) / d;
% S = (y2-y1) / d;
% l1 = (x1-x)*C + (y1-y)*S;
% l2 = (x2-x)*C + (y2-y)*S;
% % l1 = sqrt((x1-x)^2 + (y1-y)^2);
% % l2 = sqrt((x2-x)^2 + (y2-y)^2);
% 
% a = abs((y2-y1)*x - (x2-x1)*y + x2*y1 - y2*x1)/sqrt((y2-y1)^2 + (x2-x1)^2);
% phi1 = atan2d(abs(l1),a);
% if phi1 > 90
%     a = -a;
% end
% 
% g = sqrt(a^2 + h^2);
% c1 = g^2 + abs(h)*sqrt(l1^2 + g^2);
% c2 = g^2 + abs(h)*sqrt(l2^2 + g^2);
% 
% a_nu_xi  = x_int - x;
% a_nu_eta = y_int - y;
% 
% % a_nu_xi  = abs(x - x0);
% % a_nu_eta = abs(y - y0);
% 
% nu_xi  = a_nu_xi / (a);
% nu_eta = a_nu_eta / (a);

%%%%%%%%%%%%%%%%
% if x <= x_int
%     nu_xi  = abs(y2-y1) / (l2-l1);
% else
%     nu_xi  = -abs(y2-y1) / (l2-l1);
% end
% 
% if y <= y_int
%     nu_eta = abs(x2-x1) / (l2-l1);
% else
%     nu_eta = -abs(x2-x1) / (l2-l1);
% end

triGeom.a = a; triGeom.g = g; triGeom.l1 = l1; triGeom.l2 = l2;
triGeom.c1 = c1; triGeom.c2 = c2; triGeom.nu_xi = nu_xi; triGeom.nu_eta = nu_eta;

end

%%
% H_113 = atan2(a*(l2*c1-l1*c2) , c1*c2+a^2*l1*l2);
% 
% if l1 >= 0 && l2 >= 0
%     F_111 = log((sqrt(l2^2+g^2)+l2) / (sqrt(l1^2+g^2)+l1));
% elseif l1 < 0 && l2 < 0
%     F_111 = log((sqrt(l1^2+g^2)-l1) / (sqrt(l2^2+g^2)-l2));
% elseif l2 >= 0 && l1 < 0
%     F_111 = log(((sqrt(l1^2+g^2)-l1) * (sqrt(l2^2+g^2)+l2))/g^2);
% end
