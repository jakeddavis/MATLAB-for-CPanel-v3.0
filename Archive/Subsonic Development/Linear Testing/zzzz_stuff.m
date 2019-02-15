
% Jake Davis
% 02/19/2018

clc, clear, close all

% build triangle
verts = [0 0 ; 1 0 ; cosd(60) sind(60)];
[A,Q,diam] = triParams(verts);

% define 'influenced' point
P = [.75 .8 5];

% move panel local coordinate origin to center of triangle, in x-y plane
verts(:,1) = verts(:,1) - Q(1);
verts(:,2) = verts(:,2) - Q(2);
P(1) = P(1) - Q(1);
P(2) = P(2) - Q(2);
% recompute triparams
[~,Q,~] = triParams(verts);
% replot
figure
hold on
vertIDs = [1 2 3];
tri = triangulation(vertIDs,verts);
triplot(tri,'-ob')
plot(Q(1),Q(2),'*m')
plot(P(1),P(2),'*r')
axis equal
% hold off

x0 = P(1); y0 = P(2); z0 = P(3);
x1 = verts(1,1); y1 = verts(1,2);
x2 = verts(2,1); y2 = verts(2,2);
x3 = verts(3,1); y3 = verts(3,2);
verts_wrap = verts;
verts_wrap(4,:) = verts_wrap(1,:);

%%

d12 = sqrt((x2-x1)^2+(y2-y1^2));

C12 = (x2-x1)/d12;
S12 = (y2-y1)/d12;

s12_1 = (x1-x0)*C12 + (y1-y0)*S12;
s12_2 = (x2-x0)*C12 + (y2-y0)*S12;

R12 = (x0 - x1)*S12 - (y0-y1)*C12;

r1 = sqrt((x0-x1)^2 + (y0-y1)^2 + z0^2);
r2 = sqrt((x0-x2)^2 + (y0-y2)^2 + z0^2);

J12 = atan2((R12*abs(z0)*(r1*s12_1 - r2*s12_2)), ...
    (r1*r2*R12^2 + z0^2*s12_2*s12_1));

J12 = J12/4/pi

%% Starting Stuff

% h = P(3);
% F_111_sum = 0;
% H_113_sum = 0;
% H_213_sum = 0;
% H_123_sum = 0;
% for i = 1:3
%     [point(i,:),a(i),l1(i),l2(i),c1(i),c2(i),g(i),nu_xi(i),nu_eta(i)] = ...
%         triEdge_fieldPoint_geom(verts_wrap(i,:),verts_wrap(i+1,:),P);
%     plot(point(i,1),point(i,2),'*k')
%     
% %     F_111 = log((sqrt(l2(i)^2+g(i)^2)+l2(i))/(sqrt(l1(i)^2+g(i)^2)+l1(i)));
% %     F_111_sum = F_111_sum + a(i)*F_111;
% %     
% %     H_113_sum = H_113_sum + atan2(a(i) * ...
% %         (l2(i)*c1(i) - l1(i)*c2(i)),c1(i)*c2(i) + a(i)^2*l1(i)*l2(i));
%     
%     H_113_sum = H_113_sum + atan2(a(i) * ...
%         (l2(i)*c1(i) - l1(i)*c2(i)),c1(i)*c2(i) + a(i)^2*l1(i)*l2(i));
%     
%     F_111 = log((sqrt(l2(i)^2+g(i)^2)+l2(i))/(sqrt(l1(i)^2+g(i)^2)+l1(i)));
%     
%     H_213_sum = H_213_sum + nu_xi(i) * F_111;
% 
%     H_123_sum = H_123_sum + nu_eta(i) * F_111;
% end
% % H_111 = -abs(h) * H_113_sum + F_111_sum;
% % H_113_step2 = (-H_111 + F_111_sum) / h^2
% 
% H_113 = H_113_sum/abs(h);
% H_213 = -H_213_sum;
% H_123 = -H_123_sum;
% 
% I_11 = h * H_113 / (4*pi);
% I_21 = h * H_213 / (4*pi);
% I_12 = h * H_123 / (4*pi);
% 
% mu_0 = 1; mu_x = 1; mu_y = 1;
% mu = mu_0 + mu_x*x0 + mu_y*y0;
% 
% Phi = mu * I_11 + mu_x * I_21 + mu_y * I_12


%%
% % plot geometry
% figure
% hold on
% vertIDs = [1 2 3];
% tri = triangulation(vertIDs,verts);
% triplot(tri,'-ob')
% plot(Q(1),Q(2),'*m')
% plot(P(1),P(2),'*r')
% axis equal
% hold off

%%
% % compute sides and angles of triangle
% side12 = sqrt((x1-x2)^2 + (y1-y2)^2);
% side23 = sqrt((x2-x3)^2 + (y2-y3)^2);
% side31 = sqrt((x3-x1)^2 + (y3-y1)^2);
% 
% angle12 = acosd((-side12^2 + side23^2 + side31^2)/(2*side23*side31));
% angle23 = acosd((side12^2 - side23^2 + side31^2)/(2*side12*side31));
% angle31 = acosd((side12^2 + side23^2 - side31^2)/(2*side12*side23));
% 
% % compute a
% a12 = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)/sqrt((y2-y1)^2 + (x2-x1)^2);
% a23 = abs((y3-y2)*x0 - (x3-x2)*y0 + x3*y2 - y3*x2)/sqrt((y3-y2)^2 + (x3-x2)^2);
% a31 = abs((y1-y3)*x0 - (x1-x3)*y0 + x1*y3 - y1*x3)/sqrt((y1-y3)^2 + (x1-x3)^2);
% 
% m = (y1-y3)/(x1-x3);
% k = y3 - m*x3;
% % plot(0,k,'*c')
% x = (x0 + m*y0 - m*k) / (m^2 + 1);
% y = m*x + k;
% plot(x,y,'*k')







