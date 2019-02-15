
% Jake Davis
% 05/20/2018

function Phi = CPanelStuff_0520(verts,POI,mu,cntr)

x = POI(1); y = POI(2); 
% z = POI(3);
% verts = [verts(1,:) 0; verts(2,:) 0; verts(3,:) 0];
verts(4,:) = verts(1,:);

% get normal
normal = cross(verts(2,:)-verts(1,:), verts(3,:)-verts(1,:));
normal = normal/norm(normal);

% get panel center
% [~,Q,~] = triParams(verts,POI);
% cntr = [Q 0];

pjk = POI - cntr;
PN = dot(pjk,normal);

l = (verts(1,:) - cntr)/norm(verts(1,:) - cntr);    % panel x-dir
m = cross(normal,l); % panel y-dir
n = normal; % panel z-dir

% % Plot coord. system
% % figure
% hold all
% plot3([0 l(1)], [0 l(2)], [0 l(3)], 'k')
% plot3([0 m(1)], [0 m(2)], [0 m(3)], 'k')
% plot3([0 n(1)], [0 n(2)], [0 n(3)], 'k')

Phi = 0;
for i = 1:3
    
    p1 = verts(i,:);
    p2 = verts(i+1,:);
    
    a = POI - p1;
    b = POI - p2;
    s = p2 - p1;
    
    Al = dot(n,cross(s,a));
%     Al = -dot(n,cross(s,a));
% %     sa = cross(s,a);
% %     plot3(sa(1),sa(2),sa(3),'*m')
%     thing = cross(l,cross(a,s));
% %     plot3(thing(1), thing(2), thing(3), '*g')
%     
%     nu_xi = dot(s,l)/norm(s);
%     nu_eta = dot(s,m)/norm(s);
% %     l1 = dot(-a,m)*dot(s,m) + dot(-a,l)*dot(s,l);
% %     l2 = dot(-b,m)*dot(s,m) + dot(-b,l)*dot(s,l);
%     l1 = dot(-a,m)*nu_eta + dot(-a,l)*nu_xi;
%     l2 = dot(-b,m)*nu_eta + dot(-b,l)*nu_xi;
% %     nu_xi = -dot(s,l)/norm(s);
% %     nu_eta = dot(s,m)/norm(s);
%     g = sqrt(Al^2+PN^2);
%     
%     triGeom = triGeom_func(p1,p2,POI);
%     a_mine = triGeom.a;
%     l1_m = triGeom.l1; l2_m = triGeom.l2;
%     c1 = triGeom.c1; c2 = triGeom.c2;
%     num_m = a_mine*(l2_m*c1-l1_m*c2);
%     denom_m = c1*c2+a_mine^2*l1_m*l2_m;
%     nu_xi = triGeom.nu_xi; nu_eta = triGeom.nu_eta;
    
%     c1 = triGeom.c1; c2 = triGeom.c2;
    dotsm = dot(s,m);
    PA = dot(a,cross(l,cross(a,s)));
    PB = PA-Al*dot(s,m);
    num = (dot(s,m)*PN *(norm(b)*PA - norm(a)*PB));
    denom = PA*PB + PN^2*norm(a)*norm(b)*dot(s,m)^2;
    phiV = atan2(num,denom);
    
    A = norm(a);
    B = norm(b);
    S = norm(s);
    GL = 1/S * log(abs((A+B+S)/(A+B-S)));
    
    Phi = Phi + Al*GL - PN*phiV;
    
%     H_113_mine = H_113_mine + atan2(num_m,denom_m)
    
%     if l1 >= 0 && l2 >= 0
%         F_111 = log((sqrt(l2^2+g^2)+l2) / (sqrt(l1^2+g^2)+l1));
%     elseif l1 < 0 && l2 < 0
%         F_111 = log((sqrt(l1^2+g^2)-l1) / (sqrt(l2^2+g^2)-l2));
%     elseif l2 >= 0 && l1 < 0
%         F_111 = log(((sqrt(l1^2+g^2)-l1) * (sqrt(l2^2+g^2)+l2)) / g^2);
%     end
% %     disp(F_111)
% 
%     H_113_CP_tan   = H_113_CP_tan + atan2(num,denom);
%     H_213_CP = H_213_CP + F_111 * nu_xi;
%     H_123_CP = H_123_CP + F_111 * nu_eta;
    
%     phi = phi + atan2(num,denom);
end

% % H_113_CP = (1/PN)*H_113_CP_tan
% H_113_CP = (1/PN)*H_113_CP_tan;
% H_213_CP = -H_213_CP;
% H_123_CP = -H_123_CP;
% 
% I_11_CP = PN * H_113_CP / (4*pi);
% I_21_CP = PN * H_213_CP / (4*pi);
% I_12_CP = PN * H_123_CP / (4*pi);
% 
% mu_0 = mu(1); mu_x = mu(2); mu_y = mu(3);
% mu = mu_0 + mu_x*x + mu_y*y;
% 
% Phi = mu*I_11_CP + mu_x*I_21_CP + mu_y*I_12_CP;

% phi = sum(phi)/(4*pi)

end
