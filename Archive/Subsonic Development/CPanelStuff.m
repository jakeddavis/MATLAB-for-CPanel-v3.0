
% Jake Davis
% 05/04/2018

function phi = CPanelStuff(verts,POI)

% verts3D = [verts(1,:) 0; verts(2,:) 0; verts(3,:) 0];
verts3D = [verts(1,:) 0; verts(2,:) 0; verts(3,:) 0];
verts3D(4,:) = verts3D(1,:);

% get normal
normal = cross(verts3D(2,:)-verts3D(1,:), verts3D(3,:)-verts3D(1,:));
normal = normal/norm(normal);

% get panel center
[~,Q,~] = triParams(verts,POI);
cntr = [Q 0];

pjk = POI - cntr;
PN = dot(pjk,normal);

l = (verts3D(1,:) - cntr)/norm(verts3D(1,:) - cntr);    % panel x-dir
m = cross(l,normal); % panel y-dir
n = normal; % panel z-dir

% Plot coord. system
% figure
hold all
plot3([0 l(1)], [0 l(2)], [0 l(3)], 'k')
plot3([0 m(1)], [0 m(2)], [0 m(3)], 'k')
plot3([0 n(1)], [0 n(2)], [0 n(3)], 'k')

phi = 0;
H_113_mine = 0;
H_113_CP = 0;
for i = 1:3
    
    p1 = verts3D(i,:);
    p2 = verts3D(i+1,:);
    
    a = POI - p1;
    b = POI - p2;
    s = p2 - p1;
    
%     Al = dot(n,cross(s,a));
    Al = dot(n,cross(s,a));
    sa = cross(s,a);
    plot3(sa(1),sa(2),sa(3),'*m')
    
    l1 = dot(-a,m)*dot(s,m) + dot(-a,l)*dot(s,l);
    l2 = dot(-b,m)*dot(s,m) + dot(-b,l)*dot(s,l);
    nu_xi = -dot(s,l)/norm(s);
    nu_eta = dot(s,m)/norm(s);
    g = sqrt(Al^2+PN^2);
    
    triGeom = triGeom_func(p1,p2,POI);
    a_mine = triGeom.a;
    l1_m = triGeom.l1; l2_m = triGeom.l2;
    c1 = triGeom.c1; c2 = triGeom.c2;
    num_m = a_mine*(l2_m*c1-l1_m*c2);
    denom_m = c1*c2+a_mine^2*l1_m*l2_m;
%     nu_xi = triGeom.nu_xi; nu_eta = triGeom.nu_eta;
    
%     c1 = triGeom.c1; c2 = triGeom.c2;
    dotsm = dot(s,m)
    PA = dot(a,cross(l,cross(a,s)));
    PB = PA-Al*dot(s,m);
    num = (dot(s,m)*PN *(norm(b)*PA - norm(a)*PB));
    denom = PA*PB + PN^2*norm(a)*norm(b)*dot(s,m)^2;
    
    H_113_mine = H_113_mine + atan2(num_m,denom_m)
    H_113_CP   = H_113_CP + atan2(num,denom)
    
    if l1 >= 0 && l2 >= 0
        F_111 = log((sqrt(l2^2+g^2)+l2) / (sqrt(l1^2+g^2)+l1));
    elseif l1 < 0 && l2 < 0
        F_111 = log((sqrt(l1^2+g^2)-l1) / (sqrt(l2^2+g^2)-l2));
    elseif l2 >= 0 && l1 < 0
        F_111 = log(((sqrt(l1^2+g^2)-l1) * (sqrt(l2^2+g^2)+l2)) / g^2);
    end
%     disp(F_111)
    
    phi = phi + atan2(num,denom);

end

phi = sum(phi)/(4*pi)
