
% Jake Davis
% 05/21/2018

function [u, v, w] = subLinearDoubletNearFieldVel(mu,verts,P)

x = P(1); y = P(2); h = P(3);
verts(4,:) = verts(1,:);

F113 = 0;
F123 = 0;
H113 = 0;
H213 = 0;
H123 = 0;

H115 = 0;
H215 = 0;
H125 = 0;

H225 = 0;

for i = 1:3
    pnt1 = [verts(i,1) verts(i,2)];
    pnt2 = [verts(i+1,1) verts(i+1,2)];
    
    pnt1_3d = [pnt1 0];
    pnt2_3d = [pnt2 0];
    
    triGeom = triGeom_func(pnt1,pnt2,P);
    a  = triGeom.a;  g = triGeom.g;
    l1 = triGeom.l1; l2 = triGeom.l2;
    c1 = triGeom.c1; c2 = triGeom.c2;
    nu_xi = triGeom.nu_xi; nu_eta = triGeom.nu_eta;

    if l1 >= 0 && l2 >= 0
        F111 = log((sqrt(l2^2+g^2)+l2) / (sqrt(l1^2+g^2)+l1));
    elseif l1 < 0 && l2 < 0
        F111 = log((sqrt(l1^2+g^2)-l1) / (sqrt(l2^2+g^2)-l2));
    elseif l2 >= 0 && l1 < 0
        F111 = log(((sqrt(l1^2+g^2)-l1) * (sqrt(l2^2+g^2)+l2)) / g^2);
    end
    
    rho1 = norm(pnt1_3d - P);
    rho2 = norm(pnt2_3d - P);
    test1 = norm(P - pnt1_3d);
    test2 = norm(P - pnt2_3d);
    E111 = 1/rho2 - 1/rho1;
    E211 = (pnt2(1)-x)/rho2 - (pnt1(1)-x)/rho1;
    E121 = (pnt2(2)-y)/rho2 - (pnt1(2)-y)/rho1;
    
    F113_temp = (1/g^2) * (-nu_eta*E211 + nu_xi*E121);
    F113 = F113 + F113_temp;
    F123_temp = nu_eta*a*F113_temp - nu_xi*E111;
    F123 = F123 + nu_eta*F123_temp;
    
    num_mine = a*(l2*c1-l1*c2);
    denom_mine = c1*c2+a^2*l1*l2;
    H113 = H113 + atan2(num_mine , denom_mine);
    H213 = H213 + F111 * nu_xi;
    H123 = H123 + F111 * nu_eta;

    H115 = H115 + a*F113_temp;
    H215 = H215 + nu_xi*F113_temp;
    H125 = H125 + nu_eta*F113_temp;
    
    H225 = H225 + nu_xi*F123_temp;
end

% H integrals
H113 = (1/h) * H113;
H213 = -H213;
H123 = -H123;

H115 = (1/3/h^2) * (H113 + H115);
H215 = (1/3) * -H215;
H125 = (1/3) * -H125;

H225 = (1/3) * -H225;
H135 = (1/3) * (-H113 - F123);
H315 = -H135 - h^2*H115 + H113;

% Jx = (1/4/pi) * 3*h * H(M+1,N,5)
J11_x = 3*h * H215;
J21_x = -(1/h)*H113 + 3*h * H315;
J12_x = 3*h * H225;

% Jy = (1/4/pi) * 3*h * H(M,N+1,5)
J11_y = 3*h * H125;
J21_y = 3*h * H225;
% J12_y = -(h^2)*H113 + 3*h * H135_y;
J12_y = 3*h * H135;

% Jz = (1/4/pi) * (H(M,N,3) - 3*h^2*H(M,N,5))
J11_z = H113 - 3*h^2*H115;
J21_z = -H213 - 3*h^2*H215;
J12_z = -H123 - 3*h^2*H125;

% Collect J's
J11 = (1/4/pi) * [J11_x J11_y J11_z];
J21 = (1/4/pi) * [J21_x J21_y J21_z];
J12 = (1/4/pi) * [J12_x J12_y J12_z];

% get mu's
mu_0 = mu(1); mu_x = mu(2); mu_y = mu(3);
mu = mu_0 + mu_x*x + mu_y*y;

% Compute velocity vector
V = mu*J11 + mu_x*J21 + mu_y*J12;
u = V(1);
v = V(2);
w = V(3);

end
