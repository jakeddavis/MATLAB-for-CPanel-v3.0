
% Jake Davis
% 03/02/2018

function Phi = subLinearDoubletFarField(muVec,verts,A,Q,P)

mu_0 = muVec(1); mu_x = muVec(2); mu_y = muVec(3);
x0 = Q(1); y0 = Q(2);
x = P(1);  y = P(2);  z = P(3);
x1 = verts(1,1); y1 = verts(1,2);
x2 = verts(2,1); y2 = verts(2,2);
x3 = verts(3,1); y3 = verts(3,2);

mu1 = mu_0 + mu_x*x1 + mu_y*y1;
mu2 = mu_0 + mu_x*x2 + mu_y*y2;
mu3 = mu_0 + mu_x*x3 + mu_y*y3;
% mu = (A/3)*(mu1 + mu2 + mu3);
mu = (A/3)*(mu1 + mu2 + mu3);

r = ((x-x0)^2 + (y-y0)^2 + z^2)^(-3/2);

Phi = -(mu/4/pi) * (z*r);

end

%%

% % x = P_vec(1); y = P_vec(2);
% mu_0 = muVec(1); mu_x = muVec(2); mu_y = muVec(3);
% verts(4,:) = verts(1,:);
% 
% M = 1;
% N = 1;
% I_11 = I_MN_FarField(verts,P_vec,M,N);
% 
% M = 2;
% N = 1;
% I_21 = I_MN_FarField(verts,P_vec,M,N);
% 
% M = 1;
% N = 2;
% I_12 = I_MN_FarField(verts,P_vec,M,N);
% 
% % mu = mu_0 + mu_x*x + mu_y*y;
% Phi = mu_0*I_11 + mu_x*I_21 + mu_y*I_12;
% % Phi = mu_0*I_11 + mu_x*I_21 + mu_y*I_12;


%%

% G_11 = zeros(1,3);
% for i = 1:3
%     pnt1 = [verts(i,1) verts(i,2)];
%     pnt2 = [verts(i+1,1) verts(i+1,2)];
% %     x1 = pnt1(1); y1 = pnt1(2);
% %     x2 = pnt2(1); y2 = pnt2(2);
%     
%     triGeom = triGeom_Johnson(pnt1,pnt2,P_vec);
%     
% %     M = 1; N = 1;
%     G_11(i) = G_MN_func(pnt1,pnt2,M,N,triGeom);
% %     
% %     if abs(nu_eta) <= abs(nu_xi)
% %         E_120 = (x2^(1-1)*y2^(2-1)) - (x1^(1-1)*y1^(2-1));
% %         D_12 = E_120;
% %         G_11(i) = (1/(1*nu_xi)) * D_12;
% %     else
% %         E_210 = (x2^(2-1)*y2^(1-1)) - (x1^(2-1)*y1^(1-1));
% %         D_21 = E_210;
% %         G_11(i) = -(1/(1*nu_eta)) * D_21;
% %     end
%     
% end    
% 
% a = triGeom.a;
% C_11 = (1/(1+1)) * a*sum(G_11);
% E1 = C_11;

% C_MN_func = 
% C_11 = (1/(1+1)) * a*sum(G_11(i));
% E1 = C_11;
% 
% C_21 = 
% E2 = [C_21 C_12];
%     
% J_vec_11 = -(1/4/pi) * (-(1/P^2)*(E1*P_car) + (1/P^3)*(E2-3*(E2.*P_car)*P_car)...
%     + (1/P^4) * ((3/2)*E3*P_car - (15/2)*(P_car.*E4.*P_car).*P_car + 3*E4.*P_car));

% I_11 = J_vec_11;
% I_21 = 0;
% I_12 = 0;