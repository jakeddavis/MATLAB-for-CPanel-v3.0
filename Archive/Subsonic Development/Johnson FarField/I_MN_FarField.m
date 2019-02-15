
% Jake Davis
% 03/02/2018

function I_MN = I_MN_FarField(verts,P_vec,M,N)

P = norm(P_vec);
P_car = (P_vec/P)';

% J from source equations
E1 = C_MN_func(M,N,verts,P_vec);
E2 = [C_MN_func(M+1,N,verts,P_vec) C_MN_func(M,N+1,verts,P_vec) 0]';
E3 = C_MN_func(M+2,N,verts,P_vec) + C_MN_func(M,N+2,verts,P_vec);
E4 = [C_MN_func(M+2,N,verts,P_vec) C_MN_func(M+1,N+1,verts,P_vec) 0 ;...
      C_MN_func(M+1,N+1,verts,P_vec) C_MN_func(M,N+2,verts,P_vec) 0 ;...
      0 0 0];
J_vec_MN = -(1/4/pi) * (-(1/P^2)*(E1*P_car) + (1/P^3)*(E2-3*(dot(E2,P_car)).*P_car)...
    + (1/P^4) * ((3/2)*E3*P_car - (15/2)*(dot(P_car,E4*P_car))*P_car + 3*E4*P_car));

I_MN = J_vec_MN(3)

% I from doublet equations
E2 = [0 0 C_MN_func(M,N,verts,P_vec)];
E3 = 0;
E4 = [0 0 (1/2)*C_MN_func(M+1,N,verts,P_vec) ;...
      0 0 (1/2)*C_MN_func(M,N+1,verts,P_vec) ;...
      (1/2)*C_MN_func(M+1,N,verts,P_vec) (1/2)*C_MN_func(M,N+1,verts,P_vec) 0];
  
I_MN_J = (1/4/pi) * ((1/P^2)*dot(E2,P_car) + (1/P^3)*(-E3+3*(dot(P_car,E4*P_car))))

end

%%

% G_MN = zeros(1,3);
% for i = 1:3
%     
%         pnt1 = [verts(i,1) verts(i,2)];
%         pnt2 = [verts(i+1,1) verts(i+1,2)];
% %         x1 = pnt1(1); y1 = pnt1(2);
% %         x2 = pnt2(1); y2 = pnt2(2);
%     
%         triGeom = triGeom_Johnson(pnt1,pnt2,P_vec);
%         
%         G_MN(i) = G_MN_func(pnt1,pnt2,M,N,triGeom);
%     
% end

% a = triGeom.a;
% C_MN = C_MN_func(M,N,G_MN,a)
% C_MN = (1/(M+N)) * a*sum(G_MN);