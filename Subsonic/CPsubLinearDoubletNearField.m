
% Jake Davis
% 02/25/2018

function Phi = CPsubLinearDoubletNearField(mu,verts,P,cond,local)

POI = P;
x = P(1); y = P(2); h = P(3);
verts(4,:) = verts(1,:);
in_tri = cond(2);

cntr = [mean(verts(:,1)),mean(verts(:,2)),mean(verts(:,3))];
pjk = POI - cntr;
% POIloc = ...

H_113_tan = zeros(1,3);
H_213 = zeros(1,3);
H_123 = zeros(1,3);

for i = 1:3
    pnt1 = verts(i,:);
    pnt2 = verts(i+1,:);
    
    PN = dot(pjk,local(3,:));
    a = POI - pnt1;
    b = POI - pnt2;
    s = pnt2 - pnt1;
    Al = dot(local(3,:),cross(s,a));
    
    nuXi = dot(s,local(2,:)) / norm(s);
    nuEta = dot(s,local(1,:)) / norm(s);
    
    l = local(1,:); m = local(2,:);
    PA = dot(a,cross(l,cross(a,s)));
    PB = PA - Al*dot(s,m);
    numCP = dot(s,m)*PN*(norm(b)*PA - norm(a)*PB);
    denomCP = PA*PB+PN^2*norm(a)*norm(b)*dot(s,m)^2;
    H113CP = atan2(numCP,denomCP);
    
    norm(s)
    
    triGeom = triGeom_func(pnt1,pnt2,P);
    a  = triGeom.a;  g = triGeom.g;
    l1 = triGeom.l1; l2 = triGeom.l2;
    c1 = triGeom.c1; c2 = triGeom.c2;
    nu_xi = triGeom.nu_xi; nu_eta = triGeom.nu_eta;

    g = sqrt((a)^2 + (PN)^2);
    s1 = sqrt(l1^2 + g^2);
    s2 = sqrt(l2^2 + g^2);
    c1 = g^2 + abs(PN)*s1;
    c2 = g^2 + abs(PN)*s2;
    
    if l1 >= 0 && l2 >= 0
        F_111 = log((sqrt(l2^2+g^2)+l2) / (sqrt(l1^2+g^2)+l1));
    elseif l1 < 0 && l2 < 0
        F_111 = log((sqrt(l1^2+g^2)-l1) / (sqrt(l2^2+g^2)-l2));
    elseif l2 >= 0 && l1 < 0
        F_111 = log(((sqrt(l1^2+g^2)-l1) * (sqrt(l2^2+g^2)+l2)) / g^2);
    end
%     disp(F_111)
    num_mine = a*(l2*c1-l1*c2);
    denom_mine = c1*c2+a^2*l1*l2;
    H_113_tan(i) = atan2(num_mine , denom_mine);
    H_213(i) = F_111 * nu_xi;
    H_123(i) = F_111 * nu_eta;
    
    H111(i) = -abs(PN) * H_113_tan(i) + a*F_111;
    aF111(i) = a * F_111;
  
end

H_113 = (1/PN^2) * (-sum(H111) + sum(aF111));

% Check field point and compute H_113
% if (abs(h) == 0) && (in_tri == 0)
%     H_113 = 0;
% elseif (abs(h) == 0) && (in_tri == 1)
%     H_113 = (2*pi) / abs(h);
% else
%     H_111 = -abs(h)*sum(H_113_tan) + a*sum(F_111);
%     H_113 = (1/h^2) * (-H_111 + a*sum(F_111));
%     H_113_noF = (1/h)*sum(H_113_tan);
% end

H_213 = -sum(H_213);
H_123 = -sum(H_123);

% I_11 = h * H_113 / (4*pi);
% I_21 = h * H_213 / (4*pi);
% I_12 = h * H_123 / (4*pi);

I_11 = PN * H_113 / (4*pi);
I_21 = PN * H_213 / (4*pi);
I_12 = PN * H_123 / (4*pi);

mu_0 = mu(1); mu_x = mu(2); mu_y = mu(3);
mu = mu_0 + mu_x*x + mu_y*y;

Phi = mu*I_11 + mu_x*I_21 + mu_y*I_12;

end

%%


% % CONDITIONS
% 
% in_tri = cond(2); 
% dH = cond(3); delta_h = cond(4);
% dF = cond(5); delta_f = cond(6);
% 
% MXK = 3; NHK = 16;
% Kmax = MXK + NHK;
% K = 1:2:Kmax-2;
% F_11K = zeros(3,length(K));

% % Procedure 2
% 
% H_11Kstart = 0;
% sumF_11K = zeros(1,length(K));
% H_11K = zeros(1,length(K));
% H_11K(end) = H_11Kstart;
% 
% if (abs(h) < dH*delta_h) && (in_tri == 0)
%     for i = length(K):-1:1
%         if i == length(K)
%             sumF_11K(i) = sum(F_11K(:,i));
%             H_11K(i) = -a*sumF_11K(i);
%         else
%             sumF_11K(i) = sum(F_11K(:,i));
%             H_11K(i) = (1/(K(i)-4)) * (h^2*(K(i)-2)*H_11K(i+1) - a*sumF_11K(i));
%         end
%     end
% end
% 
% % ORIGINAL
% 
% H_111 = -abs(h)*sum(H_113_tan) + a*sum(F_111);
% 
% CONDITIONS
% 
% % Procedure 1
% if abs(h) >= dH*delta_h
%     H_111 = -abs(h)*sum(H_113_tan) + a*sum(F_111);
%     H_113 = (1/h^2) * (-H_111 + a*sum(F_111));
% % Procedure 2
% elseif (abs(h) < dH*delta_h) && (in_tri == 0)
%     H_113 = 0;
% elseif (abs(h) < dH*delta_h) && (in_tri == 1)
%     H_113 = (2*pi) / abs(h);
% end
% 
% % Procedure 1
% % if abs(h) >= dH*delta_h
%     H_111 = -abs(h)*sum(H_113_tan) + a*sum(F_111);
%     H_113 = (1/h^2) * (-H_111 + a*sum(F_111));
% Procedure 2

%     % CONDITIONS
%     
%     % Procedure 4
%     % THIS WILL ALWAYS BE TRUE IF CONTROL POINT IS SET OFF OF PANEL
%     if g >= delta_f*dF
%         if l1 >= 0 && l2 >= 0
%             F_111 = log((sqrt(l2^2+g^2)+l2) / (sqrt(l1^2+g^2)+l1));
%         elseif l1 < 0 && l2 < 0
%             F_111 = log((sqrt(l1^2+g^2)-l1) / (sqrt(l2^2+g^2)-l2));
%         elseif l2 >= 0 && l1 < 0
%             F_111 = log(((sqrt(l1^2+g^2)-l1) * (sqrt(l2^2+g^2)+l2))/g^2);
%         end
%     end
    
%     % Procedure 2
%     if (abs(h) < dH*delta_h) && (in_tri == 0)        
%         for j = 1:length(K)
%             if j == 1
%                 F_11K(i,j) = F_111;
%             else
%                 F_11K(i,j) = (1/(g^2*(K(j)-2))) * ((K(j)-3)*F_11K(i,j-1) - ...
%                 nu_eta*Efunc(2,1,K(j-1),pnt1,pnt2,P) + ...
%                 nu_xi*Efunc(1,2,K(j-1),pnt1,pnt2,P));
%             end
%         end
%     end
    

%     MXK = 3; NHK = 16;
%     Kstart = MXK + NHK;
%     H_11Kstart = 0;
%     K = Kstart:-2:3;
%     
%     F_11K = zeros(3,length(K));
%     for i = 1:3
%         pnt1 = [verts(i,1) verts(i,2)];
%         pnt2 = [verts(i+1,1) verts(i+1,2)];
%         for j = length(K):-1:2
%             F_11K(i,j) = (1/(g^2*(K(j)-2))) * ((K(j)-3)*F_11K(i,j+1) - ...
%                 nu_eta*Efunc(2,1,K(i,j+1),pnt1,pnt2,P) + ...
%                 nu_xi*Efunc(1,2,K(i,j+1),pnt1,pnt2,P));
%         end
%     end
%     
%     sumF_11K = zeros(1,length(K));
%     H_11K = zeros(1,length(K));
%     H_11K(1) = H_11Kstart;
%     for i = 2:length(K)
%         sumF_11K(i) = sum(F_11K(:,i));
%         H_11K(i) = (1/(K(i)-4)) * (h^2*(K(i)-2)*H_11K(i-1) - a*sumF_11K(i));
%     end
% end

%%
% % mini test
% mu_tx = [mu_x 0];
% for i = 1:2    
%     Phi_x(i) = mu_tx(i)*I_21;
% end
% DPhi_x = (Phi_x(2) - Phi_x(1)) / (mu_tx(2) - mu_tx(1));
% I_21_x = DPhi_x - mu_tx*I_11;
% 
% mu_ty = [mu_y 0];
% for i = 1:2    
%     Phi_y(i) = mu_ty(i)*I_12;
% end
% DPhi_y = (Phi_y(2) - Phi_y(1)) / (mu_ty(2) - mu_ty(1));
% I_12_y = DPhi_y - mu_ty*I_11;

%%
% F_111_sum = 0;
% H_113_sum = 0;
% H_213_sum = 0;
% H_123_sum = 0;

%     H_113_sum = H_113_sum + atan2(a(i) * ...
%         (l2(i)*c1(i) - l1(i)*c2(i)),c1(i)*c2(i) + a(i)^2*l1(i)*l2(i));
    
%     F_111 = log((sqrt(l2(i)^2+g(i)^2)+l2(i))/(sqrt(l1(i)^2+g(i)^2)+l1(i)));

%     H_213_sum = H_213_sum + nu_xi(i) * F_111;
%     H_123_sum = H_123_sum + nu_eta(i) * F_111;

% H_111 = -abs(h) * H_113_sum + F_111_sum;
% H_113_step2 = (-H_111 + F_111_sum) / h^2

% H_113 = H_113_sum/abs(h);
% H_213 = -H_213_sum;
% H_123 = -H_123_sum;

%%
%     [int(i,:),a(i),l1(i),l2(i),c1(i),c2(i),g(i),nu_xi(i),nu_eta(i)] = ...
%         triEdge_fieldPoint_geom(verts_wrap(i,:),verts_wrap(i+1,:),P);
%     plot(int(i,1),int(i,2),'*k')
%     plot([x0 int(i,1)],[y0 int(i,2)],'k')
%     plot(x0+nu_xi(i)*a(i) , y0,'*m')
%     plot([x0 x0+nu_xi(i)*a(i)],[y0 y0],'m')
% %     plot([x0+nu_xi(i)*a(i) x0+nu_xi(i)*a(i)],[y0 y0+nu_eta(i)*a(i)],'m')
%     plot([x0+nu_xi(i)*a(i) int(i,1)],[y0 int(i,2)],'m')


%     F_111 = log((sqrt(l2(i)^2+g(i)^2)+l2(i))/(sqrt(l1(i)^2+g(i)^2)+l1(i)));
%     F_111_sum = F_111_sum + a(i)*F_111;
%     
%     H_113_sum = H_113_sum + atan2(a(i) * ...
%         (l2(i)*c1(i) - l1(i)*c2(i)),c1(i)*c2(i) + a(i)^2*l1(i)*l2(i));
