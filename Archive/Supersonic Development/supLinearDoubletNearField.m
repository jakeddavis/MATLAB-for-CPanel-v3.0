
% Jake Davis
% 04/07/2018

% 09/10/2018 - added source calc

function Phi = supLinearNearField(B,mu,verts,P,cntr)

% Domain of Dependence Check
depFlag = supDomainOfDependenceCheck(B,verts,P,cntr);

% Coord. transformations to panel local scaled
transMat = coordBody2supLocal(verts,B,0,0);
PSup = transMat * (P - cntr)';
vertsSup(1,:) = transMat * (verts(1,:) - cntr)';
vertsSup(2,:) = transMat * (verts(2,:) - cntr)';
vertsSup(3,:) = transMat * (verts(3,:) - cntr)';

% Redefine stuff
x = PSup(1); y = PSup(2); z = PSup(3);
vertsSup(4,:) = vertsSup(1,:);
% verts(4,:) = verts(1,:);

% Initialize stuff
Q1 = 0;
w0 = 0;
w0_m = 0;

for i = 1:3
    pnt1 = [vertsSup(i,1) vertsSup(i,2)];
    pnt2 = [vertsSup(i+1,1) vertsSup(i+1,2)];
%     pnt1 = [verts(i,1) verts(i,2)];
%     pnt2 = [verts(i+1,1) verts(i+1,2)];
%     PSup = P;
    
    triGeom = supTriGeom_func(pnt1,pnt2,PSup);
    xm  = triGeom.xm;    xmc  = triGeom.xmc;
    ym1  = triGeom.ym1;  ym2  = triGeom.ym2;
    ym1c = triGeom.ym1c; ym2c = triGeom.ym2c;
    R1   = triGeom.R1;   R2   = triGeom.R2;
    m    = triGeom.m;    lam  = triGeom.lam;
    
    if abs(m) < 1       % subsonic edge
        Q1temp = atan2(z*xmc * (ym1c*R2-ym2c*R1) , z^2*ym1c*ym2c+xmc^2*R1*R2);
        w0temp = (m/sqrt(1-m^2)) * log((ym2c+R2*sqrt(1-m^2))/(ym1c+R1*sqrt(1-m^2)));
    elseif abs(m) > 1   % supersonic edge
        Q1temp = atan2((z*xm*(ym1*R2-ym2*R1)) , (xm^2*R1*R2+z^2*ym2*ym1));
        w0temp = (1/sqrt(1-lam^2)) * atan2((sqrt(1-lam^2)*(ym1*R2-ym2*R1)) , ...
            (ym1*ym2+(1-lam^2)*R1*R2));
    end
    
    if abs(1 - abs(m)) < 1e-4   % sonic edge (Q1 is fine)
        zr = (ym1*R2 - ym2*R1) / (ym1*ym2 + (1-lam^2)*R1*R2);
        w0temp = zr * (1 - ((1-lam^2)*zr^2)/3 + ((1-lam^2)^2*zr^4)/5 ...
            - ((1-lam^2)^3*zr^6)/7); % + . . . (atan2 expansion)
    end
    
    Q1 = Q1 + Q1temp;
    w0 = w0 + w0temp;
    if abs(m) < 1e-4 || abs(m) > 1e4
%         w0_m = 0; % do nothing
    else
        w0_m = w0_m + w0/m;
    end
end

mu_0 = mu(1); mu_x = mu(2); mu_y = mu(3);
mu = mu_0 + mu_x*x + mu_y*y;

A = -Q1;
B = -z*w0;
C = -z*w0_m;  % fix

Phi = (1/(2*pi)) * mu*A - mu_x*B - mu_y*C;

end

%% 09/09/2018 - also probably wrong

    % new
%     
%     triGeom = supTriGeom_func2(pnt1Sup,pnt2Sup,PSup);
%     nu_xi = triGeom.nu_xi; nu_eta = triGeom.nu_eta;
%     l1 = triGeom.l1; l2 = triGeom.l2;
%     a = triGeom.a; b = triGeom.b;
%     R1 = triGeom.R1; R2 = triGeom.R2;
%     F1 = triGeom.F1; F2 = triGeom.F2;
%     
    % hH113
%     if h == 0
%         hH113_temp = 0;
%     elseif R1 == 0 && R2 == 0
%         hH113_temp = -pi * sign(h*nu_xi);
%     else
%         hH113_temp = -atan2((h*a*F1), (R1*R2 + h^2*F2));
%     end
%     
%     hH113 = hH113 + hH113_temp;
%     
%     % F111
%     if b > 0
%         if R1 == 0 && R2 == 0
%             F111_temp = pi/sqrt(b);
%         else
%             F111_temp = -(1/sqrt(b) * atan2((sqrt(b)*F1), F2));
%         end
%     elseif b < 0
%         if F2 > 5*sqrt(b)*abs(F1)
%             % stuff
%         end
%         F111_temp = -sign(nu_eta)/sqrt(b) * log((sqrt(b))*R1 + abs(l1))/((sqrt(b)*R2 + abs(l2)));
%     end
%     
%     F111 = F111 + -nu_xi*F111_temp;

%% 09/08/2018 - probably wrong

% hH113 = 0;
% F111 = 0;

%     % math from relationships
%     
% %     triGeom = triGeom_func(pnt1,pnt2,P);
%     triGeom = triGeom_func(pnt1Sup,pnt2Sup,PSup);
%     a  = triGeom.a;  g = triGeom.g;
%     l1 = triGeom.l1; l2 = triGeom.l2;
%     c1 = triGeom.c1; c2 = triGeom.c2;
%     nu_xi = triGeom.nu_xi; nu_eta = triGeom.nu_eta;
% 
%     if l1 >= 0 && l2 >= 0
%         F_111 = log((sqrt(l2^2+g^2)+l2) / (sqrt(l1^2+g^2)+l1));
%     elseif l1 < 0 && l2 < 0
%         F_111 = log((sqrt(l1^2+g^2)-l1) / (sqrt(l2^2+g^2)-l2));
%     elseif l2 >= 0 && l1 < 0
%         F_111 = log(((sqrt(l1^2+g^2)-l1) * (sqrt(l2^2+g^2)+l2)) / g^2);
%     end
% %     disp(F_111)
%     num_mine = a*(l2*c1-l1*c2);
%     denom_mine = c1*c2+a^2*l1*l2;
%     H_113_tan(i) = atan2(num_mine , denom_mine);
% %     H_213(i) = F_111 * nu_xi;
% %     H_123(i) = F_111 * nu_eta;
% 
%     w0_2 = w0_2 + -nu_xi*F_111;
%     w0_2Cm = w0_2/m;


%% 09/07/2018 - PANAIR/HOPM Domain Checking

% % depFlag = 0;    % domain of dependence flag
% Pc = P - P;     % set origin at point of interest (POI)
% vertsCheck = verts - P;     % adjust panel to new CSYS
% c0 = [1 0 0];
% 
% % copmpute stuff
% [~,Q,cond] = triParams(vertsCheck,Pc);
% X0 = Q(1); Y0 = Q(2);
% x0 = dot(Q-Pc,c0);
% y0 = sqrt(norm(Q-Pc)^2 - x0^2);
% if y0 >= B * x0
%     dist = sqrt((x0+B*y0)^2/(1+B^2));
% else
%     dist = sqrt(norm(Q - Pc)^2);
% end
% 
% % Check domain of dependence for panel
% % 0 = outside; 1 = inside; 2 = intersection->further assessment req'd
% % if (cntr(1)-)
% 
% if (dot(P - cntr, c0) >= 0) && (dot(P-cntr,P-cntr) >= 0)
%     if dist > cond(1)
%         if abs(Y0) < abs(X0/B)
%             disp('inside')
%             depFlag = 1;
%         else % no 'else' needed in actual code since using flags
%             disp('outside')
%         end
%     else
%         disp('intersection')
%         depFlag = 2;
%     end
% else % no 'else' needed in actual code since using flags
%     disp('outside')
% end
% 
% verts(4,:) = verts(1,:);
% 
% if depFlag == 2
%     for i = 1:3
%         pnt1 = verts(i,:);
%         pnt2 = verts(i+1,:);
% %         pnt1 = [verts(i,1) verts(i,2)];
% %         pnt2 = [verts(i+1,1) verts(i+1,2)];
%         
%         edgeUnit = (pnt2-pnt1) / norm(pnt2-pnt1);
%         if dot(edgeUnit, edgeUnit) > 0  % subsonic edge
%             disp('subsonic edge')
%         elseif dot(edgeUnit, edgeUnit) < 0  % supersonic edge
%             disp('supersonic edge')
%         else    % sonic edge
%             error('Uh oh 2...')
%         end
%     end
% end


%% Subsonic
% 
% x = P(1); y = P(2); h = P(3);
% verts(4,:) = verts(1,:);
% in_tri = cond(2);
% 
% H_113_tan = zeros(1,3);
% H_213 = zeros(1,3);
% H_123 = zeros(1,3);
% 
% for i = 1:3
%     pnt1 = [verts(i,1) verts(i,2)];
%     pnt2 = [verts(i+1,1) verts(i+1,2)];
%     
%     triGeom = triGeom_func(pnt1,pnt2,P);
%     a  = triGeom.a;  g = triGeom.g;
%     l1 = triGeom.l1; l2 = triGeom.l2;
%     c1 = triGeom.c1; c2 = triGeom.c2;
%     nu_xi = triGeom.nu_xi; nu_eta = triGeom.nu_eta;
% 
%     if l1 >= 0 && l2 >= 0
%         F_111 = log((sqrt(l2^2+g^2)+l2) / (sqrt(l1^2+g^2)+l1));
%     elseif l1 < 0 && l2 < 0
%         F_111 = log((sqrt(l1^2+g^2)-l1) / (sqrt(l2^2+g^2)-l2));
%     elseif l2 >= 0 && l1 < 0
%         F_111 = log(((sqrt(l1^2+g^2)-l1) * (sqrt(l2^2+g^2)+l2)) / g^2);
%     end
% 
%     H_113_tan(i) = atan2(a*(l2*c1-l1*c2) , c1*c2+a^2*l1*l2);
%     H_213(i) = F_111 * nu_xi;
%     H_123(i) = F_111 * nu_eta;
% end
% 
% % Check field point and compute H_113
% if (abs(h) == 0) && (in_tri == 0)
%     H_113 = 0;
% elseif (abs(h) == 0) && (in_tri == 1)
%     H_113 = (2*pi) / abs(h);
% else
%     H_111 = -abs(h)*sum(H_113_tan) + a*sum(F_111);
%     H_113 = (1/h^2) * (-H_111 + a*sum(F_111));
% end
% 
% H_213 = -sum(H_213);
% H_123 = -sum(H_123);
% 
% I_11 = h * H_113 / (4*pi);
% I_21 = h * H_213 / (4*pi);
% I_12 = h * H_123 / (4*pi);
% 
% mu_0 = mu(1); mu_x = mu(2); mu_y = mu(3);
% mu = mu_0 + mu_x*x + mu_y*y;
% 
% Phi = mu*I_11 + mu_x*I_21 + mu_y*I_12;
% 
% end

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
