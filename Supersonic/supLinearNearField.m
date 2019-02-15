
% Jake Davis
% 10/12/2018

function [phiD, phiS] = supLinearNearField(Mach,mu,sigma,verts,P,cntr,a,b,c0,nDir)

% Domain of Dependence Check
[depFlag, ~] = supDoDCheck(Mach,verts,P,cntr,c0,nDir);

if depFlag == 0
    phiD = 0;
    phiS = 0;
else
    % Initialize
    newVerts = zeros(3,3);

    % Coordinate Transformation and Scaling
    [transMat, areaCorrect] = supCoordTrans(verts,Mach,a,b,nDir);
    PSup = transMat * (P - cntr)';
    newVerts(1,:) = transMat * (verts(1,:) - cntr)';
    newVerts(2,:) = transMat * (verts(2,:) - cntr)';
    newVerts(3,:) = transMat * (verts(3,:) - cntr)';
    
%     % PLOTS FOR THESIS
%     [A,Q,cond] = triParams(newVerts,PSup,1);
%     hold on
% %     supPlotMachCone(-1,PSup(1),500,PSup,Mach,0);
%     [xiVec_1,etaVec1_1,etaVec2_1,xiVecP_1,etaVecP1_1,etaVecP2_1] = ...
%     supPlotMachCone2(-1,PSup(1),500,PSup,Mach,1,1);
%     plot(xiVec_1,etaVec1_1,'b','LineWidth',2)
%     plot(xiVec_1,etaVec2_1,'b','LineWidth',2)
%     plot(xiVecP_1,etaVecP1_1,'--b','LineWidth',2)
%     plot(xiVecP_1,etaVecP2_1,'--b','LineWidth',2)
% % %     plot(xiVec_1,etaVec1_1,'r','LineWidth',2)
% % %     plot(xiVec_1,etaVec2_1,'r','LineWidth',2)
% % %     plot(xiVecP_1,etaVecP1_1,'--r','LineWidth',2)
% % %     plot(xiVecP_1,etaVecP2_1,'--r','LineWidth',2)
% % %     xlim([-3 3]); ylim([-1.5 2.5])
%     axis equal
%     grid on

    c0new = transMat * c0';

    % Define stuff
    x = PSup(1); y = PSup(2); z = PSup(3);
    delMag = 4;
    delSmall = 10^-delMag;
    delLarge = 10^delMag;
    
    % Initialize stuff
    Q1 = 0;
    w0 = 0;
    w0_m = 0;
    xm_w0 = 0;
    
    newVerts(end+1,:) = newVerts(1,:);
    for i = 1:3        
        
        pnt1 = newVerts(i,:);
        pnt2 = newVerts(i+1,:);

        % Original
        [Q1temp, w0temp, xm, m] = ...
            supInfCoeffs(pnt1,pnt2,PSup,c0new,delMag);
        
        % Doublet coeffs
        Q1 = Q1 + Q1temp;
        w0 = w0 + w0temp;
        if abs(m) < delSmall || abs(m) > delLarge
%             w0_m = w0_m; % do nothing
        else
            w0_m = w0_m + w0temp/m;
        end
        
        % Source Coeff
        xm_w0 = xm_w0 + xm*w0temp;
    end
    
    % Doublet calcs
    mu_0 = mu(1); mu_x = mu(2); mu_y = mu(3);
    mu = mu_0 + mu_x*x + mu_y*y;
    
    if abs(z) < delSmall
%     if abs(z) < 1e-10
        B = 0;
        C = 0;
        D = xm_w0;
    else
        B = -z*w0;
        C = -z*w0_m;
        D = xm_w0 - z*Q1;
    end
    A = -Q1;
    
    phiD = (1/(2*pi)) * (mu*A - mu_x*B - mu_y*C);
        
    % Source calcs
%     D = (xm_w0 - z*Q1) / areaCorrect;    
    phiS = (1/(2*pi)) * sigma * D;
end

% end

%%

%     % Panel-Mach cone intersection point %%%%%%%%%% not sure that this will
%     % work with any AoA or SS
%     interPnt = supMinConePanelInter(verts,P,cntr,Mach,c0,nDir);

%     interPntSup = transMat * (interPnt - cntr)';

%     plot3(interPntSup(1),interPntSup(2),interPntSup(3),'m*')

%         % Edge Intersection Testing
%         pnt1o = verts(i,:);
%         pnt2o = verts(i+1,:);
%         [interPnt] = supEdgeCheckTesting_1126(pnt1o,pnt2o,P,sqrt(Mach^2-1));
% %         pnt1 = (transMat * (interPnt - cntr)')';
%         pnt2 = (transMat * (interPnt - cntr)')';
% %         plot3(pnt1(1),pnt1(2),pnt1(3),'*r')
%         plot3(pnt2(1),pnt2(2),pnt2(3),'*r')

%         % Using edge intersections
%         [Q1temp, w0temp, xm, m] = ...
%             supInfCoeffsInters(pnt1,pnt2,PSup,interPntSup,c0new,delMag,sqrt(Mach^2-1),cntr);
        
%         % Johnson Based Scheme
%         [~, ~, xm, m] = ...
%             supInfCoeffs(pnt1,pnt2,PSup,interPntSup,c0new,delMag);
%         [Q1temp, w0temp] = supInfCoeffsJohn(pnt1,pnt2,PSup,interPntSup,delMag,sqrt(Mach^2-1));

%% transMat test

%     vertsOrig = zeros(3,3);
%     POrig = transMat \ PSup + cntr';
%     vertsOrig(1,:) = transMat \ newVerts(1,:)' + cntr';
%     vertsOrig(2,:) = transMat \ newVerts(2,:)' + cntr';
%     vertsOrig(3,:) = transMat \ newVerts(3,:)' + cntr';
%     interPntOrig = transMat \ interPntSup + cntr';
%     [~,~,~] = triParams(vertsOrig,POrig,1);
%     plot3(interPntOrig(1),interPntOrig(2),interPntOrig(3),'g*')

%% Plotting new trimmed panel

% [M,~] = size(newVerts);
% if M == 4
%     figure
%     hold on
%     plot(vertsSupTrim(1,1),vertsSupTrim(1,2),'*k')
%     plot(vertsSupTrim(2,1),vertsSupTrim(2,2),'*k')
%     plot(vertsSupTrim(3,1),vertsSupTrim(3,2),'*k')
%     plot(vertsSupTrim(4,1),vertsSupTrim(4,2),'*k')
%     title('Trimmed Local Scaled CSYS','FontSize',16)
%     grid on
% else
%     [~,~,~] = triParams(vertsSupTrim,PSup,1);
%     title('Trimmed Local Scaled CSYS','FontSize',16)
% end

%% Mach cone plotting before function made for it

% newP1 = supIntersectionPnt(PSup,vertsSup(1,:),vertsSup(2,:));
% newP2 = supIntersectionPnt(PSup,vertsSup(2,:),vertsSup(3,:));
% newP3 = supIntersectionPnt(PSup,vertsSup(3,:),vertsSup(1,:));
% xitest = linspace(-1,1.49,100);
% for i = 1:length(xitest)
%     etatest1(i) = PSup(2) + sqrt((xitest(i)-PSup(1))^2 - PSup(3)^2);
%     etatest2(i) = PSup(2) - sqrt((xitest(i)-PSup(1))^2 - PSup(3)^2);
% end
% xinew = PSup(1) - sqrt(PSup(3)^2);
% [A,Q,cond] = triParams(vertsSup,PSup,1);
% % title('Local Scaled CSYS, M = $\sqrt{2}$','FontSize',16,'Interpreter','latex')
% title('Local Scaled CSYS, M = 3','FontSize',16,'Interpreter','latex')
% hold on
% % plot(xitest,etatest1,'b','LineWidth',1)
% % plot(xitest,etatest2,'b','LineWidth',1)
% plot(xitest,etatest1,'r','LineWidth',1)
% plot(xitest,etatest2,'r','LineWidth',1)
% xlabel('x')
% ylabel('y')
% xlim([-1 2])
% ylim([-2 2])

% plot(newP1(1),newP1(2),'*m')
% plot(newP2(1),newP2(2),'*g')
% plot(newP3(1),newP3(2),'*y')
% plot(xinew(1),PSup(2),'*k')

%% Put Inf. Coeffs. calcs in one function

%     triGeom = supTriGeom_func(pnt1,pnt2,PSup,delMag);
%     xm   = triGeom.xm;   xmc  = triGeom.xmc;
%     ym1  = triGeom.ym1;  ym2  = triGeom.ym2;
%     ym1c = triGeom.ym1c; ym2c = triGeom.ym2c;
%     R1   = triGeom.R1;   R2   = triGeom.R2;
%     m    = triGeom.m;    lam  = triGeom.lam;
%     
%     if abs(m) < 1       % subsonic edge
%         Q1temp = atan2(z*xmc * (ym1c*R2-ym2c*R1) , z^2*ym1c*ym2c+xmc^2*R1*R2);
%         w0temp = (m/sqrt(1-m^2)) * log((ym2c+R2*sqrt(1-m^2))/(ym1c+R1*sqrt(1-m^2)));
%         xmSource = xmc;
%     elseif abs(m) > 1   % supersonic edge
%         Q1temp = atan2((z*xm*(ym1*R2-ym2*R1)) , (xm^2*R1*R2+z^2*ym2*ym1));
%         w0temp = (1/sqrt(1-lam^2)) * atan2((sqrt(1-lam^2)*(ym1*R2-ym2*R1)) , ...
%             (ym1*ym2+(1-lam^2)*R1*R2));
%         xmSource = xm;
%     end
%     
%     if abs(1 - abs(m)) < 10^-delMag   % sonic edge (Q1 is fine)
%         zr = (ym1*R2 - ym2*R1) / (ym1*ym2 + (1-lam^2)*R1*R2);
%         w0temp = zr * (1 - ((1-lam^2)*zr^2)/3 + ((1-lam^2)^2*zr^4)/5 ...
%             - ((1-lam^2)^3*zr^6)/7); % + . . . (atan2 expansion)
%     end

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
