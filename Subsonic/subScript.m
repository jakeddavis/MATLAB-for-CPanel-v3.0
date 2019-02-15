
% Jake Davis
% 10/06/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Testing Script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear, close all

addpath D:\Desktop\Thesis\Code\MATLAB\Generic;

alpha = 0;      % angle of attack in degrees
beta  = 0;      % side slip angle in degrees
Mach = sqrt(2);
% Mach = 0;
B = sqrt(Mach^2-1);

% % 8 tri ellipse test case
% verts = [0 0 0; 2 1 0; 2 0 -1];
% P = [4 0 0];
% 
% local = getLocalSys(verts);
% [A,Q,cond] = triParams(verts,P,1);
% hold on
% normal = getNorm(verts);

% % Supersonic testing
% % verts = [ 0 -sind(30) 0; -.5 0 0; 0 0 0];   % worst case edges
% verts = [ 0 sind(60) 0; -.5 .1 0; .5 0 0];
% % P = [1.5 .5 0];
% % P = [1.5 .5 1];
% P = [3 .5 .5];


verts = csvread('D:\Desktop\Thesis\Code\MATLAB\Supersonic\stuff.csv',0,0,[0 0 2 2]);
P = csvread('D:\Desktop\Thesis\Code\MATLAB\Supersonic\stuff.csv',3,0);


[A,Q,cond] = triParams(verts,P,1);

vertsOrig = verts;
POrig = P;
local = getLocalSys(verts,0);
verts(1,:) = local * (verts(1,:)-Q)';
verts(2,:) = local * (verts(2,:)-Q)';
verts(3,:) = local * (verts(3,:)-Q)';
P = local * (P - Q)';
% cntr = [0 0 0];

% Q = cntr;

% P = [2 0 .5];



% cntr = Q;
% cntr = [Q 0];
% start = 0.1;
start = P(1);
% varEnd = cond(1);
% varEnd = 0.7;
varEnd = 10;
interv = .1;

% verts = verts - Q;
% P = P - cntr;
% Q = [0 0];

% Define Stuff
sigma = 1;
% mu = [1 5 2];
mu = [1 0 0];



%% NEAR AND FAR FIELD COMPARISONS, POTENTIAL %%

% Preallocate
Pphi = P;
mu_0 = mu(1);
% v = -0.4:-0.1:-5*cond(1);
% v = 0.1:0.1:5*cond(1);
var = start:interv:varEnd;
r = zeros(1,length(var));
PhiNCS = zeros(length(var),1);
PhiFCS = PhiNCS;
PhiNCD = PhiNCS;
PhiFCD = PhiNCS;
PhiNLD = PhiNCS;
PhiFLD = PhiNCS;
% supPhi_near_lin_doublet = Phi_near_const_source;

uNLD = zeros(length(var),1); vNLD = uNLD; wNLD = uNLD;
uFCD = uNLD; vFCD = uNLD; wFCD = uNLD;
uNCD = uNLD; vNCD = uNLD; wNCD = uNLD;

for i = 1:length(var)
    Pphi(1) = var(i);
%     Pphi(2) = var(i);
%     Pphi(3) = var(i);
    r(i) = sqrt((Pphi(1)-Q(1))^2 + (Pphi(2)-Q(2))^2 + Pphi(3)^2);
    
    % POTENTIALS %
    
    % Near Field, Constant Source
    PhiNCS(i) = subConstantSourceNearField(sigma,verts,Pphi);
%     Phi = sourceNearField(sigma,verts,Pphi);
    Phi2(i) = CPanelStuff_0520(vertsOrig,POrig,mu,Q);
%     Phi2(i) = CPanelStuff_0520(verts,Pphi',mu,Q);
    % Far Field, Constant Source
    PhiFCS(i) = subConstantSourceFarField(sigma,A,Q,Pphi);
    
    % Near Field, Constant Doublet
    PhiNCD(i) = subConstantDoubletNearField(mu_0,verts,Pphi);
    % Far Field, Constant Doublet
    PhiFCD(i) = subConstantDoubletFarField(mu_0,A,Q,Pphi);

    % Near Field, Linear Doublet
    PhiNLD(i) = subLinearDoubletNearField(mu,verts,Pphi,cond);
%     PhiNLD(i) = CPsubLinearDoubletNearField(mu,verts,Pphi,cond,local);
    % Far Field, Linear Doublet
    PhiFLD(i) = subLinearDoubletFarField(mu,verts,A,Q,Pphi);
%     
%     % VELOCITIES %
%     
%     % Near Field, Constant Doublet - BROKEN
% %     [uNCD(i),vNCD(i),wNCD(i)] = subConstantDoubletNearFieldVel(mu_0,verts,Pphi);
%     % Far Field, Constant Doublet
%     [uFCD(i),vFCD(i),wFCD(i)] = subConstantDoubletFarFieldVel(mu_0,A,Q,Pphi);
%     
%     % Near Field, Linear Doublet
%     [uNLD(i),vNLD(i),wNLD(i)] = subLinearDoubletNearFieldVel(mu,verts,Pphi);
% 
%     % CPanel validation stuff
%     phi = CPanelStuff_0520(verts,Pphi,mu,Q);
    
end

%% PLOTTING, POTENTIALS

% Subsonic Doublet
figure()
hold on
plot(r,PhiFCD,'o-k','LineWidth',1)
plot(r,PhiFLD,'x-k','LineWidth',1)
plot(r,PhiNCD,'o-r','LineWidth',1)
plot(r,PhiNLD,'x-r','LineWidth',1)
lgd = legend('Constant Far','Linear Far','Constant Near','Linear Near','Location','southeast');
lgd.FontSize = 12;
title('Doublets','FontSize',16)
xlabel('Distance','FontSize',13,'FontWeight','bold')
ylabel('$\Phi$','Interpreter','latex','FontSize',18,'FontWeight','bold')
% ylim([-.06 0])
set(get(gca,'ylabel'),'rotation',0)
grid on
hold off

% Source
figure()
hold on
plot(r,PhiNCS,'r')
plot(r,PhiFCS,'b')
plot(r,Phi2,'k')
legend('Near','Far')
title('Sources'); xlabel('r'); ylabel('Phi')
hold off
grid on

%% NEAR AND FAR FIELD COMPARISONS, SUBSONIC VELOCITY %%

Pvel_x = P;
Pvel_y = P;
Pvel_z = P;
eps = 0.001;
velVar = start:eps:5*cond(1);
rFD = zeros(1,length(velVar));
% uNLD_FD = zeros(length(velVar),1); vNLD_FD = uNLD_FD; wNLD_FD = uNLD_FD;
% uFCD_FD = uNLD_FD; vFCD_FD = uNLD_FD; wFCD_FD = uNLD_FD;

% Compute u
uNLD_FD = zeros(length(velVar),1);
for i = 1:length(velVar)
    Pvel_x(1) = velVar(i);
    rFD(i) = sqrt((Pvel_x(1)-Q(1))^2 + (Pvel_x(2)-Q(2))^2 + Pvel_x(3)^2);
%     PhiNCD(i) = subConstantDoubletNearField(mu_0,verts,Pvel_x);
    PhiNLD(i) = subLinearDoubletNearField(mu,verts,Pvel_x,cond);
    if i == 1
        % do nothing
    else
        uNLD_FD(i-1) = (PhiNLD(i) - PhiNLD(i-1)) / eps;
    end
end

% Compute v
vNLD_FD = uNLD_FD;
for i = 1:length(velVar)
    Pvel_y(2) = velVar(i);
    rFD(i) = sqrt((Pvel_y(1)-Q(1))^2 + (Pvel_y(2)-Q(2))^2 + Pvel_y(3)^2);
%     PhiNCD(i) = subConstantDoubletNearField(mu_0,verts,Pvel_y);
    PhiNLD(i) = subLinearDoubletNearField(mu,verts,Pvel_y,cond);
    if i == 1
        % do nothing
    else
        vNLD_FD(i-1) = (PhiNLD(i) - PhiNLD(i-1)) / eps;
    end
end

% Compute w
wNLD_FD = uNLD_FD;
for i = 1:length(velVar)
    Pvel_z(3) = velVar(i);
    rFD(i) = sqrt((Pvel_z(1)-Q(1))^2 + (Pvel_z(2)-Q(2))^2 + Pvel_z(3)^2);
%     PhiNCD(i) = subConstantDoubletNearField(mu_0,verts,Pvel_z);
    PhiNLD(i) = subLinearDoubletNearField(mu,verts,Pvel_z,cond);
    if i == 1
        % do nothing
    else
        wNLD_FD(i-1) = (PhiNLD(i) - PhiNLD(i-1)) / eps;
    end
end

%% PLOTTING, VELOCITIES

% % u velocity
% figure
% hold on
% plot(rFD,uNLD_FD,'-r','LineWidth',1)
% plot(r,uFCD,'o-k','LineWidth',1)
% plot(r,-uNLD,'x-r','LineWidth',1)
% lgd = legend('uNLD FD','uFCD','uNLD','Location','northeast');
% lgd.FontSize = 12;
% title('u Velocity','FontSize',16)
% xlabel('Distance','FontSize',13,'FontWeight','bold')
% ylabel('u','FontSize',18,'FontWeight','bold')
% set(get(gca,'ylabel'),'rotation',0)
% grid on
% hold off

% % v velocity
% figure
% hold on
% plot(rFD,vNLD_FD,'-r','LineWidth',1)
% plot(r,vFCD,'o-k','LineWidth',1)
% plot(r,-vNLD,'x-r','LineWidth',1)
% lgd = legend('vNCD FD','vFCD','vNLD','Location','northeast');
% lgd.FontSize = 12;
% title('v Velocity','FontSize',16)
% xlabel('Distance','FontSize',13,'FontWeight','bold')
% ylabel('v','FontSize',18,'FontWeight','bold')
% set(get(gca,'ylabel'),'rotation',0)
% grid on
% hold off

% w velocity
figure
hold on
plot(rFD,wNLD_FD,'-r','LineWidth',1)
plot(r,wFCD,'o-k','LineWidth',1)
plot(r,-wNLD,'x-r','LineWidth',1)
lgd = legend('wNLD FD','wFCD','wNLD','Location','northeast');
lgd.FontSize = 12;
title('w Velocity','FontSize',16)
xlabel('Distance','FontSize',13,'FontWeight','bold')
ylabel('w','FontSize',18,'FontWeight','bold')
set(get(gca,'ylabel'),'rotation',0)
grid on
hold off

%%

% delta_h = 0.01;
% cond(4) = delta_h;
% cond(6) = cond(4);
