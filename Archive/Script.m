
% Jake Davis
% 02/03/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Testing Script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear, close all

% addpath D:\jake_davis\Thesis\Code\InputFiles\InputFile;

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

% Supersonic testing
% verts = [ 0 -sind(30) 0; -.5 0 0; 0 0 0];   % worst case edges
verts = [ 0 sind(60) 0; -.5 .1 0; .5 0 0];
% P = [1.5 .5 0];
P = [1.5 .5 1];
[A,Q,cond] = triParams(verts,P,1);
% hold on
% 
% Mach = sqrt(2);
% B = sqrt(Mach^2-1);
% xitest = linspace(-1,1.49,100);
% for i = 1:length(xitest)
%     etatest1(i) = P(2) + sqrt((xitest(i)-P(1))^2 - P(3)^2)/B;
%     etatest2(i) = P(2) - sqrt((xitest(i)-P(1))^2 - P(3)^2)/B;
% end
% plot(xitest,etatest1,'b','LineWidth',1)
% plot(xitest,etatest2,'b','LineWidth',1)

Mach = 3;
B = sqrt(Mach^2-1);
xitest = linspace(-1,1.49,100);
for i = 1:length(xitest)
    etatest1(i) = P(2) + sqrt((xitest(i)-P(1))^2 - P(3)^2)/B;
    etatest2(i) = P(2) - sqrt((xitest(i)-P(1))^2 - P(3)^2)/B;
end
plot(xitest,etatest1,'r','LineWidth',1)
plot(xitest,etatest2,'r','LineWidth',1)

title('Reference CSYS','FontSize',16)
xlabel('x')
ylabel('y')
xlim([-1 2])
ylim([-2 2])
% lgdPlot(1) = plot(nan, nan, strcat('-b'),'LineWidth',1);
% lgdStr{1} = 'M = $\sqrt{2}$';
% lgdPlot(2) = plot(nan, nan, strcat('-r'),'LineWidth',1);
% lgdStr{2} = 'M = 3.0';
lgdPlot(1) = plot(nan, nan, strcat('-r'),'LineWidth',1);
lgdStr{1} = 'M = 3.0';
legend(lgdPlot, lgdStr,'Location','southeast');
set(legend,'FontSize',14,'Interpreter','latex')

verts = [ 0 sind(60) 0; -.5 .1 0; .5 0 0];
% P = [.8512 .4 .1];
% P = [.8512 0.02075 .1];
% P = [.75 .2 0.5];
% P = [1.1 .8 0.5];
% P = [1.35 0.1 0.5];
% P = [0.05 .9 0.5];
% P = [2 .8 0.5];
% [A,Q,cond] = triParams(verts,P,1);

cntr = Q;
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
    
%     % Near Field, Constant Source
%     PhiNCS(i) = -subConstantSourceNearField(sigma,verts,Pphi);
%     % Far Field, Constant Source
%     PhiFCS(i) = subConstantSourceFarField(sigma,A,Q,Pphi);
    
%     % Near Field, Constant Doublet
%     PhiNCD(i) = subConstantDoubletNearField(mu_0,verts,Pphi);
%     % Far Field, Constant Doublet
%     PhiFCD(i) = subConstantDoubletFarField(mu_0,A,Q,Pphi);

    % Near Field, Linear Doublet
%     PhiNLD(i) = subLinearDoubletNearField(mu,verts,Pphi,cond,local);
%     PhiNLD(i) = CPsubLinearDoubletNearField(mu,verts,Pphi,cond,local);
%     % Far Field, Linear Doublet
%     PhiFLD(i) = subLinearDoubletFarField(mu,verts,A,Q,Pphi);
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
    
    [phiD_SS(i), phiS_SS(i)] = supLinearNearField(Mach,mu,sigma,verts,Pphi,cntr);
end

%% PLOTTING, POTENTIALS

% Supersonic Doublet
figure()
hold on
plot(r,phiD_SS,'x-k','LineWidth',1)
title('SS Linear Doublet','FontSize',16)
xlabel('Distance','FontSize',13,'FontWeight','bold')
ylabel('$\Phi$','Interpreter','latex','FontSize',18,'FontWeight','bold')
% ylim([-.06 0])
set(get(gca,'ylabel'),'rotation',0)
grid on
hold off
% 
% % Supersonic Source
% figure()
% hold on
% plot(r,phiS_SS,'x-k','LineWidth',1)
% title('SS Constant Source','FontSize',16)
% xlabel('Distance','FontSize',13,'FontWeight','bold')
% ylabel('$\Phi$','Interpreter','latex','FontSize',18,'FontWeight','bold')
% % ylim([-.06 0])
% set(get(gca,'ylabel'),'rotation',0)
% grid on
% hold off

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
