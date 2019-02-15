
% Jake Davis
% 10/06/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Testing Script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear
close all

addpath D:\Desktop\Thesis\Code\MATLAB\Generic;

%% Notes

% Subsonic edge
%    - The entire edge is always inside the Mach cone from the leading point of the edge. Thus
%      if the Mach cone intersects the edge upstream of the edge trailing point, then the 
%      intersection point is the new leading point of the edge, so the Mach cone emenates from
%      this new point

% Supersonic edge
%    - No matter what, if the DoD Mach cone intersects the edge, the respective R will be 0.
%      PRESUMABLY, the signs of the ym's don't change whether using the original edge point,
%      or the calculated intersection point. Thus the calculations are the same for both.

% Generic tri1: verts = [ 0 sind(60) 0; -.5 .1 0; .5 0 0];
%    - When subsonic edge is cut, intersection point is new LEADING pnt

% Generic tri2: verts = [ 0 sind(60) 0; -.5 .1 0; .5 0 0];
%    - When subsonic edge is cut, intersection point is new TRAILING pnt

%% Other

% % Worst Case Edges
% Mach = sqrt(2);
% % Mach = 2;
% verts = [0 -.51 0; -.5 0 0; .1 .1 0];   % worst case edges for M = sqrt(2)
% % verts = [.001 -.5001 0; -.5 0 0; -.001 0 0];   % worst case edges for M = sqrt(2)
% % verts = [.1 -.5001 0; -.5 0 0; -.001 0 0]; % sub parallel edge only
% % verts = [.001 -.5001 0; -.5 .2 0; -.001 0 0]; % sup perp edge only
% % P = [-.4 .15 0.1];
% % P = [-.2 .15 0.1];
% P = [0 .15 0.1];
% % P = [.5 -.2 0.05];
% nDir = 1;

% % Symmetric Tri
% Mach = sqrt(2);
% verts = [ -.5 .5 0; -.5 -.5 0; .5 0 0];
% P = [0.1 0 0.05];
% nDir = 1;

% % 'Subsonic' and 'Supersonic' Tri
% Mach = sqrt(2);
% verts = [-.25 0 0;  .25 -.25 0; .25 .25 0]; % Subsonic edge
% % verts = [-.25 0 0; -.05 .25 0; -.05 -.25 0];
% nDir = -1;
% P = [-.3 .1 .1];

% % Johnson Comparison
% Mach = sqrt(2);
% verts = [ 0 sind(60) 0; -.5 .1 0; .5 0 0];      % generic tri
% nDir = 1;
% 
% % P = [1.5 .3 .1];    % fully inside
% % P = [.3 1 .1];      % sup edges intersection
% % P = [0 -.1 .1];     % sub edges intersection
% P = [0 .6 .1];      % pure edge intersection

%% CPanel Output Tests

% format long
% verts = csvread('stuff.csv',0,0,[0 0 2 2]);
% P = csvread('stuff.csv',3,0);
% % P = [2 0 .5];
% format short
% % nDir = 1;
% nDir = -1;
% % Mach = sqrt(2);
% Mach = 1.745;

% x1 = verts(1,1); z1 = verts(1,3);
% x2 = verts(3,1); z2 = verts(3,3);
% m = (z2-z1) / (x2-x1)
% theta = atand(m)

%% Traveling Tests

% % verts = [ 0 sind(60) 0; -.5 .1 0; .5 0 0];      % generic tri1
% % verts = [ 0 sind(60) 0; -.3 .1 0; .5 0 0];      % generic tri1 - steeper sup edge
% verts = [ 0 sind(60) .1; -.5 .1 -.2; .5 0 .1];    % generic tri, off plane
% % verts = [ 0 sind(60) 0; -.5 -.1 0; .5 0 0];      % generic tri2
% nDir = 1;
% Mach = sqrt(2);
% % Mach = 2;
% % Mach = 3;
% 
% % P = [-.1 1 .2];
% % P = [-.3 .6 .2];    % wedge travel
% % P = [.2 .6 .2];    % wedge
% P = [.2 -.4 .2];    % wedge
% % P = [0 -.2 .2];
% % P = [-.5 .1 .2];

% P = [0 .1 .2];

% P = [0 .3 .2];
% P = [1 .1 .2];     % fully inside

% x-Travel
% % P = [0 .1 .2];            % bottom
% P = [-.2 .5 .1];            % center
% % P = [.6 .5 .1];            % center

% % P = [-.4 .2 .1];            % center
% % % P = [.8 .2 .1];            % center
% % P = [0 .5 .1];            % center
% P = [0 .5 -.2];            % center

% P = [0 .9 .2];            % top

% y-Travel
% P = [1 -.5 .05];
% P = [.75 -.5 .05];
% P = [.75 0 .05];

% z-Travel
% P = [.6 0.4 -.5];

%% THESIS - DOD Check

% verts = [ 0 sind(60) 0; -.5 .1 0; .5 0 0];      % generic tri1
% P = [.3 .65 .2];    % wedge
% Mach = 2;
% nDir = -1;
% 
% % Stuff
% [A,Q,cond] = triParams(verts,P,1);
% supPlotMachCone(-1,P(1),500,P,Mach,0)
% lgdPlot(1) = plot(nan, nan, strcat('-k'),'LineWidth',2);
% lgdStr{1} = sprintf('DOD intersection\nwith panel plane');
% lgdPlot(2) = plot(nan, nan, strcat('--k'),'LineWidth',2);
% lgdStr{2} = sprintf('DOD in plane\nof control point');
% lgdPlot(3) = plot(nan, nan, strcat('sr'),'MarkerSize',10,'MarkerFaceColor','r');
% lgdStr{3} = sprintf('Control point');
% lgdPlot(4) = plot(nan, nan, strcat('ok'),'MarkerSize',10,'MarkerFaceColor','k');
% lgdStr{4} = sprintf('Panel center');
% legend(lgdPlot, lgdStr,'Location','northeast')
% set(legend,'FontSize',16,'Interpreter','latex')
% legend boxoff
% set(gcf,'position',[0,0,1000,600])
% axis equal

%% THESIS - Transformation Example

% verts = [ 0 sind(60) 0; -.5 .1 0; .5 0 0];      % generic tri1
% P = [1.25 .3 .2];    % wedge
% Mach1 = sqrt(2);
% Mach2 = 3;
% nDir = 1;
% 
% % Combined - ref csys
% [A,Q,cond] = triParams(verts,P,1);
% hold on
% [xiVec_1,etaVec1_1,etaVec2_1,xiVecP_1,etaVecP1_1,etaVecP2_1] = ...
%     supPlotMachCone2(-1,P(1),500,P,Mach1,0);
% plot(xiVec_1,etaVec1_1,'b','LineWidth',2)
% plot(xiVec_1,etaVec2_1,'b','LineWidth',2)
% plot(xiVecP_1,etaVecP1_1,'--b','LineWidth',2)
% plot(xiVecP_1,etaVecP2_1,'--b','LineWidth',2)
% [xiVec_2,etaVec1_2,etaVec2_2,xiVecP_2,etaVecP1_2,etaVecP2_2] = ...
%     supPlotMachCone2(-1,P(1),500,P,Mach2,0);
% plot(xiVec_2,etaVec1_2,'r','LineWidth',2)
% plot(xiVec_2,etaVec2_2,'r','LineWidth',2)
% plot(xiVecP_2,etaVecP1_2,'--r','LineWidth',2)
% plot(xiVecP_2,etaVecP2_2,'--r','LineWidth',2)
% lgdPlot(1) = plot(nan, nan, strcat('-b'),'LineWidth',2);
% lgdStr{1} = '$M = \sqrt{2}$';
% lgdPlot(2) = plot(nan, nan, strcat('-r'),'LineWidth',2);
% lgdStr{2} = sprintf('$M = 3$');
% legend(lgdPlot, lgdStr,'Location','northeast')
% set(legend,'FontSize',16,'Interpreter','latex')
% xlim([-3 3]); ylim([-1.5 2.5])
% legend boxoff
% axis equal
% 
% % Other two plots done seperately later in code
% % Mach = sqrt(2);
% Mach = 3;

%% THESIS - Edge Mach cones

% Mach = sqrt(2);
% 
% % sonic edge plot
% verts1 = [0 -.5 0; -.5 0 0; 0 0 0];   % worst case edges for M = sqrt(2)
% supEdgePlotCones(verts1(2,:),verts1(1,:),.5,.5,Mach)
% axis equal
% 
% verts2 = [0 sind(60) 0; -.5 .1 0; .5 0 0];      % generic tri1
% % sup edge
% supEdgePlotCones(verts2(1,:),verts2(2,:),.5,.5,Mach)
% axis equal
% % sub edge
% supEdgePlotCones(verts2(2,:),verts2(3,:),0.5,1.5,Mach)
% lgdPlot(1) = plot(nan, nan, strcat('-ok'),'LineWidth',2,'MarkerSize',7,'MarkerFaceColor','w','MarkerEdgeColor','k');
% lgdStr{1} = 'Edge';
% lgdPlot(2) = plot(nan, nan, strcat('-b'),'LineWidth',2);
% lgdStr{2} = 'Pnt 1 Mach cone';
% lgdPlot(3) = plot(nan, nan, strcat('-r'),'LineWidth',2);
% lgdStr{3} = 'Pnt 2 Mach cone';
% legend(lgdPlot, lgdStr,'Location','northwest')
% set(legend,'FontSize',16,'Interpreter','latex')
% axis equal; xlim([-1.7 1.5])
% legend boxoff

%%

% [A,Q,cond] = triParams(verts,P,1);
% [xiVec_P,etaVec1_P,etaVec2_P,xiVecP_P,etaVecP1_P,etaVecP2_P] = ...
%     supPlotMachCone2(-1,P(1),500,P,Mach,0,1);
% plot(xiVec_P,etaVec1_P,'b','LineWidth',2)
% plot(xiVec_P,etaVec2_P,'b','LineWidth',2)
% plot(xiVecP_P,etaVecP1_P,'--b','LineWidth',2)
% plot(xiVecP_P,etaVecP2_P,'--b','LineWidth',2)
% axis equal
% grid on

%% z = 0 Testing

verts = [ 0 sind(60) 0; -.5 .1 0; .5 0 0];      % generic tri1
nDir = 1;
Mach = sqrt(2);
% P = [1 .4 .00001];    % wedge
% P = [1 .4 .000000001];    % wedge
P = [.7 .4 .000000001];    % wedge

[A,Q,cond] = triParams(verts,P,1);
supPlotMachCone(-.75,P(1),50,P,Mach,0)
axis equal

%% NEAR AND FAR FIELD COMPARISONS, POTENTIAL %%

% Define Stuff
sigma = 1;
% mu = [1 0 0];
mu = [1 1 1];
a = 0;
b = 0;
c0 = [cosd(a)*cosd(b) -sind(b) sind(a)*cosd(b)];

start = P(1);
% start = P(2);
% start = P(3);
varEnd = 5;
interv = .1;
var = start:interv:varEnd;

Pphi = P;
mu_0 = mu(1);
r = zeros(1,length(var));
cntr = Q;

phiDN = zeros(1, length(var));
phiSN = zeros(1, length(var));
phiDF = zeros(1, length(var));
phiSF = zeros(1, length(var));
for i = 1:length(var)
    Pphi(1) = var(i);
%     Pphi(2) = var(i);
%     Pphi(3) = var(i);

    r(i) = sqrt((Pphi(1)-Q(1))^2 + (Pphi(2)-Q(2))^2 + Pphi(3)^2);
    
    % POTENTIALS %
    
    % Near field
    [phiDN(i), phiSN(i)] = supLinearNearField(Mach,mu,sigma,verts,Pphi,cntr,a,b,c0,nDir);
    
    % Far Field
    [phiDF(i), phiSF(i)] = supLinearFarField(mu,sigma,verts,Pphi,Q,Mach,nDir);
end

%% PLOTTING, POTENTIALS

% Supersonic Doublet

figure()
hold on
% Near Field
plot(var,phiDN,'x-k','LineWidth',1)
% Far Field
plot(var,phiDF,'o-r','LineWidth',1)
legend('Linear Dub', 'Point Dub','Location','southeast')
title('SS Linear Doublet','FontSize',16)
xlabel('X Position','FontSize',13,'FontWeight','bold')
ylabel('$\Phi$','Interpreter','latex','FontSize',18,'FontWeight','bold')
set(get(gca,'ylabel'),'rotation',0)
grid on
hold off

% Supersonic Source

figure()
hold on
% Near Field
plot(var,phiSN,'x-k','LineWidth',1)
% Far Field
plot(var,phiSF,'o-r','LineWidth',1)
legend('Constant Src', 'Point Src')
title('SS Constant Source','FontSize',16)
xlabel('X Position','FontSize',13,'FontWeight','bold')
ylabel('$\Phi$','Interpreter','latex','FontSize',18,'FontWeight','bold')
set(get(gca,'ylabel'),'rotation',0)
grid on
hold off

%% CPanel 'Ellipse' Test Cases

% % CPanel Test Case 1 - Aft point, fwd-bot-star tri
% Mach = sqrt(2);
% verts = [-2 0 0; 0 2*tand(15) 0; 0 0 -2*tand(15)];
% P = [2 0 0];

% % CPanel Test Case 2 - Port pnt, fwd-bot-star tri
% Mach = sqrt(2);
% verts = [-2 0 0; 0 2*tand(15) 0; 0 0 -2*tand(15)];
% P = [0 -2*tand(15) 0];

% % CPanel Test Case 3 - Shared starboard point, fwd-bot-star tri
% Mach = sqrt(2);
% verts = [-2 0 0; 0 2*tand(15) 0; 0 0 -2*tand(15)];
% P = [0 2*tand(15)-0.000001 0];

% % CPanel Test Case 4 - Shared forward point, fwd-bot-star tri
% Mach = sqrt(2);
% verts = [-2 0 0; 0 2*tand(15) 0; 0 0 -2*tand(15)];
% P = [-2+0.000001 0 0];

% % CPanel Test Case 5 - Shared aft point, aft-bot-star tri
% Mach = sqrt(2);
% verts = [0 2*tand(15) 0; 2 0 0; 0 0 -2*tand(15)];
% P = [2-0.000001 0 0];

%% CPanel 'Diamond Test Cases

% % Stopped for intersection with floating edge
% Mach = sqrt(2);
% verts = [0 -2 0; 0 0 0; 1 0 0.2679];
% P = [1 0 -.26789];
% % P = [1 0 .4];
% nDir = -1;

% % Stopped for intersection with floating edge
% Mach = sqrt(2);
% verts = [0 -2 0; 0 0 0; 1 0 -0.2679];
% P = [1 0 .267898];
% % P = [1 0 -.45];
% nDir = 1;

% % Stopped at end of DOI check
% Mach = sqrt(2);
% verts = [1 0 -.2679; 2 2 0; 2 0 0];
% P = [1 0 .26789];

% % Stopped in panel intersection check
% Mach = sqrt(2);
% verts = [1 0 -.2679; 2 0 0; 2 -2 0];
% P = [1 0 -.267898];
% % P = [1.998 0 0];

% Mach = sqrt(2);
% verts = [2 0 0; 2 -.5 0; 1.6263 0 -.1001];
% P = [1.3737 0 -.16778];
% nDir = 1;

% Mach = sqrt(2);
% verts = [1.5 -.65715 -.13396; 1.0 -.5 -.2679; 1.3737 0 -.16778];
% P = [1.3737 0 -.16778];
% nDir = 1;

% Mach = sqrt(2);
% verts = [2 2 0; 1.6263 2 -.100133; 1.5 2 0];
% P = [1.6263 1.99999 -.100133];
% nDir = 1;

% Mach = sqrt(2);
% verts = [0 -.5 0; 0 0 0; .374111 0 -.03170555];
% P = [1.9999 1 0];
% nDir = 1;

% Mach = sqrt(2);
% verts = [2 0 0; 2 -.2575 0; 1.74769 0 -.0213828];
% P = [2.5 0 0];
% nDir = 1;

%%

%%% Symmetric Triangle %%%

% verts = [ -.5 .5 0; -.5 -.5 0; .5 0 0];         % symmetric tri
% verts = [ -.5 .5 -.1; -.5 -.5 -.1; .5 0 0];         % symmetric tri offplane
% Mach = sqrt(2);
% P = [0 .7 .1];          % top pnt intersection
% P = [0 -.7 .1];         % bot pnt intersection
% P = [-.45 .6 .05];
% P = [-.45 -.6 .05];

%%% Vertical Freestream Triangle %%%
% % Mach = 3;
% Mach = sqrt(2);
% verts = [-.5 0 -.5; 0 0 .3; .5 0 -.5];
% % verts = [-.5 -.1 -.5; 0 0 .3; .5 .5 -.5];
% P = [.2 .1 .1];

%%

% Geometry Definition

% Panel

% verts = [ 0 -sind(60) 0; -.5 .1 0; .5 0 0];     % generic tri flipped
% verts = [ .001 -.5001 0; -.5 0 0; -.001 0 0];   % worst case edges for M = sqrt(2)

% Field Point (w.r.t. generic tri, and Mach = sqrt(2))
% Mach = sqrt(2);

% P = [.8 .5 .1];             % 3rd vertex cut off
% P = [3 .5 .5];              % Fully in
% P = [.01 .9 0];             % Just outside
% P = [1.5 .25 0.0001];       % z = 0, fully inside
% P = [.1 .6 .001];           % z = 0, 2nd pnt in only
% P = [-.25 .3 .1];           % 2nd pnt in only
% P = [0 .3 .1];              % 2nd pnt in only
% P = [.75 .3 .1];              % 2nd pnt in only
% P = [-.2 .4 .1];            % edge1 only, POI over panel

% Field Point (w.r.t. generic tri, and Mach = 3)
% Mach = 3;

% P = [.4 .5 .05];            % pure edge
% P = [.2 .5 .05];            % pure edge, pnt above panel

% P = [0.37 .5 .05];            % pure edge, pnt above panel
% P = [.2 .5 -.05];            % pure edge, pnt below panel

%%

% plot(r,-phiDN,'x-k','LineWidth',1)
% plot(r,phiDF,'o-r','LineWidth',1)

% plot(r,phiSN,'x-k','LineWidth',1)
% plot(r,phiSF,'o-r','LineWidth',1)
