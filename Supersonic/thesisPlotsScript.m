
% Jake Davis
% 01/05/2019

%%%%%%%%%%%%%%%%%%%%%%%%% Thesis Script-Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear
close all

addpath D:\Desktop\Thesis\Code\MATLAB\Generic;

%% FAR FIELD COMPARISONS

% verts = [ 0 sind(60) 0; -.5 .1 0; .5 0 0];      % generic tri1
% nDir = 1;
% Mach = sqrt(2);
% P = [.8 sind(60) .15];
% 
% start = P(1);
% varEnd = 3;
% interv = .1;
% var = start:interv:varEnd;
% 
% [~,Q1,~] = triParams(verts,P,0);
% P(2) = Q1(2)-.1;
% 
% % Starting State
% hold on
% [~,~,~] = triParams(verts,P,1);
% supPlotMachCone(0,P(1),500,P,Mach,0);
% % Ending State
% P2 = P;
% P2(1) = var(end);
% plot3(P2(1),P2(2),P2(3),'sr','MarkerSize',10,'MarkerFaceColor','r')
% % [~,~,~] = triParams(verts,P2,1);
% supPlotMachCone(2,P2(1),500,P2,Mach,0);
% xlim([-.7 3.3])
% axis equal
% hold off

%% TRAVEL PLOTS

verts = [ 0 sind(60) 0; -.4 .1 0; .5 0 0];
nDir = 1;
Mach = sqrt(2);

% GENERIC PLOT

P = [.15 0.5 .15]; % starts outside then wegde
hold on
supPlotMachCone(-.6,P(1),500,P,Mach,0);
[~,Q,cond] = triParams(verts,P,1);
supPlot3Dcone(-.6,P(1),P,Mach)
% plot3(P(1),P(2),P(3),'sk','MarkerSize',8,'MarkerFaceColor','r')

% X-TRAVEL

% P = [-.35 0.5 .15]; % starts outside then wegde
% % P = [-.25 -0.1 .15]; % bottom
% % P = [.3 1 .15]; % top
start = P(1);
varEnd = 2;
interv = .05;

% hold on
% % supPlotMachCone(-.6,P(1),500,P,Mach,0);
% [~,Q,cond] = triParams(verts,P,1);
% % axis equal

% Y-TRAVEL

% P = [.6 -1 .15];
% start = P(2);
% varEnd = 1.5;
% interv = .05;
% % supPlotMachCone3(-.5,P(1),500,P,Mach,0,'top');
% hold on
% [~,Q,cond] = triParams(verts,P,1);
% % axis equal

% Z-TRAVEL

% P = [.8 .2 -1.3]; % starts in wegde
% start = P(3);
% varEnd = 1.3;
% interv = .05;
% hold on
% % supPlotMachCone(-.5,P(1),500,P,Mach,0);
% [~,Q,cond] = triParams(verts,P,1);
% axis equal

%% CALCULATE POTENTIAL %%

% Define Stuff
sigma = 1;
% mu = [1 0 0];
mu = [1 1 1];
a = 0;
b = 0;
c0 = [cosd(a)*cosd(b) -sind(b) sind(a)*cosd(b)];

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

%% FAR FIELD COMPARISONS

% xPlot = (var - cntr(1)) / (2*cond(1));
% travelPlots(xPlot,phiDN,phiSN,'Panel Diams.',phiDF,phiSF)

%% X-TRAVEL

% % data
% xPlot = (var - cntr(1)) / (2*cond(1));
% travelPlots2(xPlot,phiDN,phiSN,'Panel Diams.')

% % % state 1
% P1 = [-.025 P(2) P(3)];
% supPlotMachCone(-.45,P1(1),500,P1,Mach,0);
% plot3(P1(1),P1(2),P1(3),'sk','MarkerSize',10,'MarkerFaceColor','m')
% % state 2
% P2 = [0.15 P(2) P(3)];
% supPlotMachCone(-.32,P2(1),500,P2,Mach,0);
% plot3(P2(1),P2(2),P2(3),'sk','MarkerSize',10,'MarkerFaceColor','[0 0.7 0]')
% % state 3
% P3 = [0.365 P(2) P(3)];
% % P2 = [0.25 P(2) P(3)];
% supPlotMachCone(-.17,P3(1),500,P3,Mach,0);
% plot3(P3(1),P3(2),P3(3),'sk','MarkerSize',10,'MarkerFaceColor','b')
% % state 4
% P4 = [0.5 P(2) P(3)];
% supPlotMachCone(-.1,P4(1),500,P4,Mach,0);
% plot3(P4(1),P4(2),P4(3),'sk','MarkerSize',10,'MarkerFaceColor','y')
% state 5
% P5 = [1.1 P(2) P(3)];
% supPlotMachCone(.45,P5(1),500,P5,Mach,0);
% plot3(P5(1),P5(2),P5(3),'sk','MarkerSize',10,'MarkerFaceColor','c')
% xlim([-.4 1.2])
% axis equal
% x0=1;
% y0=3;
% width=9;
% height=7;
% set(gcf,'units','inches','position',[x0,y0,width,height])

% ax = gca;
% set(gca,'xticklabel',[],'yticklabel',[],'zticklabel',[])
% set(gca,'xtick',[],'ytick',[],'ztick',[])
% set(gca,'Color','w')
% 
% x0=1;
% y0=3;
% width=1.3;
% height=1.3;
% set(gcf,'units','inches','position',[x0,y0,width,height])
% 
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% axis equal
% axis off
% 
% % print(gcf,'figs\5_xPosState6.png','-dpng','-r600');


% 
% lgdPlot(1) = plot(nan, nan, strcat('sk'),'MarkerSize',10,'MarkerFaceColor','r');
% lgdStr{1} = sprintf('State 1');
% lgdPlot(2) = plot(nan, nan, strcat('sk'),'MarkerSize',10,'MarkerFaceColor','m');
% lgdStr{2} = sprintf('State 2');
% lgdPlot(3) = plot(nan, nan, strcat('sk'),'MarkerSize',10,'MarkerFaceColor','[0 0.7 0]');
% lgdStr{3} = sprintf('State 3');
% lgdPlot(4) = plot(nan, nan, strcat('sk'),'MarkerSize',10,'MarkerFaceColor','b');
% lgdStr{4} = sprintf('State 4');
% lgdPlot(5) = plot(nan, nan, strcat('sk'),'MarkerSize',10,'MarkerFaceColor','y');
% lgdStr{5} = sprintf('State 5');
% lgdPlot(6) = plot(nan, nan, strcat('sk'),'MarkerSize',10,'MarkerFaceColor','c');
% lgdStr{6} = sprintf('State 6');
% legend(lgdPlot, lgdStr,'Location','north')
% set(legend,'FontSize',14,'Interpreter','latex','Orientation','horizontal')
% legend boxoff
% hold off

% print(gcf,'figs\5_xPosGeom.png','-dpng','-r300');

%% Y-TRAVEL

% % data
% xPlot = (var - cntr(2)) / (2*cond(1));
% travelPlots2(xPlot,phiDN,phiSN,'Panel Diams.')
% axis equal

% state 1
% P1 = [P(1) -0.5 P(3)];
% supPlotMachCone3(-.45,P1(1),500,P1,Mach,0,'top');
% plot3(P1(1),P1(2),P1(3),'sk','MarkerSize',10,'MarkerFaceColor','m')
% state 2
% P2 = [P(1) 0.1 P(3)];
% supPlotMachCone3(-.2,P2(1),500,P2,Mach,0);
% plot3(P2(1),P2(2),P2(3),'sk','MarkerSize',10,'MarkerFaceColor','[0 0.7 0]')
% % state 3
% P3 = [P(1) 0.5 P(3)];
% supPlotMachCone3(0,P3(1),500,P3,Mach,0);
% plot3(P3(1),P3(2),P3(3),'sk','MarkerSize',10,'MarkerFaceColor','b')
% state 4
% P4 = [P(1) 1.2 P(3)];
% supPlotMachCone3(-.4,P4(1),500,P4,Mach,0,'bot');
% plot3(P4(1),P4(2),P4(3),'sk','MarkerSize',10,'MarkerFaceColor','y')
% ylim([-1.2 1.3])
% axis equal


% ax = gca;
% set(gca,'xticklabel',[],'yticklabel',[],'zticklabel',[])
% set(gca,'xtick',[],'ytick',[],'ztick',[])
% set(gca,'Color','w')
% 
% ylim([-.05 1.3])
% 
% x0=1;
% y0=3;
% width=1.3;
% height=1.3;
% set(gcf,'units','inches','position',[x0,y0,width,height])
% 
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% axis equal
% axis off

% print(gcf,'figs\5_yPosState5.png','-dpng','-r600');


% 
% lgdPlot(1) = plot(nan, nan, strcat('sk'),'MarkerSize',10,'MarkerFaceColor','r');
% lgdStr{1} = sprintf('State 1');
% lgdPlot(2) = plot(nan, nan, strcat('sk'),'MarkerSize',10,'MarkerFaceColor','m');
% lgdStr{2} = sprintf('State 2');
% lgdPlot(3) = plot(nan, nan, strcat('sk'),'MarkerSize',10,'MarkerFaceColor','[0 0.7 0]');
% lgdStr{3} = sprintf('State 3');
% lgdPlot(4) = plot(nan, nan, strcat('sk'),'MarkerSize',10,'MarkerFaceColor','b');
% lgdStr{4} = sprintf('State 4');
% lgdPlot(5) = plot(nan, nan, strcat('sk'),'MarkerSize',10,'MarkerFaceColor','y');
% lgdStr{5} = sprintf('State 5');
% legend(lgdPlot, lgdStr,'Location','west')
% set(legend,'FontSize',12,'Interpreter','latex')
% legend boxoff
% hold off

% print(gcf,'figs\5_yPosGeom.png','-dpng','-r300');

%% Z-TRAVEL

% % data
% xPlot = (var - cntr(3)) / (2*cond(1));
% travelPlots2(xPlot,phiDN,phiSN,'Panel Diams.')
% axis equal

% % state 1
% P1 = [P(1) P(2) -1.3];
% supPlotMachCone(-.6,P1(1),500,P1,Mach,0);
% supPlot3Dcone(-.6,P1(1),P1,Mach,1)
% plot3(P1(1),P1(2),P1(3),'sk','MarkerSize',8,'MarkerFaceColor','r')
% % state 2
% P2 = [P(1) P(2) -.9];
% supPlotMachCone(-.6,P2(1),500,P2,Mach,0);
% supPlot3Dcone(-.6,P2(1),P2,Mach,1)
% plot3(P2(1),P2(2),P2(3),'sk','MarkerSize',8,'MarkerFaceColor','m')
% % state 3
% P3 = [P(1) P(2) -.47];
% supPlotMachCone(-.4,P3(1),500,P3,Mach,0);
% supPlot3Dcone(-.4,P3(1),P3,Mach,1)
% plot3(P3(1),P3(2),P3(3),'sk','MarkerSize',8,'MarkerFaceColor','[0 0.7 0]')
% % state 4
% P4 = [P(1) P(2) -.33];
% supPlotMachCone(-.25,P4(1),500,P4,Mach,0);
% supPlot3Dcone(-.25,P4(1),P4,Mach,1)
% plot3(P4(1),P4(2),P4(3),'sk','MarkerSize',8,'MarkerFaceColor','b')
% % state 4
% P4 = [P(1) P(2) 0];
% supPlotMachCone(-.25,P4(1),500,P4,Mach,0);
% supPlot3Dcone(-.25,P4(1),P4,Mach,1)
% plot3(P4(1),P4(2),P4(3),'sk','MarkerSize',8,'MarkerFaceColor','y')
% 
% x0=1;
% y0=3;
% width=1.5;
% height=1.5;
% set(gcf,'units','inches','position',[x0,y0,width,height])
% 
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
% ax.View = [35 30];
ax.View = [30 30];
axis equal
axis off

% print(gcf,'figs\5_zPosState5.png','-dpng','-r600');
% print(gcf,'figs\5_unitTestGeneral.png','-dpng','-r600');

% 
% lgdPlot(1) = plot(nan, nan, strcat('k'),'LineWidth',2);
% lgdStr{1} = sprintf('State 1');
% lgdPlot(2) = plot(nan, nan, strcat('b'),'LineWidth',2);
% lgdStr{2} = sprintf('State 2');
% lgdPlot(3) = plot(nan, nan, strcat('r'),'LineWidth',2);
% lgdStr{3} = sprintf('State 3');
% lgdPlot(4) = plot(nan, nan,'Color','[0 0.35 0]','LineWidth',2);
% lgdStr{4} = sprintf('State 4');
% lgdPlot(5) = plot(nan, nan,strcat('--k'),'LineWidth',2);
% lgdStr{5} = sprintf('State 5');
% legend(lgdPlot, lgdStr,'Location','west')
% set(legend,'FontSize',12,'Interpreter','latex')
% legend boxoff
% hold off
% 
% % print(gcf,'figs\5_zPosGeom.png','-dpng','-r300');

%% PLOTTING, POTENTIALS

% Supersonic Doublet

% figure
% hold on
% % Near Field
% plot(var,phiDN,'x-k','LineWidth',1)
% % % Far Field
% % plot(var,phiDF,'o-r','LineWidth',1)
% % legend('Linear Dub', 'Point Dub','Location','southeast')
% title('SS Linear Doublet','FontSize',16)
% xlabel('X Position','FontSize',13,'FontWeight','bold')
% ylabel('$\Phi$','Interpreter','latex','FontSize',18,'FontWeight','bold')
% set(get(gca,'ylabel'),'rotation',0)
% grid on
% hold off
% 
% % Supersonic Source
% 
% figure
% hold on
% % Near Field
% plot(var,phiSN,'x-k','LineWidth',1)
% % % Far Field
% % plot(var,phiSF,'o-r','LineWidth',1)
% % legend('Constant Src', 'Point Src')
% title('SS Constant Source','FontSize',16)
% xlabel('X Position','FontSize',13,'FontWeight','bold')
% ylabel('$\Phi$','Interpreter','latex','FontSize',18,'FontWeight','bold')
% set(get(gca,'ylabel'),'rotation',0)
% grid on
% hold off
