
% Jake Davis
% 01/09/2019

clc, clear, close all

addpath D:\Desktop\Thesis\Code\MATLAB\Generic;

%%

mu = [1 0 0];
% mu = [1 1 1];
alpha = 0;      % angle of attack in degrees
beta  = 0;      % side slip angle in degrees
% verts = [ 0 sind(60) 0; -.5 .1 0; .5 0 0];

verts = csvread('stuff.csv',0,0,[0 0 2 2]);
P = csvread('stuff.csv',3,0);
% P = [2 0 .5];

% verts = [ .5 sind(60) 0; 0 .3 0; 1 0 0];
% P = [.9 .5 .05];

[A,Q,cond] = triParams(verts,P,1);
axis equal

local = getLocalSys(verts,0);
verts(1,:) = local * (verts(1,:)-Q)';
verts(2,:) = local * (verts(2,:)-Q)';
verts(3,:) = local * (verts(3,:)-Q)';
P = local * (P - Q)';

start = P(3);
varEnd = 7;
interv = .1;
var = start:interv:varEnd;

% verts = verts - Q;
% P = P - Q;
% Q = [0 0 0];

%% NEAR AND FAR FIELD COMPARISONS, POTENTIAL %%

% Preallocate
Pphi = P;
mu_0 = mu(1);
PhiNCD = zeros(length(var),1);
PhiFCD = PhiNCD;
PhiNLD = PhiNCD;
PhiFLD = PhiNCD;

for i = 1:length(var)
    
%     Pphi(1) = var(i);
%     Pphi(2) = var(i);
    Pphi(3) = var(i);
    
    % POTENTIALS %
    
    % Near Field, Constant Doublet
    PhiNCD(i) = subConstantDoubletNearField(mu_0,verts,Pphi);
    % Far Field, Constant Doublet
    PhiFCD(i) = subConstantDoubletFarField(mu_0,A,Q,Pphi);

    % Near Field, Linear Doublet
    PhiNLD(i) = subLinearDoubletNearField(mu,verts,Pphi,cond);
    % Far Field, Linear Doublet
    PhiFLD(i) = subLinearDoubletFarField(mu,verts,A,Q,Pphi);

end

%% PLOTTING, POTENTIALS

xPlot = (var - Q(3)) / (2*cond(1));

figure
hold on
plot(xPlot,PhiFCD,'o-k','LineWidth',1)
plot(xPlot,PhiFLD,'x-k','LineWidth',1)
plot(xPlot,PhiNCD,'o-r','LineWidth',1)
plot(xPlot,PhiNLD,'x-r','LineWidth',1)
lgd = legend('Constant Far','Linear Far','Constant Near','Linear Near','Location','southeast');
lgd.FontSize = 12;
title('Doublets','FontSize',16)
xlabel('Distance','FontSize',13,'FontWeight','bold')
ylabel('$\Phi$','Interpreter','latex','FontSize',18,'FontWeight','bold')
% ylim([-.06 0])
set(get(gca,'ylabel'),'rotation',0)
grid on
hold off
