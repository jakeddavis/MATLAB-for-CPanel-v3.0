
% Jake Davis
% 11/11/2018

clc, clear
close all

gam = 1.4;

%% My Test Case - alpha = 0 deg

% % My test case
% M1 = sqrt(2);
% alpha = 0;
% eps = 5;     % deg
% theta1 = eps;
% theta2 = 2 * eps;    % deg
% 
% % Compression turns (symmetric top & bottom)
% [ M2,P2_P1 ] = Oblique_Solver( M1,theta1,gam );
% Cp2 = computeCp(gam,M1,P2_P1,1.0);
% Cp4 = Cp2;
% 
% % Expansion turns (symmetric top & bottom)
% [ M3,P3_P2 ] = Expansion_Solver( M2,theta2,gam );
% Cp3 = computeCp(gam,M1,P3_P2,P2_P1);
% Cp5 = Cp3;
% 
% % Coefficients
% Cn = 0.5*(Cp4 - Cp2) + 0.5*(Cp5 - Cp3);
% 
% Ca2 = 0.5*tand(theta1)*Cp2;
% Ca3 = -0.5*tand(theta1)*Cp3;
% Ca4 = 0.5*tand(theta1)*Cp4;
% Ca5 = -0.5*tand(theta1)*Cp5;
% 
% Ca = Ca2 + Ca3 + Ca4 + Ca5;

%% MARCAP Test Case - alpha = 0 deg

% % MARCAP test case
% M1 = 1.745;
% alpha = 0;
% theta1 = 6.005;         % deg
% theta2 = 2 * theta1;    % deg
% 
% % Compression turns (symmetric top & bottom)
% [ M2,P2_P1 ] = Oblique_Solver( M1,theta1,gam );
% Cp2 = computeCp(gam,M1,P2_P1,1.0);
% Cp4 = Cp2;
% 
% % Expansion turns (symmetric top & bottom)
% [ M3,P3_P2 ] = Expansion_Solver( M2,theta2,gam );
% Cp3 = computeCp(gam,M1,P3_P2,P2_P1);
% Cp5 = Cp3;

% % Coefficients
% Cn = 0.5*(Cp4 - Cp2) + 0.5*(Cp5 - Cp3);
% 
% Ca2 = 0.5*tand(theta1)*Cp2;
% Ca3 = -0.5*tand(theta1)*Cp3;
% Ca4 = 0.5*tand(theta1)*Cp4;
% Ca5 = -0.5*tand(theta1)*Cp5;
% 
% Ca = Ca2 + Ca3 + Ca4 + Ca5;

%% MARCAP Test Case - alpha = 2 & 5 deg

% MARCAP test case
M1 = 1.745;
alpha = 2;
% alpha = 5;
eps = 6.005;         % deg
theta1 = eps - alpha;
theta2 = 2 * eps;    % deg
theta3 = eps + alpha;

% Compression turn top - forward top ramp
[ M2,P2_P1 ] = Oblique_Solver( M1,theta1,gam );
Cp2 = computeCp(gam,M1,P2_P1,1.0)

% Expansion turn top - aft top ramp
[ M3,P3_P2 ] = Expansion_Solver( M2,theta2,gam );
Cp3 = computeCp(gam,M1,P3_P2,P2_P1)

% Compression turn bottom - forward bottom ramp
[ M4,P4_P1 ] = Oblique_Solver( M1,theta3,gam );
Cp4 = computeCp(gam,M1,P4_P1,1.0)

% Expansion turn bottom - aft bottom ramp
[ M5,P5_P4 ] = Expansion_Solver( M4,theta2,gam );
Cp5 = computeCp(gam,M1,P5_P4,P4_P1)

% Coefficients
Cn = 0.5*(Cp4 - Cp2) + 0.5*(Cp5 - Cp3);

Ca2 = 0.5*tand(theta1)*Cp2;
Ca3 = -0.5*tand(theta1)*Cp3;
Ca4 = 0.5*tand(theta1)*Cp4;
Ca5 = -0.5*tand(theta1)*Cp5;
Ca = Ca2 + Ca3 + Ca4 + Ca5;

Cl = Cn*cosd(alpha) - Ca*sind(alpha)
Cd = Cn*sind(alpha) + Ca*cosd(alpha)

%% 405 Notes Case

% % Conditions
% M1 = 2.4;
% alpha = 15;
% eps = 5;     % deg
% theta1 = alpha - eps;
% theta2 = 2 * eps;    % deg
% theta3 = alpha;
% 
% % Expansion turn - ramp 2
% [ M2,P2_P1 ] = Expansion_Solver( M1,theta1,gam );
% Cp2 = computeCp(gam,M1,P2_P1,1.0);
% 
% % Expansion turn - ramp 3
% [ M3,P3_P2 ] = Expansion_Solver( M2,theta2,gam );
% Cp3 = computeCp(gam,M1,P3_P2,P2_P1);
% 
% % Compression turn - ramp 4
% [ M4,P4_P1 ] = Oblique_Solver( M1,theta3,gam );
% Cp4 = computeCp(gam,M1,P4_P1,1.0);
% 
% % Coefficients
% Cn = 0.5*(Cp4 - Cp2) + 0.5*(Cp4 - Cp3);
% Ca = 0.5*tand(eps)*(Cp2 - Cp3);
% 
% Cl = Cn*cosd(alpha) - Ca*sind(alpha);
% Cd = Cn*sind(alpha) + Ca*cosd(alpha);

%% THESIS - Entropy Analysis (alpha = 0)

% % Wedge w/ varying thickness and Mach number
% % Plot: (x,y1,y2) = (Mach,Cp,entropy)
% 
% % Inputs
% Mach = 1.5:0.25:3.5;
% theta = 5; % deg;
% nMach = length(Mach);
% 
% % Outputs
% M = zeros(1,nMach);
% Cp = zeros(1,nMach);
% ds = zeros(1,nMach);
% 
% for i = 1:nMach
%     % Compression turn (symmetric top & bottom)
%     [ M(i),Cp(i),ds(i) ] = wedgeSolver(Mach(i),theta);
% end
% 
% figure
% % hold on
% 
% % Cp vs. Mach
% line(Mach,Cp,'Color','b')
% xlabel('Mach')
% ylabel('Cp')
% ax1 = gca; % current axes position
% ax1Pos = ax1.Position;
% 
% % ds vs. Mach
% ax2xLabel = 'ds';
% ax2 = axes('Position',ax1Pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none');
% set(get(ax2,'XLabel'),'String','ds')
% ax2.YColor = 'w';
% line(ds,Mach,'Parent',ax2,'Color','w')
% 
% grid on

%%

% % Inputs
% Mach = 1.25:0.25:3.5;
% theta = 5:5:20; % deg
% nTheta = length(theta);
% nMach = length(Mach);
% 
% % Outputs
% Cp = zeros(nTheta,nMach);
% ds = zeros(nTheta,nMach);
% 
% for i = 1:nTheta
%     for j = 1:nMach
%         % Compression turn (symmetric top & bottom)
%         [ ~,Cp(i,j),ds(i,j) ] = wedgeSolver(Mach(j),theta(i));
%     end
% end

% figure
% % hold on
% 
% % Cp vs. Mach
% line(Mach,Cp(1,:),'Color','b')
% ax1 = gca; % current axes position
% ax1Pos = ax1.Position;
% 
% % ds vs. Mach
% ax2 = axes('Position',ax1Pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none');
% ax2.YColor = 'w';
% line(ds(1,:),Cp(1,:),'Parent',ax2,'Color','r')
% 
% grid on




%% entropy calc test

% % Inputs
% alpha = 0; % deg
% M1 = 2.4;
% theta = 10;
% 
% [ M2,Cp,ds ] = wedgeSolver( M1,theta )

%%

% Cp = Expansion_Solver( M1,theta1,gam );

% [Cp2, Cp3] = compExpansCalc(M1, v1, alpha, theta1, theta2);
