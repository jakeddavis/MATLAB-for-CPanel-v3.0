
% Jake Davis
% 01/04/2019

clc, clear
close all

%% THESIS - Mach Variation (alpha = 0)

% % Wedge w/ varying Mach number
% % Plot: (x1,x2,y) = (Mach,entropy,Cp)
% 
% % Inputs
% mach = 1.5:0.25:3.5;
% theta = 5; % deg;
% nMach = length(mach);
% 
% % Outputs
% M = zeros(1,nMach);
% MaCp = zeros(1,nMach);
% MaDs = zeros(1,nMach);
% 
% for i = 1:nMach
%     % Compression turn (symmetric top & bottom)
%     [ M(i),MaCp(i),MaDs(i) ] = wedgeSolver(mach(i),theta);
% end
% 
% save('D:\Desktop\Thesis\Code\MATLAB\Thesis Plots\mach','mach')
% save('D:\Desktop\Thesis\Code\MATLAB\Thesis Plots\MaCp','MaCp')
% save('D:\Desktop\Thesis\Code\MATLAB\Thesis Plots\MaDs','MaDs')

%% THESIS - Theta Variation (alpha = 0)

% Wedge w/ varying wedge half angle
% Plot: (x1,x2,y) = (Mach,entropy,Cp)

% Use alpha as a proxy for half angle, only take top surface of wedge

% Inputs
mach = 2.0;
eps = 5; % deg;
alpha = 2.5:-2.5:-17.5;
theta = eps - alpha;
nTheta = length(theta);

% Outputs
M = zeros(1,nTheta);
ThCp = zeros(1,nTheta);
ThDs = zeros(1,nTheta); % J/kg K

for i = 1:nTheta
    % Compression turn (symmetric top & bottom)
    [ M(i),ThCp(i),ThDs(i) ] = wedgeSolver(mach,theta(i));
end

save('D:\Desktop\Thesis\Code\MATLAB\Thesis Plots\Theta Variation\theta','theta')
save('D:\Desktop\Thesis\Code\MATLAB\Thesis Plots\Theta Variation\ThCp','ThCp')
save('D:\Desktop\Thesis\Code\MATLAB\Thesis Plots\Theta Variation\ThDs','ThDs')

%%

% figure
% % hold on
% 
% % Cp vs. Mach
% line(Mach,Cp,'Color','b','LineStyle','none','Marker','o','MarkerSize',8)
% xlabel('Mach')
% ylabel('Cp')
% ax1 = gca; % current axes position
% ax1Pos = ax1.Position;
% grid on
% 
% % ds vs. Mach
% ax2xLabel = 'ds';
% ax2 = axes('Position',ax1Pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none');
% set(get(ax2,'XLabel'),'String','ds')
% ax2.YColor = 'w';
% line(ds,Mach,'Parent',ax2,'LineStyle','none')
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
