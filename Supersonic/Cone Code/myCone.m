
% Jake Davis
% 01/07/2019

clc, clear, close all

%% Development Case (MARCAP Basic)

% % Cone Dimensions
% origin = [0 0 0];
% length = 1.0; % ft
% theta = 10; % cone half angle in degrees
% % Cone Meshing
% phi = -15; % radial panel spacing in degrees
% lSpace = 5; % # regions along length
% 
% % no base on the cones (commented out)
% myConeFunc(origin,length,theta,phi,lSpace)

%% Mach Variation

% % Cone Dimensions
% origin = [0 0 0];
% length = 1.0; % ft
% theta = 10; % degrees
% 
% % Cone Meshing
% phi = -10; % radial panel spacing in degrees
% lSpace = 15; % # regions along length
% 
% myConeFunc(origin,length,theta,phi,lSpace)

%% Theta Variation

% % Cone Dimensions
% origin = [0 0 0];
% length = 1.0; % ft
% theta = 2.5:2.5:25; % degrees
% thetaNames = {'2_50','5_00','7_50','10_00','12_50','15_00','17_50','20_00','22_50','25_00',};
% 
% % Cone Meshing
% phi = -10; % radial panel spacing in degrees
% lSpace = 15; % # regions along length
% 
% [~,lengthTh] = size(thetaNames);
% for i = 1:lengthTh
%     filename = ['coneTh',thetaNames{i},'.tri'];
%     myConeFunc(origin,length,theta(i),phi,lSpace,filename)
% end

%% Az. Variation Test - New function format

% Cone Dimensions
origin = [0 0 0];
length = 1.0; % ft
theta = 10; % degrees

% Cone Meshing
phi = 30; % radial panel spacing in degrees
lSpace = 15; % # regions along length

myConeFunc(origin,length,theta,phi,lSpace,1)

% myConeFunc(origin,length,theta,phiSpace,lSpace,writeFile,varargin)

%% Mach Variation 1

% % Cone Dimensions
% origin = [0 0 0];
% length = 1.0; % ft
% theta = 10; % degrees
% 
% % Cone Meshing
% phi = -15; % radial panel spacing in degrees
% lSpace = 10; % # regions along length
% 
% myConeFunc(origin,length,theta,phi,lSpace)
