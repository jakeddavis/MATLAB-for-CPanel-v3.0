
% Jake Davis
% 01/09/2019

clc, clear, close all

% addpath('D:\Desktop\Thesis\Code\MATLAB\Mesh Generation\Dist Mesh\distmesh_2d')

%% Input Data

% Geometry
a = 2;
b = 1;
c = 1;

h0 = 0.4; % 1100 panels

% Flow
Vinf  = 20;  % ft/s
alpha = 15; % degrees
beta  = 10; % degrees
% alpha = 0; % degrees
% beta  = 0; % degrees
windDir = [cosd(alpha)*cosd(beta) -sind(beta) sind(alpha)*cosd(beta)];
V = Vinf * windDir;

%% Create Ellipse

[x,y,z] = ellipsoid(0,0,0,a,b,c,50);

% % .tri file name
% pre = '500';
% geom = 'ellipse211.tri';
% 
% % read tri data
% fid = fopen([pre,geom],'r');
% nData = sscanf(fgetl(fid),'%i');
% % nData = fscanf(fid,'%i %i');
% nNodes = nData(1);
% nPanels = nData(2);
% nodes = fscanf(fid,'%f %f %f',[3,nNodes])';
% conn = fscanf(fid,'%i %i %i',[3,nPanels])' + 1;
% fclose(fid);
% 
% x = nodes(:,1);
% y = nodes(:,2);
% z = nodes(:,3);

%% Solve Potential Flow

[pot, u, v, w] = ellipsoidflow(x, y, z, V, a, b, c);

Vmag = sqrt(u.^2 + v.^2 + w.^2);
Cp = 1 - (Vmag/Vinf).^2;

figure
% surf(x,y,z,Cp)
surf(x,y,z,pot)
colorbar
axis equal

% % figure
% trisurf(conn,x,y,z,Cp)
% shading interp
% axis equal
% colorbar

max(max(pot))
min(min(pot))

max(max(Cp))
min(min(Cp))

