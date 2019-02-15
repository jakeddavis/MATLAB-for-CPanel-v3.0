
% Jake Davis
% 09/28/2018

% Ellipse Exact Solution Script

clc, clear, close all

% freestream velocity
V = [20 0 0]';

% generate ellipsoid
a = 2;
b = 1;
c = 1;
[x,y,z] = ellipsoid(0,0,0,a,b,c,50);

% plot ellipsoid
figure
surf(x,y,z)
axis equal

[pot, u, v, w] = ellipsoidflow(x, y, z, V, a, b, c);

Vmag = sqrt(u^2 + v^2 + w^2);

figure
surf(x,y,z,pot)
% surf(x,y,z,Vmag)
colorbar
axis equal

max(max(pot))
min(min(pot))

