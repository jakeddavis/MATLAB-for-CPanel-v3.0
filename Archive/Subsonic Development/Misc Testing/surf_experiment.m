
% Equation of a surface experimenting

clc, clear, close all

a = 0;
b = 0;

[x,y] = meshgrid(0:.05:1,0:.05:1);

% x = 0:.05:1;
% y = 0:.05:1;

z = a*x.^2 + b*y.^2;

% for i = 1:length(x)
%     z(i) = a*x(i)^2 + b*y(i)^2;
% end

surf(x,y,z)