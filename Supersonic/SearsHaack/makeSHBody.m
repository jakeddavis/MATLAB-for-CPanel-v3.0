clc, clear, close all

% clear all
format compact

% Rs = [0.02 0.05 0.1 0.15 0.2 0.25 0.3];
% Rs = 0.02;
Rs = 0.037879;

L = 1.0;
% n = 51;
% m = 30;
n = 151;
m = 40;

for k=1:length(Rs);


Rmax = Rs(k); %0.25;
% L = 1.0;
% n=201;
% m=101;

x = linspace(0,1,n);

r = Rmax * (4*x.*(1-x)).^.75;

figure(1);
plot(x*L,r)
axis equal


thetas = fliplr(linspace(0,2*pi,m));

ycirc = cos(thetas);
zcirc = sin(thetas);

figure(2)
plot(ycirc,zcirc);
axis equal

pts = nan(n,m,3);

for i=1:n
  for j=1:m
    pts(i,j,1) = L*x(i);
    pts(i,j,2) = r(i)*ycirc(j);
    pts(i,j,3) = r(i)*zcirc(j);
  end
end


c.name = 'SHBody';
c.group = 0;
c.type = 1;
c.ncross = n;
c.npt = m;
c.pts = pts;

fname = ['SHBody_', num2str(Rmax), '.hrm'];

hrmwrite( c, fname);

figure(3)
hrmplot(c)

Dwave = 9*pi^3*Rmax^4/(2*L^2)


end


CD_C3D = [6.9332786e-06
  0.00027928616
  0.0040614218
  0.017767534
  0.046465767
  0.09180539
  0.15140754];

Rs = [0.02 0.05 0.1 0.15 0.2 0.25 0.3];

for k=1:length(Rs);
  Rmax = Rs(k); %0.25;
  L = 1.0;
  Dwave(k) = 9*pi^3*Rmax^4/(2*L^2);
end