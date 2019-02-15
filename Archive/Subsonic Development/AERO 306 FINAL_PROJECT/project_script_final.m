% 306

% Jake Davis, Andrew Mercier, Jake Stone, AERO 306-01, Panel Code Project

clc, clear all, close all

%% Inputs

Uinf = 10;              % freestream velocity   [m/s]
alpha = 0;              % angle of attack       [rad]
c = 2;                  % airfoil chord length  [m]
M = 4;                  % NACA specification
P = 4;                  % NACA specification
TT = 12;                % NACA specifiaction
npanels = 100;          % number of panels
K = 0.1;                % time varying constant
dt = K*(c/Uinf);        % time step
t_end = 2;              % end of time

%% RUN UNSTEADY WITH SHED VORTICES MOVING AT THIER LOCAL VELOCITIES

[Cl,Cm] = unsteady_vortex_local_v(Uinf,alpha,M,P,TT,c,npanels,dt,t_end);

figure (5)
t = 0:dt:t_end;
t_plot = t*Uinf/c;
plot(t_plot,Cl)
title('Cl vs t*Uinf/c')
xlabel('t*Uinf/c')
ylabel('Cl')

figure (6)
t = 0:dt:t_end;
t_plot = t*Uinf/c;
plot(t_plot,Cm)
title('Cm vs t*Uinf/c')
xlabel('t*Uinf/c')
ylabel('Cm')
