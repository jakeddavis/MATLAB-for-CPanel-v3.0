% 306

% Jake Davis, Andrew Mercier, Jake Stone, AERO 306-01, Panel Code Project

% UNSTEADY FLOW VORTEX CODE

function [Cl,Cm] = unsteady_vortex_local_v(Uinf,alpha,M,P,TT,c,npanels,dt,t_end)

% OUTPUTS

% Cl ------------ vector of Cl's as time approaches t_end
% Cm ------------ vector of Cm_c/4 as time approaches t_end

% INPUTS

% Uinf ---------- freestream velocity [m/s]
% alpha --------- angle of attack [rad]
% M,P,TT -------- 4-dig NACA airfoil specification
% c ------------- airfoil chord length [m]
% npanels ------- number of panels
% dt ------------ time step
% t_end --------- time to be run to

% some vectors and matrices are preallocated by setting them equal to
% others that have already been defined by zeros() to decrease run time

%% GET AIRFOIL POINTS

xc = 1 - cos(pi*linspace(0, 1, npanels/2 + 1)/2);
[xu, yu, xl, yl] = naca4series(M,P,TT,xc,c);
xpts = cat(2, xl(end: - 1:2), xu);
ypts = cat(2, yl(end: - 1:2), yu);

% preallocate vectors
xcolo = zeros(npanels, 1);
ycolo = xcolo;
nx = xcolo;
ny = xcolo;
tx = xcolo;
ty = xcolo;
len = xcolo;

for i = 1:npanels
    % find collocation points
    xcolo(i) = (xpts(i + 1) + xpts(i))/2;
    ycolo(i) = (ypts(i + 1) + ypts(i))/2;
    
    % find unit tangent
    tx(i) = xpts(i + 1) - xpts(i);
    ty(i) = ypts(i + 1) - ypts(i);
    len(i) = sqrt(tx(i)^2 + ty(i)^2);
    tx(i) = tx(i)/len(i);
    ty(i) = ty(i)/len(i);
    
    % find unit normal
    nx(i) = -ty(i);
    ny(i) = tx(i);
end

%% BUILD COEFFICIENT MATRIX AND STEADY COMPONENT OF RIGHT HAND SIDE

% preallocate matrices
As = zeros(npanels + 1, npanels);   % steady state coefficient matrix
bs = zeros(npanels + 1, 1);         % steady state rhs
steady_b = zeros(npanels + 2, 1);      % steady state part of b for later use

% compute and build steady state components of A and b

for i = 1:npanels % each collocation point
    for j = 1:npanels % panel vortex
        
        % coefficient matrix
        if (i == j)
            As(i,j) = -0.5;
        else
            [utemp, vtemp] = line_vortex_constant_2d(1, [xpts(j),xpts(j + 1)],...
                [ypts(j), ypts(j + 1)], xcolo(i), ycolo(i));
            As(i,j) = tx(i)*utemp + ty(i)*vtemp;
        end
    end
    
    % right hand side
    bs(i) = -(tx(i)*(Uinf*cos(alpha)) + ty(i)*(Uinf*sin(alpha)));
end

steady_b(1:end-1) = bs;

% compute and build unsteady components of A

% preallocate matrices
A = zeros(npanels + 2, npanels + 1);    % final coefficeint matrix
b = zeros(npanels + 2, 1);              % final right hand side

% fill new A with values computed from steady state
for i = 1:npanels
    for j = 1:npanels
        A(i,j) = As(i,j);
    end
    b(i) = bs(i);
end

% enforce Kutta condition
i = npanels + 1;
A(i, 1) = 1;
A(i, end-1) = 1;

% fill the last row of A with the lengths of each panel
j = npanels + 2;
A(j,1:end-1) = len;
A(end,end) = 1;

% compute induced velocity on each panel caused by newest shed wake vortex

% average of the tangent vectors of the upper and lower trailing edges
t_te_x = (tx(end) - tx(1))/2;
t_te_y = (ty(end) - ty(1))/2;

% x and y locations of the starting/newly shed vortex
xstart = c + Uinf*dt*t_te_x;
ystart = Uinf*dt*t_te_y;

for i = 1:npanels
    [utemp, vtemp] = point_vortex(1,xstart,ystart,xcolo(i),ycolo(i));
    A(i,end) = tx(i)*utemp + ty(i)*vtemp;
end

%% PROPOGATE THROUGH TIME

% define time to be propogated through
time = (0:dt:t_end)';

% preallocate

% x and y locations of the shed wake vortices
x_loc_wv = double.empty(0);
y_loc_wv = double.empty(0);
% final x and y locations of shed wake vortices 
x_loc_wv_final = double.empty(0);
y_loc_wv_final = double.empty(0);

% induced velocities in the x and y caused by flow over the airfoil
u_induced_af = double.empty(0);
v_induced_af = double.empty(0);
% change in x and y of vortex locations caused by flow over the airfoil
xpush_af = double.empty(0);
ypush_af = double.empty(0);

% induced velocities in the x and y caused by wake vortices
u_induced_wv = double.empty(0);
v_induced_wv = double.empty(0);
% change in x and y of vortex locations caused by wake vortices
xpush_wv = double.empty(0);
ypush_wv = double.empty(0);

% change in wake vortex locations caused by freestream flow
% there is no change in the y-direction
xpush_Uinf = Uinf*dt;

b_wv = b;       % sums of every wake vortex's influence on each panel
wv = zeros(npanels,1);
gamma = bs;                         % circulation of each panel
GAMMA_tot = zeros(length(time),1);  % total circulation for each time step
GAMMA_wv = GAMMA_tot;               % strength of each shed vortex
Cl = GAMMA_tot;                     % Cl for each time step
Cm = GAMMA_tot;                     % Cm for each time step

% ITERATE THROUGH TIME

for t = 1:length(time)
    
    % COMPUTE LOCATIONS OF WAKE VORTICES

    if t == 1
        % at time t = 1, the location of the first shed vortex is given as
        x_loc_wv(1) = xstart;
        y_loc_wv(1) = ystart;
    else
        
        % compute wake vortex movement due to flow around the AIRFOIL
        
        for tt = 1:t-1
            % initialize 
            u_induced_af(tt) = 0;
            v_induced_af(tt) = 0;
            xpush_af(tt) = 0;
            ypush_af(tt) = 0;
            for i = 1:npanels
                [temp_u_induced_af,temp_v_induced_af] = line_vortex_constant_2d(gamma(i),...
                    [xpts(i),xpts(i+1)],[ypts(i),ypts(i+1)],x_loc_wv(tt),y_loc_wv(tt));
                u_induced_af(tt) = u_induced_af(tt) + temp_u_induced_af;
                v_induced_af(tt) = v_induced_af(tt) + temp_v_induced_af;
            end
            xpush_af(tt) = u_induced_af(tt)*dt;
            ypush_af(tt) = v_induced_af(tt)*dt;
        end
        
        % compute wake vortex movement due to OTHER VORTICES
        
        for tt = 1:t-1
            % initialize
            u_induced_wv(tt) = 0;
            v_induced_wv(tt) = 0;
            xpush_wv(tt) = 0;
            ypush_wv(tt) = 0;
            for w = 1:t-1
                % assume a vortex can't induce velocity on itself
                if tt == w
                    temp_u_induced_wv = 0;
                    temp_v_induced_wv = 0;
                else
                    [temp_u_induced_wv, temp_v_induced_wv] = point_vortex(GAMMA_wv(w),...
                        x_loc_wv_final(w),y_loc_wv_final(w),x_loc_wv_final(tt),y_loc_wv_final(tt));
                end
                u_induced_wv(tt) = u_induced_wv(tt) + temp_u_induced_wv;
                v_induced_wv(tt) = v_induced_wv(tt) + temp_v_induced_wv;
            end
            xpush_wv(tt) = u_induced_wv(tt)*dt;
            ypush_wv(tt) = v_induced_wv(tt)*dt;
        end
        xpush_wv_final = fliplr(xpush_wv)';
        ypush_wv_final = fliplr(ypush_wv)';
        
    % compute FINAL VORTEX LOCATIONS
        
        x_loc_wv(t) = 0;
        y_loc_wv(t) = 0;
        x_loc_wv_temp = x_loc_wv;
        y_loc_wv_temp = y_loc_wv;
        for tt = 2:t
            x_loc_wv(tt) = x_loc_wv_temp(tt-1) + xpush_Uinf + xpush_af(tt-1) + xpush_wv_final(tt-1);
            y_loc_wv(tt) = y_loc_wv_temp(tt-1) + ypush_af(tt-1) + ypush_wv_final(tt-1);
        end
    end    
             
    x_loc_wv_final = fliplr(x_loc_wv)';
    y_loc_wv_final = fliplr(y_loc_wv)';
    
    % ITERATE THROUGH PANELS
    
    for i = 1:npanels
        % fill b_wv with the sums of the induced velocity on each panel
        % caused by the all of the shed wake vortices
        for tt = 1:t
            [u_wv_temp,v_wv_temp] = point_vortex(GAMMA_wv(tt),x_loc_wv_final(tt),...
                y_loc_wv_final(tt),xcolo(i),ycolo(i));
            wv(tt) = tx(i)*u_wv_temp + ty(i)*v_wv_temp;
        end
        % sum of all shed wake vortices
        b_wv(i) = -sum(wv);
    end
    
    % build right hand side
    b = steady_b + b_wv;
    b(end) = -sum(GAMMA_wv);
    
    % calculate circulation at each panel and magnitude of new wv
    gamma = A\b;
    
    % get new wake vortex magnitudes
    GAMMA_wv(t) = gamma(end);
    gamma(end) = [];
    
    % COMPUTE PROPERTIES
    
    % calculate total circulation
    GAMMA_tot(t) = -sum(GAMMA_wv);
    
    % calculate lift coefficient
    Cl(t) = (2*GAMMA_tot(t))/(Uinf*c);
    
    % calculate quarter chord moment coefficient
    
    Cl_panel = (2*gamma.*len)/(Uinf*c); % Vector of Cl's of each panel
    Cl_norm = Cl_panel.*abs(ny);        % Component of Cl's in norm. dir.
    Cm_v = zeros(npanels,1);            % vector of 1/4 chord mom. coeff.
    for k = 1:length(xcolo)             % finds Cm's
        if xcolo(k) > 0.25*c
            % if the colocation point is to the right of the quarter chord,
            % then the sign of the moment about the quarter chord caused by
            % the lift of each panel needs to be switched
            Cm_v(k) = -Cl_norm(k)*(abs((0.25*c) - xcolo(k)));
        else
            Cm_v(k) = Cl_norm(k)*(abs((0.25*c) - xcolo(k)));
        end
    end
    Cm(t) = sum(Cm_v)/c;  % compute Cm for current time step
    
    % DISPLAY CURRENT TIME STEP
    disp(['time step = ',num2str(t)])
end

%% Plot Stuff

figure(1)
plot(x_loc_wv_final-c/c,fliplr(GAMMA_wv))
title('Vortex Strength vs d/c')
xlabel('d/c')
ylabel('Vortex Stength')

% airfoil and shed vortices
figure(2)
hold on
plot(xpts, ypts, '-b')
title('Airfoil with Vortex Strength Vizualiztion')
xlabel('x direction (m)')
ylabel('y direction (m)')
% wake vortices
for i = 1:length(x_loc_wv_final)
    plot(x_loc_wv_final(i), y_loc_wv_final(i),'mo','markersize',1500*abs(GAMMA_wv(i)))
end
axis equal

figure(3)
hold on
plot(xpts, ypts, '-b')
plot(x_loc_wv_final,y_loc_wv_final)
title('Linear Plot of Final Vortex Locations')
xlabel('x direction (m)')
ylabel('y direction (m)')
axis equal

figure(4)
hold on
plot(xpts, ypts, '-b')
plot(x_loc_wv_final,y_loc_wv_final,'r*')
title('Point Plot of Final Vortex Locations')
xlabel('x direction (m)')
ylabel('y direction (m)')
axis equal

end
