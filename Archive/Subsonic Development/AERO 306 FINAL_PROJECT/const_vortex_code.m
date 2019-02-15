
%306

function [Cl,Cm_c4] = const_vortex_code(Uinf,alpha,M,P,TT,c,npanels)

    %% Get Airfoil Points
    xc = 1 - cos(pi*linspace(0, 1, npanels/2 + 1)/2);
    %[xu, yu, xl, yl] = NACA_4dig(M,P,TT,xc,c);
    [xu, yu, xl, yl] = naca4series(M,P,TT,xc,c);
    xpts = cat(2, xl(end: - 1:2), xu);
    ypts = cat(2, yl(end: - 1:2), yu);

    xcolo = zeros(npanels, 1);
    ycolo = zeros(npanels, 1);
    nx = zeros(npanels, 1);
    ny = zeros(npanels, 1);
    tx = zeros(npanels, 1);
    ty = zeros(npanels, 1);
    len = zeros(npanels, 1);
    
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

    %% Find Strength of panels

    % pre-allocate stuff
    %gamma = zeros(npanels, 1);
    A = zeros(npanels + 1, npanels);
    b = zeros(npanels + 1, 1);

    % build coefficient matrix & rhs

    for i = 1:npanels % each collocation point
        for j = 1:npanels % panel vortex
            
            % coefficient matrix
            if (i == j)
                A(i,j) = -0.5;
            else
                [utemp, vtemp] = line_vortex_constant_2d(1, [xpts(j),xpts(j + 1)],...
                    [ypts(j), ypts(j + 1)], xcolo(i), ycolo(i));
                A(i,j) = tx(i)*utemp + ty(i)*vtemp;
            end
        end
  
        % right hand side
        b(i) = -(tx(i)*(Uinf*cos(alpha)) + ty(i)*(Uinf*sin(alpha)));
    end

    % enforce Kutta condition
    i = npanels + 1;
    A(i, 1) = 1;
    A(i, end) = 1;
    b(i) = 0;

    % calculate circulation at each panel
    gamma = A\b;

    %% Compute Properties
    
    GAMMA = dot(gamma,len);
    
    % calculate lift coefficient
    Cl = (2*GAMMA)/(Uinf*c);
    
    % calculate quarter chord moment coefficient
    Cl_panel = (2*gamma.*len)/(Uinf.*c); 
    Cl_norm = Cl_panel.*abs(ny);
    Cm_c4_v = zeros(npanels,1);
    for i = 1:length(xcolo)
        if xcolo(i) > 0.25*c
            Cm_c4_v(i) = -Cl_norm(i)*(abs((0.25*c) - xcolo(i)));
        else
            Cm_c4_v(i) = Cl_norm(i)*(abs((0.25*c) - xcolo(i)));
        end
    end
    
    Cm_c4 = sum(Cm_c4_v)/c;
    
    % velocity vectors


end

