% 306

% Marshall's NACA function

%function [xu, yu, xl, yl] = naca4series(m,p,t,n)
function [xu, yu, xl, yl] = naca4series(m,p,t,xc,c)
    
   % generate actual m,p values
    mm = m/100;
    pp = p/10;
    
    % generate the initial x-locations
    %xc = linspace(0,1,n);
    
    % generate the camber line
    yc = zeros(1, length(xc));
    thetac = zeros(1, length(xc));
    for i = 1:length(xc)
        if (xc(i) < pp)
            yc(i) = (mm/pp^2)*(2*pp*xc(i) - xc(i)^2);
            thetac(i) = atan((2*mm/pp^2)*(pp - xc(i)));
        else
            yc(i) = (mm/(1-pp)^2)*((1 - 2*pp) + 2*pp*xc(i) - xc(i)^2);
            thetac(i) = atan((2*mm/(1-pp)^2)*(pp - xc(i)));
        end
    end
    
    % calculate the thickness
    yt = (t/20.0)*(0.29690*sqrt(xc)-0.12600*xc-0.35160*xc.^2+0.28430*xc.^3-0.10150*xc.^4);
    
    % calculate the upper coordinates
    xu = c*(xc - yt.*sin(thetac));
    yu = c*(yc + yt.*cos(thetac));
    
    % calculate the lower coordinates
    xl = c*(xc + yt.*sin(thetac));
    yl = c*(yc - yt.*cos(thetac));
end