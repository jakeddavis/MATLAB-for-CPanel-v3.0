
% Jake Davis
% 09/09/2018

function [depFlag, pntInterData] = supDomainOfDependenceCheck(B,verts,P,PSup,cntr)

% Define Stuff
depFlag = 0;    % initialize domain of dependence flag
% 0 = completely outside -> no influence
% 1 = completely inside
% 2 = intersection with at least one point inside
% 3 = intersection with no points inside, or outside
x = P(1); y = P(2); z = P(3);
x0 = cntr(1); y0 = cntr(2); z0 = cntr(3);
RBcntr = (x-x0)^2 - B^2*(y-y0)^2 - B^2*(z-z0)^2;
c0 = [1 0 0];

% compute stuff
Pc = P - P;     % set origin at point of interest (POI)
vertsCheck = verts - P;     % adjust panel to new CSYS

[~,Q,cond] = triParams(vertsCheck,Pc,0);
x0 = dot(Q-Pc,c0);
y0 = sqrt(norm(Q-Pc)^2 - x0^2);
if y0 >= B * x0
    dist = sqrt((x0+B*y0)^2/(1+B^2));
else
    dist = sqrt(norm(Q - Pc)^2);
end

PHyp = [PSup(1) - sqrt(PSup(3)^2), PSup(2), 0];
pntInterData = zeros(3,1);
if dot(P - cntr, c0) >= 0   % panel cntr is upstream of CP
    if dist > cond(1)   % cntr of panel is further from domain than panel rad.
        if sign(RBcntr) == 1     % cntr is inside domain
            depFlag = 1;
            disp('inside')
        elseif sign(RBcntr) == -1   % don't actually need
            disp('outside')      % cntr is outside domain
        else
            error('uh oh...')
        end
    else
        count = 0;
        for i = 1:3
            xi = verts(i,1);
            yi = verts(i,2);
            zi = verts(i,3);
            RBverts = (x-xi)^2 - B^2*(y-yi)^2 - B^2*(z-zi)^2;
            if sign(RBverts) == 1   % at least one corner is in domain
                count = count + 1;
                pntInterData(i) = 1;
            else
                pntInterData(i) = 0;
            end
        end
        
        if count == 3
            depFlag = 1;
            disp('inside')
        elseif count == 1 || count == 2
            depFlag = 2;
            disp('intersection with pnt inside')
        else
            depFlag = 3;
            disp('pure edge intersection, or outside')
        end
    end
else
    disp('outside')
end

end
