
% Jake Davis
% 10/30/2018

function interFlag = supConeTriInter(verts,cntr,P,M,c0,nDir)

% Panel normal
n = getNorm(verts);
% n = nDir * getNorm(verts);

interPnt = supMinConePanelInter(verts,P,cntr,M,c0,nDir);
% plot3(interPnt(1),interPnt(2),interPnt(3),'g*')

AB = verts(2,:) - verts(1,:);
AC = verts(3,:) - verts(1,:);
% AP = P - verts(1,:);
AP = interPnt - verts(1,:);

gamma = (dot(cross(AB,AP),n)) / dot(n,n);
beta = (dot(cross(AP,AC),n)) / dot(n,n);
alpha = 1 - gamma - beta;

if (alpha >= 0) && (beta >= 0) && (gamma >= 0)
    interFlag = 1;
else
%     if interPnt(1) > cntr(1)
    if dot(interPnt-cntr,c0) >= 0
        interFlag = 1;
    else
        interFlag = 0;
    end
end

end

%% Replaced by a function

% w = P - cntr;
% if dot(n,w) > 0
%     u = [sqrt(1-(1/M)^2) 0 1/M];
% elseif dot(n,w) < 0
%     u = [sqrt(1-(1/M)^2) 0 -1/M];
% else
%     error('Uh oh...')
% end
% s = -dot(n,w) / dot(n,u);
% 
% interPnt = P + s*u;
% plot3(interPnt(1),interPnt(2),interPnt(3),'m*')

%% MATLAB's bary function testing

% vertIDs = [1 2 3];
% tri3d = triangulation(vertIDs,verts);
% B = cartesianToBarycentric(tri3d,1,interPnt)


%%

% u = [sqrt(1-(1/M)^2) 0 1/M];
% ut = [sqrt(1-(1/M)^2) 0 -1/M];
