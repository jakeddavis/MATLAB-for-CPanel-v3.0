
% Jake Davis
% 11/02/2018

function pnt = supMinConePanelInter(verts,P,cntr,M,c0,nDir)

% Transformation matrix for panel on z=0
transLocal = supGetLocalSys(verts,c0,nDir);
Plocal = transLocal * (P - cntr)';

nOrig = nDir * getNorm(verts);

% vertsFlat(1,:) = transLocal * (verts(1,:) - cntr)';
% vertsFlat(2,:) = transLocal * (verts(2,:) - cntr)';
% vertsFlat(3,:) = transLocal * (verts(3,:) - cntr)';
% hold on
% [~,~,~] = triParams(vertsFlat,Plocal,1);
% plot3(Plocal(1),Plocal(2),Plocal(3),'m*')

z0Flag = 0;
h = dot(P - cntr,transLocal(3,:));
if sign(nOrig(3)) == sign(h)
    u = [-sqrt(1-(1/M)^2) 0 -1/M];
elseif sign(nOrig(3)) ~= sign(h)
    if nOrig(1) == 0 && nOrig(3) == 0
        z0Flag = 1;
    end
    u = [-sqrt(1-(1/M)^2) 0 1/M];
else
    error('Uh oh...')
end

% h = dot(P - cntr,transLocal(3,:));
% if h > 0
%     u = [-sqrt(1-(1/M)^2) 0 -1/M];
% elseif h < 0
%     u = [-sqrt(1-(1/M)^2) 0 1/M];
% else
%     error('Uh oh...')
% end

% Transform z-dir of u, and other stuff
if z0Flag == 0
    u(3) = transLocal(3,:) * u';
end
u = u / norm(u);
n = [0 0 1];
w = Plocal';
s = -dot(n,w) / dot(n,u);

% Intersection point in local CSYS
pntFlat = Plocal' + s*u;

% Intersection point in global CSYS
pnt = ((transLocal' * pntFlat') + cntr')';

% pause

end

%% Plot in global CSYS

% vertsOld(1,:) = transLocal' * vertsFlat(1,:)' + cntr';
% vertsOld(2,:) = transLocal' * vertsFlat(2,:)' + cntr';
% vertsOld(3,:) = transLocal' * vertsFlat(3,:)' + cntr';
% [~,~,~] = triParams(verts,P,1);
% plot3(pnt(1),pnt(2),pnt(3),'r*')

%% Plot in local CSYS

% vertsFlat(1,:) = transLocal * (verts(1,:) - cntr)';
% vertsFlat(2,:) = transLocal * (verts(2,:) - cntr)';
% vertsFlat(3,:) = transLocal * (verts(3,:) - cntr)';
% hold on
% [~,~,~] = triParams(vertsFlat,Pflat,1);
% plot3(pntFlat(1),pntFlat(2),pntFlat(3),'m*')

%% pre cleanup - 11/2

% n = getNorm(verts);
% transLocal = supGetLocalSys(verts);
% % % transLocal = getLocalSys(verts,[1 0 0]);
% 
% P = transLocal * (P - cntr)';
% verts(1,:) = transLocal * (verts(1,:) - cntr)';
% verts(2,:) = transLocal * (verts(2,:) - cntr)';
% verts(3,:) = transLocal * (verts(3,:) - cntr)';
% hold on
% [~,~,~] = triParams(verts,P,1);
% 
% newFS = transLocal * [1 0 0]';
% newFS = newFS / norm(newFS);
% 
% zTrans = eye(3);
% zTrans(3,:) = transLocal(3,:);
% 
% if P(3) > 0
%     u = [sqrt(1-(1/M)^2) 0 1/M]
% elseif P(3) < 0
%     u = [sqrt(1-(1/M)^2) 0 -1/M]
% else
%     error('Uh oh...')
% end
% 
% u(3) = transLocal(3,:) * u'
% u = u / norm(u)
% n = [0 0 1];
% w = P';
% 
% % u = transLocal * u'
% % u = u/norm(u)
% 
% s = -dot(n,w) / dot(n,u);
% 
% pnt = P' + s*u;
% 
% plot3(pnt(1),pnt(2),pnt(3),'m*')

%%

% n = getNorm(verts);
% transLocal = supGetLocalSys(verts);
% % % transLocal = getLocalSys(verts,[1 0 0]);
% 
% P = transLocal * (P - cntr)';
% verts(1,:) = transLocal * (verts(1,:) - cntr)';
% verts(2,:) = transLocal * (verts(2,:) - cntr)';
% verts(3,:) = transLocal * (verts(3,:) - cntr)';
% hold on
% [~,Q,~] = triParams(verts,P,1);
% 
% newFS = transLocal * [1 0 0]';
% 
% w = P - Q';
% if dot(n,w) > 0
%     u = [sqrt(1-(1/M)^2) 0 1/M]
% elseif dot(n,w) < 0
%     u = [sqrt(1-(1/M)^2) 0 -1/M]
% else
%     error('Uh oh...')
% end
% 
% u = transLocal * u'
% u = u/norm(u)
% s = -dot(n,w) / dot(n,u);
% 
% pnt = P + s*u;
% plot3(pnt(1),pnt(2),pnt(3),'m*')
