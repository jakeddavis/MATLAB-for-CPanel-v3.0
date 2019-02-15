
% Jake Davis
% 10/30/2018

function pnt = supMinConePanelInter(verts,P,cntr,M)

n = getNorm(verts);
transLocal = supGetLocalSys(verts);
% % transLocal = getLocalSys(verts,[1 0 0]);

P = transLocal * (P - cntr)';
verts(1,:) = transLocal * (verts(1,:) - cntr)';
verts(2,:) = transLocal * (verts(2,:) - cntr)';
verts(3,:) = transLocal * (verts(3,:) - cntr)';
hold on
[~,Q,~] = triParams(verts,P,1);

newFS = transLocal * [1 0 0]';

w = P - Q';
if dot(n,w) > 0
    u = [sqrt(1-(1/M)^2) 0 1/M]
elseif dot(n,w) < 0
    u = [sqrt(1-(1/M)^2) 0 -1/M]
else
    error('Uh oh...')
end

u = transLocal * u'
u = u/norm(u)
s = -dot(n,w) / dot(n,u);

pnt = P + s*u;
plot3(pnt(1),pnt(2),pnt(3),'m*')

end