
% Jake Davis
% 10/30/2018

function pnt = supMinConePanelInterOrig(verts,P,cntr,M)

n = getNorm(verts);
w = P - cntr;
if dot(n,w) > 0
    u = [sqrt(1-(1/M)^2) 0 1/M];
elseif dot(n,w) < 0
    u = [sqrt(1-(1/M)^2) 0 -1/M];
else
    error('Uh oh...')
end

s = -dot(n,w) / dot(n,u);

pnt = P + s*u;
% plot3(pnt(1),pnt(2),pnt(3),'m*')

end