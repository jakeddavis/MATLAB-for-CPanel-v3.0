
% Jake Davis
% 05/25/2018

function [u,v,w] = subConstantDoubletNearFieldVel(mu,verts,P)

x = P(1); y = P(2); z = P(3);
verts(4,:) = verts(1,:);

uStuff = 0;
vStuff = 0;
wStuff = 0;
for i = 1:3
    pnt1 = [verts(i,1) verts(i,2)];
    pnt2 = [verts(i+1,1) verts(i+1,2)];
    x1 = pnt1(1); y1 = pnt1(2);
    x2 = pnt2(1); y2 = pnt2(2);
    
    triGeom = triGeom_func(pnt1,pnt2,P);
    r1 = triGeom.r1; r2 = triGeom.r2;
%     d1 = triGeom.d1; d2 = triGeom.d2;

%     u(i) = (z*(y1-y2)*(r1+r2)) / (r1*r2*(r1*r2-((x-x1)*(x-x2)+(y-y1)*(y-y2)+z^2)));
%     v(i) = (z*(x2-x1)*(r1+r2)) / (r1*r2*(r1*r2-((x-x1)*(x-x2)+(y-y1)*(y-y2)+z^2)));
%     w(i) = (((x-x2)*(y-y1)-(x-x1)*(y-y2))*(r1+r2)) / (r1*r2*(r1*r2-((x-x1)*(x-x2)+(y-y1)*(y-y2)+z^2)));

    uStuff = uStuff + (z*(y1-y2)*(r1+r2)) / (r1*r2*(r1*r2-((x-x1)*(x-x2)+(y-y1)*(y-y2)+z^2)));
    vStuff = vStuff + (z*(x2-x1)*(r1+r2)) / (r1*r2*(r1*r2-((x-x1)*(x-x2)+(y-y1)*(y-y2)+z^2)));
    wStuff = wStuff + (((x-x2)*(y-y1)-(x-x1)*(y-y2))*(r1+r2)) / (r1*r2*(r1*r2-((x-x1)*(x-x2)+(y-y1)*(y-y2)+z^2)));
end

% u = (mu/(4*pi)) * sum(u);
% v = (mu/(4*pi)) * sum(v);
% w = (mu/(4*pi)) * sum(w);

u = (mu/(4*pi)) * uStuff;
v = (mu/(4*pi)) * vStuff;
w = (mu/(4*pi)) * wStuff;

end
