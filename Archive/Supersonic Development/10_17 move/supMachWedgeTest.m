
% Jake Davis
% 10/16/2018

function edgeFlags = supMachWedgeTest(verts,Mach,P,cntr)

transMat = supCoordTrans(verts,Mach,0,0);
PSup = transMat * (P - cntr)';
newVerts(1,:) = transMat * (verts(1,:) - cntr)';
newVerts(2,:) = transMat * (verts(2,:) - cntr)';
newVerts(3,:) = transMat * (verts(3,:) - cntr)';
[~,~,~] = triParams(newVerts,PSup,1);
supPlotMachCone(-1,PSup(1),100,PSup,Mach,1)

% newVerts = verts;
% PSup = P;

x = PSup(1); y = PSup(2); z = PSup(3);
newVerts(end+1,:) = newVerts(1,:);
for i = 1:3
    pnt1 = newVerts(i,:);
    pnt2 = newVerts(i+1,:);
    x1 = pnt1(1); y1 = pnt1(2);
    x2 = pnt2(1); y2 = pnt2(2);
    m  = (y2-y1) / (x2-x1);
    xm = (x-x1) - (y-y1)/m;
    xmc = m*xm;
    if abs(m) < 1
        zm = xmc^2 + z^2;
    elseif abs(m) > 1
        zm = xm^2 + z^2/m^2;
    end
    stuff = zm - z^2
    if zm - z^2 >= 0
        edgeFlags(i) = 1;
    else
        edgeFlags(i) = 0;
    end
end

end