
% Jake Davis
% 03/02/2018

function C_MN = C_MN_func(M,N,verts,P_vec)

G_MN = zeros(1,3);
for i = 1:3
    pnt1 = [verts(i,1) verts(i,2)];
    pnt2 = [verts(i+1,1) verts(i+1,2)];
    %         x1 = pnt1(1); y1 = pnt1(2);
    %         x2 = pnt2(1); y2 = pnt2(2);
    
    triGeom = triGeom_Johnson(pnt1,pnt2,P_vec);
    
    G_MN(i) = G_MN_func(pnt1,pnt2,M,N,triGeom);
end

a = triGeom.a;
C_MN = (1/(M+N)) * a*sum(G_MN);

end