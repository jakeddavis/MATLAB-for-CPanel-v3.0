
% Jake Davis
% 10/10/2018

function [Q1,w0] = supEdgeInfVertexSup(mySign,lam,R,xm,ym,s,sm,z,delSmall)

if R > delSmall
    Q1 = mySign * atan((-z*ym) / (xm*R));
    w0 = mySign * (1/sqrt(1-lam^2)) * atan(-ym / (R*sqrt(1-lam^2)));
elseif R < delSmall
    Q1 = sign(z*(s-lam*sm)) * pi/2;
    w0 = sign(s-lam*sm) * (pi/2)*sqrt(1-lam^2);
end

end
