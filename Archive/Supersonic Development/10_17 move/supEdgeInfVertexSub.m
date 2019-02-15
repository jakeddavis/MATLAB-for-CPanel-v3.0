
% Jake Davis
% 10/10/2018

function [Q1,w0] = supEdgeInfVertexSub(mySign,m,R,xmc,ymc,s,sm,z,delSmall)

Q1 = 0;
w0 = 0;
if R > delSmall
    Q1 = mySign * sign(z) * atan((xmc*R) / (abs(z)*ymc));
    w0 = mySign * (m/(2*sqrt(1-m^2))) * log(((sm-m*s+R*sqrt(1-m^2))) / ...
        (sm-m*s-R*sqrt(1-m^2)));
end

end
