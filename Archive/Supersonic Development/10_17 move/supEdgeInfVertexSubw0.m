
% Jake Davis
% 10/10/2018

function w0 = supEdgeInfVertexSubw0(m,R,s,sm,delSmall)

w0 = 0;
if R > delSmall
    w0 = (m/(2*sqrt(1-m^2))) * log(((sm-m*s+R*sqrt(1-m^2))) / ...
        (sm-m*s-R*sqrt(1-m^2)));
end

end
