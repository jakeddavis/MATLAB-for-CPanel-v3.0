
% Jake Davis
% 02/03/2018

function [r1,r2,r3] = rFunc(verts,x,y,z)

r = zeros(length(verts),1);
for i = 1:length(verts)
    r(i) = sqrt((x-verts(i,1))^2 + (y-verts(i,2))^2 + z^2);
end

r1 = r(1);
r2 = r(2);
r3 = r(3);

end