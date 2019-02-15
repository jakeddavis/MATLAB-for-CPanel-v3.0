
% Jake Davis
% 02/03/2018

function [e1,e2,e3] = eFunc(verts,x,z)

e = zeros(length(verts),1);
for i = 1:length(verts)
    e(i) = (x-verts(i,1))^2 + z^2;
end

e1 = e(1);
e2 = e(2);
e3 = e(3);

end