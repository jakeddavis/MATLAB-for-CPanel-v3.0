
% Jake Davis
% 02/03/2017

function [m12,m23,m31] = mFunc(verts)

k = 1;
[pts,~] = size(verts);
m = zeros(pts,1);
for i = 1:pts-1
    for j = 2:pts
        if i ~= j
            m(k) = (verts(j,2) - verts(i,2)) / (verts(j,1) - verts(i,1));
%             m(k) = sqrt((verts(j,1) - verts(i,1))^2 + (verts(j,2) - verts(i,2))^2);
            k = k + 1;
        end
    end
end

m12 = m(1);
m31 = m(2);
m23 = m(3);

end