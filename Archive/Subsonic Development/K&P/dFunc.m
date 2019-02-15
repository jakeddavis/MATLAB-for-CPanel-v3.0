
% Jake Davis
% 02/03/2018

% Compute d in source influence near field equation

function [d12,d23,d31] = dFunc(verts)

k = 1;
[pts,~] = size(verts);
d = zeros(pts,1);
for i = 1:pts-1
    for j = 2:pts
        if i ~= j
            d(k) = sqrt((verts(j,1) - verts(i,1))^2 + (verts(j,2) - verts(i,2))^2);
            k = k + 1;
        end
    end
end

d12 = d(1);
d31 = d(2);
d23 = d(3);

end