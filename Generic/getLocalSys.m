
% Jake Davis
% 09/28/2018

function local = getLocalSys(verts,FSdir)

p1 = verts(1,:);
Q = [mean(verts(:,1)),mean(verts(:,2)),mean(verts(:,3))];
normal = getNorm(verts);

unit = zeros(1,3);
local = zeros(3,3);
if FSdir == 1
    unit = [1 0 0];
else
    for i = 1:3
        unit(i) = p1(i) - Q(i);
    end
end
local(1,:) = unit/norm(unit);
local(2,:) = cross(normal,local(1,:));
local(3,:) = normal;

end