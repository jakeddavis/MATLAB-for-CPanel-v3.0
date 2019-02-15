
% Jake Davis
% 09/28/2018

function normal = getNorm(verts)

normal = cross(verts(2,:)-verts(1,:), verts(3,:)-verts(1,:))';
normal = normal / norm(normal);  % normalize

end