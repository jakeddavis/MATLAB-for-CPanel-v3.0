
% Jake Davis
% 02/03/2018

function [h1,h2,h3] = hFunc(verts,x,y)

h = zeros(length(verts),1);
for i = 1:length(verts)
    h(i) = (x-verts(i,1)) * (y-verts(i,2));
end

h1 = h(1);
h2 = h(2);
h3 = h(3);

end