
% Jake Davis
% 10/18/2018

function trans = global2wind(a,b)

trans = [cosd(a)*cosd(b) -sind(b) sind(a)*cosd(b);...
         cosd(a)*sind(b) cosd(b) sind(a)*sind(b);...
         -sind(a) 0 cosd(a)];

end