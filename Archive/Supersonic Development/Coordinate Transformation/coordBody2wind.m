
% Jake Davis
% 04/20/2018   aayyyooooo

function windVec = coordBody2wind(a,b,bodyVec)

% a = angle of attack in degrees
% b = angle of sideslip in degrees

a1 = cosd(a)*cosd(b);
a2 = cosd(a)*sind(b);
a3 = sind(a);
a4 = -sind(b);
a5 = cosd(b);
a6 = 0;
a7 = -sind(a)*cosd(b);
a8 = -sind(a)*sind(b);
a9 = cosd(a);
Ar = [a1 a2 a3; a4 a5 a6; a7 a8 a9];

windVec = Ar .* bodyVec';

end
