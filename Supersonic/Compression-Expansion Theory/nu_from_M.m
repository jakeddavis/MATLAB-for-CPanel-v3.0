% Jake Davis, HW 2, Problem 1, AERO 405-01

% Calculates the Prandtl-Meyer angle for a given Mach number

function nu = nu_from_M( M,gam )

   nu = sqrt((gam+1)/(gam-1))*atand(sqrt((M^2-1)*(gam-1)/(gam+1)))-atand(sqrt(M^2-1));

end

