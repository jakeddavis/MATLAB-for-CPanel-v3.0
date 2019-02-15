function [ M2 ] = M2_from_M1( M1,gam )

   M2 = sqrt((1 + ((gam-1)/2)*M1^2)/(gam*M1^2 - (gam-1)/2));

end

