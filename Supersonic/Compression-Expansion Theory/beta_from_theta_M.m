
function beta = beta_from_theta_M( theta, M, del, gamma )

lam = sqrt((M^2-1)^2-3*(1+((gamma-1)/2)*M^2)*(1+((gamma+1)/2)*M^2)*(tand(theta)^2));

chi = ((M^2-1)^3-9*(1+((gamma-1)/2)*M^2)*(1+((gamma-1)/2)*M^2+((gamma+1)/4)*M^4)*(tand(theta)^2))/(lam^3);

beta = atand((M^2-1+2*lam*cos((4*pi*del+acos(chi))/3))/(3*(1+M^2*(gamma-1)/2)*tand(theta)));

end

