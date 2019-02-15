function M = Mach_from_nu_CPG( nu )

   nu_max = 90*(sqrt(6)-1);
   y = (nu/nu_max)^(2/3);

   A = 1.3604;
   B = 0.0962;
   C = -0.5127;
   D = -0.6722;
   E = -0.3278;
   
   M = (1 + A*y + B*y^2 + C*y^3)/(1 + D*y + E*y^2);

end

