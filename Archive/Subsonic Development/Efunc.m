
% Jake Davis
% 04/02/2018

function E = Efunc(M,N,K,pnt1,pnt2,P)

x = P(1); y = P(2); z = P(3);
x1 = pnt1(1,1); y1 = pnt1(1,2);
x2 = pnt2(1,1); y2 = pnt2(1,2);

rho2 = sqrt((x2-x)^2 + (y2-y)^2 + z^2);
rho1 = sqrt((x1-x)^2 + (y1-y)^2 + z^2);

E2 = ((x2-x)^(M-1) * (y2-y)^(N-1)) / rho2^(K-2);
E1 = ((x1-x)^(M-1) * (y1-y)^(N-1)) / rho1^(K-2);

E = E2 - E1;

end
