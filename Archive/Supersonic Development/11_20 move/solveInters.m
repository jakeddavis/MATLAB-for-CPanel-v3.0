
% Jake Davis
% 10/14/2018

syms m n l x y z xi xi1 eta eta1 zeta zeta1 B

% line = m*(xi-xi1) + eta1;
% hyper = y + sqrt((xi-x)^2/B^2 - z^2);
% % hyper = y + sqrt((xi-x)^2 - z^2);
% 
% s = solve(line == hyper, xi)


% Standard Solution

% line1 = m*(xi-xi1) + eta1;
% line2 = (zeta-zeta1)/n + eta1;
% hyper = y + sqrt((xi-x)^2/B^2 - (zeta-z)^2);
% eqn1 = line1 == hyper;
% eqn2 = line2 == hyper;
% [xiSol,zetaSol] = solve([eqn1 eqn2], [xi zeta])

% xi2 = xi1

% xi = xi1;
% line = (zeta-zeta1)/n + eta1;
% hyper = y + sqrt((xi1-x)^2/B^2 - (zeta-z)^2);
% eqn = line == hyper;
% zetaSol = solve(eqn, zeta)

% eta2 = eta1

% eta = eta1;
% line = (zeta-zeta1)/l + xi1;
% hyper = sqrt(B^2 * ((eta1-y)^2)+(zeta-z)^2);
% eqn = line == hyper;
% zetaSol = solve(eqn, zeta)

% zeta2 = zeta1

zeta = zeta1;
line = m*(xi-xi1) + eta1;
hyper = y + sqrt((xi-x)^2/B^2 - (zeta1-z)^2);
eqn = line == hyper;
xiSol = solve(eqn, xi)
