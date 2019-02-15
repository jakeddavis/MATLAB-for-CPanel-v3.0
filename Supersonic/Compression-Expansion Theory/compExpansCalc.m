
% Jake Davis
% 11/16/2018

function [Cp2, Cp3] = compExpansCalc(M1, v1, alpha, theta1, theta2)

a = v1 / M1;    % ft/s
gam = 1.4;

theta1 = alpha - theta1;

% Compression turn
[ ~,~,~,M2,P2_P1,T2_T1] = Oblique_Solver( M1,theta1,gam );
Cp2 = (2 / (gam*M1^2)) * (P2_P1 - 1);

% v2 = a * M2;
% % Cp2 = 1 - (v2/v1)^2;
% Cp2 = 2 * (theta1*pi/180) / sqrt(M2^2 - 1);

% Expansion turn
[ ~,~,M3,~,~ ] = Expansion_Solver( M2,theta2,gam );
v3 = a * M3;
% Cp3 = 1 - (v3/v1)^2;
Cp3 = 2 * (-theta1*pi/180) / sqrt(M3^2 - 1);

% Display
format long
fprintf('Freestream Paramaters:\n')
fprintf(['M1 = ',num2str(M1), '\n'])
fprintf(['v1 = ',num2str(v1), '\n'])

fprintf('\n')
fprintf('Compression Turn:\n')
fprintf(['M2 = ',num2str(M2), '\n'])
fprintf(['v2 = ',num2str(v2), '\n'])
fprintf(['Cp2 = ',num2str(Cp2), '\n'])

fprintf('\n')
fprintf('Expansion Turn:\n')
fprintf(['M3 = ',num2str(M3), '\n'])
fprintf(['v3 = ',num2str(v3), '\n'])
fprintf(['Cp3 = ',num2str(Cp3), '\n'])
format short

end