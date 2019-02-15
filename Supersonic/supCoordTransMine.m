
% Jake Davis
% 11/29/2018

function transMat = supCoordTransMine(verts,M,a,b,nDir)

%%%% NOTE
% verts = panel vertices in 3-D
% M = freestream Mach number
% a = angle of attack in degrees
% b = angle of sideslip in degrees

% Super-/Sub-sonic switch
s = sign(1-M^2);
% If M > 1, s = -1, so beta = B
beta = sqrt(s*(1-M^2));

% Freestream direction in reference CSYS
c_r = [cosd(a)*cosd(b) -sind(b) sind(a)*cosd(b)]';

% Reference to wind transformation (aka compressible)
Gam_c = [cosd(a)*cosd(b) -sind(b) sind(a)*cosd(b);...
         cosd(a)*sind(b) cosd(b) sind(a)*sind(b);...
         -sind(a) 0 cosd(a)];

% Scaling for Mach number stuff
C = [1 0 0; 0 s*beta^2 0; 0 0 s*beta^2];
B = [s*beta^2 0 0; 0 1 0; 0 0 1];

% Combine Mach number scaling with wind-dir transformation
C_r = Gam_c' * C * Gam_c;

% Panel normal in reference CSYS
n_r = nDir * getNorm(verts);

% Other stuff
v_r = cross(n_r,c_r) / norm(cross(n_r,c_r)); % = transLocal(2,:)
r = sign(dot(n_r, B*n_r)); % SUBINCLINED: r = 1 % superinclined r = -1
u_r = cross(v_r,n_r); % = transLocal(1,:)

% Columns of A transpose
% col1 = (1/sqrt(abs(dot(n_r, B*n_r))))*C_r*u_r;
% col2 = (r*s/beta) * C_r*v_r;
% col3 = beta*n_r/sqrt(abs(dot(n_r, B*n_r)));

col1 = C_r*u_r;
col2 = (r*s/beta) * C_r*v_r;
col3 = beta*n_r;

transMat = [col1, col2, col3]';

end

%%

% ref2wind = [cosd(a)*cosd(b) -sind(b) sind(a)*cosd(b);...
%          cosd(a)*sind(b) cosd(b) sind(a)*sind(b);...
%          -sind(a) 0 cosd(a)];
% 
% beta = sqrt(M^2 - 1);
% B = [1 0 0; 0 beta 0; 0 0 beta];
% 
% transLocal = supGetLocalSys(verts,[1 0 0],nDir);
% 
% transMat = transLocal * ref2wind * B;


%%

% transLocal = supGetLocalSys(verts,[1 0 0],nDir);
