
% Jake Davis
% 12/02/2018

function transMat = supCoordTransPilot(verts,M,a,b,nDir)

%%%% NOTE
% verts = panel vertices in 3-D
% M = freestream Mach number
% a = angle of attack in degrees
% b = angle of sideslip in degrees

% Reference to wind transformation (aka compressible)
Gam_c = [cosd(a)*cosd(b) -sind(b) sind(a)*cosd(b);...
         cosd(a)*sind(b) cosd(b) sind(a)*sind(b);...
         sind(a) 0 cosd(a)];

% Panel normal in reference CSYS
n_r = nDir * getNorm(verts);

Bt = sqrt(M^2-1);

s = -1;
nx_c = dot(Gam_c(1,:),n_r);
ny_c = dot(Gam_c(2,:),n_r);
nz_c = dot(Gam_c(3,:),n_r);
% n_c = [nx_c ny_c nz_c];
gam = sqrt(ny_c^2 + nz_c^2);
nc_wind = [-Bt^2*nx_c, ny_c nz_c];
den = sqrt(norm(dot(n_r,nc_wind)));

% nct = Bt^2*nx_c^2 + ny_c^2 + nz_c^2;

% nc_wind = [(-B*nx_c)^2 ny_c nz_c];
% sub-/super-inclined check
% dot(n_c,nc_wind)

A1 = [1, 0, 0; 0, Bt*nz_c/gam, -Bt*ny_c/gam; 0, Bt*ny_c/gam, Bt*nz_c/gam];
A2 = [gam/den, 0, -s*Bt*nx_c/den; 0, 1, 0; Bt*nx_c/den, 0, gam/den];
At = A2 * A1 * Gam_c;
sin(dot(n_r,nc_wind))

transMat = At;

end

%%

% % Super-/Sub-sonic switch
% s = sign(1-M^2);
% % If M > 1, s = -1, so beta = B
% beta = sqrt(s*(1-M^2));
% 
% % Freestream direction in reference CSYS
% c_r = [cosd(a)*cosd(b) -sind(b) sind(a)*cosd(b)]';
% 
% % Reference to wind transformation (aka compressible)
% Gam_c = [cosd(a)*cosd(b) -sind(b) sind(a)*cosd(b);...
%          cosd(a)*sind(b) cosd(b) sind(a)*sind(b);...
%          -sind(a) 0 cosd(a)];
% 
% % Scaling for Mach number stuff
% C = [1 0 0; 0 s*beta^2 0; 0 0 s*beta^2];
% B = [s*beta^2 0 0; 0 1 0; 0 0 1];
% 
% % Combine Mach number scaling with wind-dir transformation
% C_r = Gam_c' * C * Gam_c;
% 
% % Panel normal in reference CSYS
% n_r = nDir * getNorm(verts);
% 
% % Other stuff
% v_r = cross(n_r,c_r) / norm(cross(n_r,c_r)); % = transLocal(2,:)
% r = sign(dot(n_r, B*n_r)); % SUBINCLINED: r = 1 % superinclined r = -1
% u_r = cross(v_r,n_r); % = transLocal(1,:)
% 
% % Columns of A transpose
% % col1 = (1/sqrt(abs(dot(n_r, B*n_r))))*C_r*u_r;
% % col2 = (r*s/beta) * C_r*v_r;
% % col3 = beta*n_r/sqrt(abs(dot(n_r, B*n_r)));
% 
% col1 = (1/sqrt(abs(dot(n_r, B*n_r))))*C_r*u_r;
% col2 = (r*s/beta) * C_r*v_r;
% col3 = beta*n_r/sqrt(abs(dot(n_r, B*n_r)));
% 
% A = [col1, col2, col3]';
% 
% nxCo = dot(Gam_c(1,:),n_r);
% areaCorrect = beta / sqrt(1 - M^2*nxCo^2);
