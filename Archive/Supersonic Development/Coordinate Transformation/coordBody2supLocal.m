
% Jake Davis
% 04/20/2018   aayyyooooo

function transMat = coordBody2supLocal(verts,B,a,b)
% function localSupVec = coordBody2supLocal(verts,B,bodyVec,a,b,cntr)

% NOTE
% a = angle of attack in degrees
% b = angle of sideslip in degrees

[~,Q,~] = triParams(verts,[0 0 0],0);


% vertIDs = [1 2 3];
% tri = triangulation(vertIDs,verts);
% n_ref = faceNormal(tri, vertIDs);
% [nx,ny,nz] = surfnorm(verts(:,1), verts(:,2), verts(:,3));
% n_ref = [nx, ny, nz]';
n_ref = cross(verts(2,:)-verts(1,:), verts(3,:)-verts(2,:));
n_ref = n_ref / norm(n_ref);

Ar = [cosd(a)*cosd(b) cosd(a)*sind(b) sind(a);...
     -sind(b) cosd(b) 0;...
     -sind(a)*cosd(b) -sind(a)*sind(b) cosd(a)];

nx_c = dot(Ar(1,:),n_ref);
ny_c = dot(Ar(2,:),n_ref);
nz_c = dot(Ar(3,:),n_ref);

n_wind = [nx_c ny_c nz_c];
% nc_wind = [(-B*nx_c)^2 ny_c nz_c];
nc_wind = [-B^2*nx_c, ny_c nz_c];

hold on
plot3([Q(1) Q(1)+n_ref(1)], [Q(2) Q(2)+n_ref(2)], [Q(3) Q(3)+n_ref(3)])
plot3([Q(1) Q(1)+nc_wind(1)], [Q(2) Q(2)+nc_wind(2)], [Q(3) Q(3)+nc_wind(3)])

% sub-/super-inclined check
dot(n_ref,nc_wind)
test = -B^2*nx_c^2 + ny_c^2 + nz_c^2

gam = sqrt(ny_c^2 + nz_c^2);

A1 = [1 0 0;...
      0 B*nz_c/gam -B*ny_c/gam;...
      0 B*ny_c/gam B*nz_c/gam];

s = -1;
A2 = [gam/sqrt(abs(dot(n_wind,nc_wind))) 0 -s*B*nx_c/sqrt(abs(dot(n_wind,nc_wind)));...
      0 sin(dot(n_wind,nc_wind)) 0;...
      B*nx_c/sqrt(abs(dot(n_wind,nc_wind))) 0 gam/sqrt(abs(dot(n_wind,nc_wind)))];

transMat = A2 * A1 * Ar;
% localSupVec = A2 * A1 * Ar * (bodyVec-cntr)';

end

%%

% Bmat = [B^2 0 0; 0 1 0; 0 0 1];

% n = stuff;
% n_c_T = Ar \ Bmat * Ar * n';     % row vector
% n_c = n_c_T';
% nx_c = n_c(1); ny_c = n_c(2); nz_c = n_c(3);
