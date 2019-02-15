
% Jake Davis
% 11/01/2018

function transLocal = supGetLocalSys(verts,FS,nDir)

n = nDir * getNorm(verts)';
transLocal = zeros(3,3);

transLocal(2,:) = cross(n,FS) / norm(cross(n,FS));
transLocal(1,:) = cross(transLocal(2,:), n);
transLocal(3,:) = n;

% transLocal(2,:) = cross(FS,n) / norm(cross(FS,n));
% transLocal(1,:) = cross(n, transLocal(2,:));
% transLocal(3,:) = n;

end

%% supersonic transformation - projection onto z = 0

% n_r = n';
% c_r = [1 0 0]';
% C_r = [1 0 0; 0 -1 0; 0 0 -1];
% v_r = (cross(n_r,c_r) / norm(cross(n_r,c_r)));
% r = sign(dot(n_r,n_r));
% u_r = cross(v_r,n_r);
% s = -1;
% beta = 1;
% 
% 
% transOld = [(1/sqrt(abs(dot(n_r,n_r))))*C_r*u_r, r*s/beta * C_r*v_r ,...
%     beta*n_r/sqrt(abs(dot(n_r,n_r)))]'
% 
% transOldt = transOld'
