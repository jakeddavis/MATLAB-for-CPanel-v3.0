% build triangle
verts = [0 .25 ; 1 0 ; cosd(60) sind(60)];
[A,Q,diam] = triParams(verts);
P = [1 -5 1];
muVec = [1 5 0];
% mu = mean(muVec);
% mu = muVec(1);
mu = 10;
sigma = 1;

% move panel local coordinate origin to center of triangle, in x-y plane
verts(:,1) = verts(:,1) - Q(1);
verts(:,2) = verts(:,2) - Q(2);
P(1) = P(1) - Q(1);
P(2) = P(2) - Q(2);
% recompute triparams
[~,Q,~] = triParams(verts);

[mu] = bary(verts,muVec)