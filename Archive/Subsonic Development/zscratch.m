
clc, clear, close all

% Create figure

figure
hold all; grid on

% Geometry definition

verts = [.2 .3 ; .5 0 ; .8 .5];
POI = [1 .3 .5];
[A,cntr,cond] = triParams(verts,POI);

pnt1 = [verts(1,:),0];
pnt2 = [verts(2,:),0];
pnt3 = [verts(3,:),0];
cntr = [cntr, 0];

plot3(cntr(1),cntr(2),cntr(3),'*k')
plot3([cntr(1) POI(1)],[cntr(2) POI(2)],[cntr(3) POI(3)],'r')

% Panel geometry

a = POI - pnt1;
b = POI - pnt2;
s = pnt2 - pnt1;
pjk = POI - cntr;
n = cross(pnt2-pnt1, pnt3-pnt1);
n = n/norm(n);
l = (pnt1 - cntr)/norm(pnt1 - cntr);    % x-dir
m = cross(l,n);

% Plot coord. system
plot3([0 l(1)], [0 l(2)], [0 l(3)], 'k')
plot3([0 m(1)], [0 m(2)], [0 m(3)], 'k')
plot3([0 n(1)], [0 n(2)], [0 n(3)], 'k')
% plot3([cntr(1) l(1)], [cntr(2) l(2)], [cntr(3) l(3)], 'k')
% plot3([cntr(1) m(1)], [cntr(2) m(2)], [cntr(3) m(3)], 'k')
% plot3([cntr(1) n(1)], [cntr(2) n(2)], [cntr(3) n(3)], 'k')

PN = dot(pjk,n);

as = cross(a,s); % normal to plane of s and a
las = cross(l,as);
plot3([cntr(1) as(1)],[cntr(2) as(2)],[cntr(3) as(3)],'m')
plot3([cntr(1) las(1)],[cntr(2) las(2)],[cntr(3) las(3)],'g')
Al = dot(n,cross(s,a));

% plot3([cntr(1) Al(1)],[cntr(2) Al(2)],[cntr(3) Al(3)],'k')


% legend
% legend('edge','normal','POI')
