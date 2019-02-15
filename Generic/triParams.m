
% Jake Davis
% 02/03/2018

% Compute geometry characteristics

function [A,Q,cond] = triParams(verts,P,plotFlag)

x = P(1); y = P(2);
x1 = verts(1,1); y1 = verts(1,2);
x2 = verts(2,1); y2 = verts(2,2);
x3 = verts(3,1); y3 = verts(3,2);
triVerts = verts;
triVerts(:,3) = [];

% area of triangle
% A = abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))/2;
a = verts(2,:) - verts(1,:);
b = verts(3,:) - verts(1,:);
theta = acos(dot(a,b) / norm(a)*norm(b));
A = 0.5*norm(a)*norm(b)*sin(theta);

% find centroid (influencing point)
Q = [mean(verts(:,1)),mean(verts(:,2)),mean(verts(:,3))];
% Q = [mean(verts(:,1)),mean(verts(:,2))];

a = norm(verts(1,:));
b = norm(verts(2,:));
c = norm(verts(3,:));

verts_wrap = verts;
verts_wrap(4,:) = verts(1,:);
diams = zeros(3,1);
cond = zeros(1,2);
for i = 1:3
    pnt1 = verts_wrap(i,:);
    pnt2 = verts_wrap(i+1,:);
    x1 = pnt1(1); y1 = pnt1(2); z1 = pnt1(3);
    x2 = pnt2(1); y2 = pnt2(2); z2 = pnt2(3);
    
    % find triangle diameter
    diams(i) = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2);
%     test(i) = norm(pnt2 - pnt1)
end
cond(1)   = max(diams)/2;

vertIDs = [1 2 3];
tri2d = triangulation(vertIDs,triVerts);
P_bary_coords = cartesianToBarycentric(tri2d,1,[x y]);
if any(P_bary_coords(:)<=0) == 1
    % field point (x,y,0) outside of panel
    cond(2) = 0;
else
    % field point (x,y,0) inside of panel
    cond(2) = 1;
end

% Plot panel geometry
if plotFlag == 1
%     figure
    hold on
    tri3d = triangulation(vertIDs,verts);
    t = trisurf(tri3d);
    t.FaceColor = [0 .7 .7];
%     plot3(Q(1),Q(2),Q(3),'ok','MarkerSize',10,'MarkerFaceColor','k')
    plot3(Q(1),Q(2),Q(3),'ok','MarkerSize',5,'MarkerFaceColor','k')
    plot3(P(1),P(2),P(3),'sk','MarkerSize',10,'MarkerFaceColor','r')
%     lgdPlot(1) = plot(nan, nan, strcat('-r'),'LineWidth',1);
%     lgdStr{1} = 'M = 3.0';
%     legend(lgdPlot, lgdStr,'Location','northeast');
%     set(legend,'FontSize',14,'Interpreter','latex')
%     title('Geometry','FontSize',16)
%     xlabel('x','FontSize',13,'FontWeight','bold')
%     ylabel('y','FontSize',13,'FontWeight','bold')
%     zlabel('z','FontSize',13,'FontWeight','bold')
%     set(get(gca,'ylabel'),'rotation',0)
%     set(get(gca,'zlabel'),'rotation',0)
%     grid on
end
    
end

%%

% % Plot panel geometry
% if plotFlag == 1
%     figure
%     hold on
%     tri3d = triangulation(vertIDs,verts);
%     trisurf(tri3d)
%     % triplot(tri,'-ob','LineWidth',1.5)
%     plot3(Q(1),Q(2),Q(3),'ok','MarkerSize',10,'MarkerFaceColor','k')
%     % plot3(Q(1),Q(2),0,'*k')
%     plot3(P(1),P(2),P(3),'*k')
% %     lgd = legend('Panel Geometry','Panel Centroid','Field Point');
% %     lgd.FontSize = 12;
%     lgdPlot(1) = plot(nan, nan, strcat('-r'),'LineWidth',1);
%     lgdStr{1} = 'M = 3.0';
%     legend(lgdPlot, lgdStr,'Location','southeast');
%     set(legend,'FontSize',14,'Interpreter','latex')
% %     title('Geometry','FontSize',16)
% %     xlabel('x','FontSize',13,'FontWeight','bold')
% %     ylabel('y','FontSize',13,'FontWeight','bold')
% %     zlabel('z','FontSize',13,'FontWeight','bold')
% %     set(get(gca,'ylabel'),'rotation',0)
% %     set(get(gca,'zlabel'),'rotation',0)
%     axis equal
% %     grid on
% end

%%

%  'LineWidth',2,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor',[.49 1 .63],...
%     'MarkerSize',10)

%%

% dHs = zeros(3,1);
% dFs = zeros(3,1);
% cond = zeros(1,6);

%     % find min distance from P to the perimeter of the triangle
%     [dHs(i),dFs(i)] = project_point_to_line_segment(pnt1,pnt2,P);

% cond(3)   = min(dHs);
% cond(5)   = min(dFs);
% disp(dHs)

%%

% % find triangle diameter
% k = 1;
% [pts,~] = size(verts);
% diams = zeros(pts,1);
% for i = 1:pts-1
%     for j = 2:pts
%         if i ~= j
%             diams(k) = sqrt((verts(j,1) - verts(i,1))^2 + (verts(j,2) - verts(i,2))^2);
%             k = k + 1;
%         end
%     end
% end
    
% diam = max(diams);

%%

%     xd = verts(i,1); yd = verts(i,2);
%     k = i - 1;
%     dHs(i+k) = sqrt((xd-x)^2 + (yd-y)^2);
%     dHs(2*i) = abs((y2-y1)*x - (x2-x1)*y + x2*y1 - y2*x1) / ...
%         sqrt((y2-y1)^2 + (x2-x1)^2);
