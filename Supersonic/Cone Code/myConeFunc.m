
% Jake Davis
% 01/07/2019

function myConeFunc(origin,length,theta,phiSpace,lSpace,writeFile,varargin)

% origin - cone tip
% length - cone length
% theta - cone half angle
% phi - spacing between circular points
% lSpace - spacing along length of cone

%% BUILD CONTROLING LINES OF MESH

% Find streamwise mesh divisions
spacing = length / lSpace;
for i = 1:lSpace
    circles.Locs(i) = origin(1) + i*spacing;
    circles.Diams(i) = 2 * circles.Locs(i) * tand(theta);
end

%% COMPUTE AND NUMBER NODAL POINTS

figure
hold on
plot3(origin(1),origin(2),origin(3),'*k')
plot3(origin(1)+length,origin(2),origin(3),'*k')

% nCircPnts = abs(360/phi);
nCircPnts = phiSpace;
phi = -360/phiSpace;
k = 1;
% starting point (origin)
nodeStart.Pnt = origin;
nodeStart.Index = k;
for i = 1:lSpace
    xCenter = origin(1) + i*spacing;
    for j = 1:nCircPnts
        % cell where each row is a circle, and each elements is a point
        circPnts{i,j} = [xCenter cosd(j*phi)*circles.Diams(i)/2 sind(j*phi)*circles.Diams(i)/2];
        plot3(xCenter,cosd(j*phi)*circles.Diams(i)/2,sind(j*phi)*circles.Diams(i)/2,'*k')
        
        k = k + 1;
        node.Pnt = circPnts{i,j};
        node.Index = k;
        nodesCell{i,j} = node;
    end
end
axis equal
% ending point (center of cone base)
nodeEnd.Pnt = [origin(1)+length,origin(2),origin(3)];
nodeEnd.Index = k + 1;

%% BUILD PANEL CONNECTIVITY

% origin to first division
i = 1;
while i < nCircPnts
    connectivity{i} = [nodeStart.Index, nodesCell{1,i}.Index, nodesCell{1,i+1}.Index];
    i = i + 1;
end
connectivity{i} = [nodeStart.Index, nodesCell{1,end}.Index, nodesCell{1,1}.Index];

% remaining regions (continue using i for connectivity index)
i = i + 1;
for j = 2:lSpace
    k = 1;
    while k < nCircPnts
        connectivity{i} = [nodesCell{j-1,k}.Index, nodesCell{j,k}.Index, nodesCell{j,k+1}.Index];
        connectivity{i+1} = [nodesCell{j-1,k}.Index, nodesCell{j,k+1}.Index, nodesCell{j-1,k+1}.Index];
        i = i + 2;
        k = k + 1;
    end
    connectivity{i} = [nodesCell{j-1,k}.Index, nodesCell{j,k}.Index, nodesCell{j,1}.Index];
    connectivity{i+1} = [nodesCell{j-1,k}.Index, nodesCell{j,1}.Index, nodesCell{j-1,1}.Index];
    i = i + 2;
end

% % cone base
% k = 1;
% while k < nCircPnts
%     connectivity{i} = [nodeEnd.Index, nodesCell{end,k}.Index, nodesCell{end,k+1}.Index];
%     k = k + 1;
%     i = i + 1;
% end
% connectivity{i} = [nodeEnd.Index, nodesCell{end,k}.Index, nodesCell{end,1}.Index];

%% BUILD MATRICES FOR OUTPUT AND PLOT

% triangulation matrices for origin to first division
[~,nPanels] = size(connectivity);
for i = 1:nPanels
    outCon(i,:) = connectivity{i};
end
outPnt(1,:) = origin;
k = 1;
for i = 1:lSpace
    for j = 1:nCircPnts
        outPnt(k+1,:) = nodesCell{i,j}.Pnt;
        k = k + 1;
    end
end
% outPnt(end+1,:) = nodeEnd.Pnt;
[nNodes,~] = size(outPnt);

tri3d = triangulation(outCon,outPnt);
trisurf(tri3d)

%% WRITE TRI FILE

if writeFile == true
    if size(varargin) > 0
        fNameIn = varargin{1};
        filename = fNameIn;
    else
        filename = ['cone_',num2str(nPanels),'tri.tri'];
    end
    
    fid = fopen(filename,'w');
    % write nNodes and nPanels
    fprintf(fid,'%i %i\r\n',nNodes,nPanels);
    % write node points
    for i = 1:nNodes
        fprintf(fid,'\t%12.8f \t%12.8f \t%12.8f\r\n',outPnt(i,1),outPnt(i,2),outPnt(i,3));
    end
    % write connectivity
    for i = 1:nPanels
        fprintf(fid,'%i %i %i\r\n',outCon(i,1),outCon(i,2),outCon(i,3));
    end
    % write surface IDs (all 1's)
    for i = 1:nPanels
        fprintf(fid,'%i\r\n',1);
    end
    
    fclose(fid);
end

% if size(varargin) > 0
%     fNameIn = varargin{1};
%     filename = fNameIn;
% else
%     filename = ['cone_',num2str(nPanels),'tri.tri'];
% end

% fid = fopen(filename,'w');
% % write nNodes and nPanels
% fprintf(fid,'%i %i\r\n',nNodes,nPanels);
% % write node points
% for i = 1:nNodes
%     fprintf(fid,'\t%12.8f \t%12.8f \t%12.8f\r\n',outPnt(i,1),outPnt(i,2),outPnt(i,3));
% end
% % write connectivity
% for i = 1:nPanels
%     fprintf(fid,'%i %i %i\r\n',outCon(i,1),outCon(i,2),outCon(i,3));
% end
% % write surface IDs (all 1's)
% for i = 1:nPanels
%     fprintf(fid,'%i\r\n',1);
% end
% 
% fclose(fid);

end







