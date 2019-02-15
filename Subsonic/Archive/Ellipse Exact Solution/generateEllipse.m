
% Jake Davis
% 01/09/2018

function generateEllipse(a,b,c,h0)

%% DIST MESH

% Generate points and connectivity
fd = @(p) p(:,1).^2/a^2 + p(:,2).^2/b^2 + p(:,3).^2/c^2 - 1; % Distance function
fh = @huniform; % mesh density function
bbox = [-a-.3,-b-.3,-c-.3; a+.3,b+.3,c+.3]; % region definition
[nodes,conn] = distmeshsurface(fd,fh,h0,bbox);
[nNodes,~] = size(nodes);
[nPanels,~] = size(conn);

% Write node points
fnamePrefix = ['ellipse_',num2str(nPanels),'tri'];

nodesFname = [fnamePrefix,'_nodes.txt'];
nodesFid = fopen(nodesFname,'w');
for i = 1:nNodes
    fprintf(nodesFid,'%12.8f %12.8f %12.8f\r\n',nodes(i,1),nodes(i,2),nodes(i,3));
end
fclose(nodesFid);
% Write connectivity
connFname = [fnamePrefix,'_elements.txt'];
connFid = fopen(connFname,'w');
for i = 1:nPanels
    fprintf(connFid,'%i %i %i\r\n',conn(i,1),conn(i,2),conn(i,3));
end
fclose(connFid);

triangulation_orient(fnamePrefix)
table_merge ( nodesFname, .001 )

%% WRITE TRI FILE

connNewFid = fopen([fnamePrefix,'_orient_elements.txt'],'r');
connNew = fscanf(connNewFid,'%i %i %i',[3 nPanels])';

triFid = fopen([fnamePrefix,'.tri'],'w');
% write nNodes and nPanels
fprintf(triFid,'%i %i\r\n',nNodes,nPanels);
% write node points
for i = 1:nNodes
    fprintf(triFid,'\t%12.8f \t%12.8f \t%12.8f\r\n',nodes(i,1),nodes(i,2),nodes(i,3));
end
% write connectivity
for i = 1:nPanels
    fprintf(triFid,'%i %i %i\r\n',connNew(i,1),connNew(i,2),connNew(i,3));
end
% write surface IDs (all 1's)
for i = 1:nPanels
    fprintf(triFid,'%i\r\n',1);
end

fclose(triFid);
fclose(connNewFid);

end
