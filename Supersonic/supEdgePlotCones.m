
% Jake Davis
% 01/01/2019

function supEdgePlotCones(pnt1,pnt2,max1,max2,Mach)

figure
hold on
edgeVecX = [pnt1(1) pnt2(1)];
edgeVecY = [pnt1(2) pnt2(2)];
plot(edgeVecX,edgeVecY,'-ok','LineWidth',2,'MarkerSize',7,'MarkerFaceColor','w','MarkerEdgeColor','k')
[~,~,~,xiVecP_1,etaVecP1_1,etaVecP2_1] = ...
    supPlotMachCone2(pnt1(1),max1,500,pnt1,Mach,0,0);
[~,~,~,xiVecP_2,etaVecP1_2,etaVecP2_2] = ...
    supPlotMachCone2(pnt2(1),max2,500,pnt2,Mach,0,0);
plot(xiVecP_1,etaVecP1_1,'-b','LineWidth',2)
plot(xiVecP_1,etaVecP2_1,'-b','LineWidth',2)
plot(xiVecP_2,etaVecP1_2,'-r','LineWidth',2)
plot(xiVecP_2,etaVecP2_2,'-r','LineWidth',2)

end