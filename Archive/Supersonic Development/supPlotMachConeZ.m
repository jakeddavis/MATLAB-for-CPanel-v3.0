
% Jake Davis
% 10/06/2018

function supPlotMachConeZ(min,max,n,P,verts,Mach,trans)

xiVec = linspace(min,max,n);
etaVec1 = zeros(1,length(xiVec));
etaVec2 = zeros(1,length(xiVec));

if trans == 1
    for i = 1:length(xiVec)
        etaVec1(i) = P(2) + sqrt((xiVec(i)-P(1))^2 - P(3)^2);
        etaVec2(i) = P(2) - sqrt((xiVec(i)-P(1))^2 - P(3)^2);
    end
else
    B = sqrt(Mach^2-1);
    for i = 1:length(xiVec)
%         etaVec1(i) = P(2) + sqrt((xiVec(i)-P(1))^2/B^2 - P(3)^2);
%         etaVec2(i) = P(2) - sqrt((xiVec(i)-P(1))^2/B^2 - P(3)^2);
        etaVec1(i) = P(2) + sqrt((xiVec(i)-P(1))^2/B^2 - P(3)^2);
        etaVec2(i) = P(2) - sqrt((xiVec(i)-P(1))^2/B^2 - P(3)^2);
    end
end

hold on
plot(xiVec,etaVec1,'b','LineWidth',1)
plot(xiVec,etaVec2,'b','LineWidth',1)

end

%%

%         etaVec1(i) = P(2) + sqrt((xiVec(i)-P(1))^2 + P(3)^2);
%         etaVec2(i) = P(2) - sqrt((xiVec(i)-P(1))^2 + P(3)^2);

%% plot stuff with legend

% title('Reference CSYS','FontSize',16)
% xlabel('x')
% ylabel('y')
% xlim([-1 2])
% ylim([-2 2])
% lgdPlot(1) = plot(nan, nan, strcat('-r'),'LineWidth',1);
% lgdStr{1} = 'M = 3.0';
% legend(lgdPlot, lgdStr,'Location','southeast');
% set(legend,'FontSize',14,'Interpreter','latex')
