
% Jake Davis
% 10/06/2018

function supPlotMachCone(min,max,n,P,Mach,trans,varargin)

xiVecP = linspace(min,max,n);
% etaVec1 = zeros(1,length(xiVec));
% etaVec2 = zeros(1,length(xiVec));
zVec = zeros(1,length(xiVecP));
zVec(1,:) = P(3);

if trans == 1
%     pntMin = [P(1) - sqrt(P(3)^2), P(2), 0];
    xiMax = P(1) - sqrt(P(3)^2);
    xiVec = linspace(min,xiMax,n);
    etaVec1 = zeros(1,length(xiVec));
    etaVec2 = zeros(1,length(xiVec));
    for i = 1:length(xiVec)
        % panel intersection
        etaVec1(i) = P(2) + sqrt((xiVec(i)-P(1))^2 - P(3)^2);
        etaVec2(i) = P(2) - sqrt((xiVec(i)-P(1))^2 - P(3)^2);
        % plane of P
        etaVecP1(i) = P(2) + sqrt((xiVecP(i)-P(1))^2);
        etaVecP2(i) = P(2) - sqrt((xiVecP(i)-P(1))^2);
    end
else
    B = sqrt(Mach^2-1);
    xiMax = P(1) - sqrt(B^2*P(3)^2);
    xiVec = linspace(min,xiMax,n);
    etaVec1 = zeros(1,length(xiVec));
    etaVec2 = zeros(1,length(xiVec));
    for i = 1:length(xiVec)
        etaVec1(i) = P(2) + sqrt((xiVec(i)-P(1))^2/B^2 - P(3)^2);
        etaVec2(i) = P(2) - sqrt((xiVec(i)-P(1))^2/B^2 - P(3)^2);
        
        etaVecP1(i) = P(2) + sqrt((xiVecP(i)-P(1))^2/B^2);
        etaVecP2(i) = P(2) - sqrt((xiVecP(i)-P(1))^2/B^2);
    end
%     pntMin = [P(1) - sqrt(B^2*P(3)^2), P(2), 0];
end

if length(varargin) > 0
    if length(varargin) > 1
        color = varargin{2};
        plot(xiVec,etaVec1,'Color',color,'LineWidth',2)
        plot(xiVec,etaVec2,'Color',color,'LineWidth',2)
    else
        plot(xiVec,etaVec1,'k','LineWidth',2)
        plot(xiVec,etaVec2,'k','LineWidth',2)
    end
else
    plot(xiVec,etaVec1,'k','LineWidth',2)
    plot(xiVec,etaVec2,'k','LineWidth',2)
    plot3(xiVecP,etaVecP1,zVec,'--k','LineWidth',2)
    plot3(xiVecP,etaVecP2,zVec,'--k','LineWidth',2)
end

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
