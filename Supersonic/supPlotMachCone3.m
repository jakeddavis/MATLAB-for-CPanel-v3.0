
% Jake Davis
% 10/06/2018

function supPlotMachCone3(min,max,n,P,Mach,trans,varargin)

both = true;
top = false;
bot = false;
if size(varargin) > 0
    both = false;
    if strcmp(varargin{1},'top') == 1
        top = true;
    elseif strcmp(varargin{1},'bot') == 1
        bot = true;
    else
        both = true;
    end
end

xiVecP = linspace(min,max,n);
% etaVec1 = zeros(1,length(xiVec));
% etaVec2 = zeros(1,length(xiVec));

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

hold on

% Panel-Mach Cone Intersection
if both == true
    if size(varargin) > 0
        plot(xiVec,etaVec1,'-','Color',[.6 .6 .6],'LineWidth',2)
        plot(xiVec,etaVec2,'-','Color',[.6 .6 .6],'LineWidth',2)
        plot(xiVecP,etaVecP1,'--','Color',[.6 .6 .6],'LineWidth',2)
        plot(xiVecP,etaVecP2,'--','Color',[.6 .6 .6],'LineWidth',2)
    else
        plot(xiVec,etaVec1,'k','LineWidth',2)
        plot(xiVec,etaVec2,'k','LineWidth',2)
        plot(xiVecP,etaVecP1,'--k','LineWidth',2)
        plot(xiVecP,etaVecP2,'--k','LineWidth',2)
    end
elseif top == true
    plot(xiVec,etaVec1,'k','LineWidth',2)
    plot(xiVecP,etaVecP1,'--k','LineWidth',2)
else
    plot(xiVec,etaVec2,'k','LineWidth',2)
    plot(xiVecP,etaVecP2,'--k','LineWidth',2)
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
