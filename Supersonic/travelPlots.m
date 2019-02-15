
% Jake Davis
% 01/05/2019

function travelPlots(vec,phiD,phiS,xAxLab,varargin)

% Doublet Plot
figure(1)
plot(vec,phiD,'-xb','LineWidth',1,'MarkerSize',7)
% Source Plot
figure(2)
plot(vec,phiS,'-xb','LineWidth',1,'MarkerSize',7)

if size(varargin) > 0
    % Doublet Plot
    figure(1)
    hold on
    phiDfar = varargin{1};
    plot(vec,phiDfar,'-or','LineWidth',1,'MarkerSize',5)
    hold off
    
    % Source Plot
    figure(2)
    hold on
    phiSfar = varargin{2};
    plot(vec,phiSfar,'-or','LineWidth',1,'MarkerSize',5)
    hold off
end

% Doublet Plot
figure(1)
hold on
set(gca,'box','off')
xlabel(xAxLab,'Interpreter','latex','FontSize',12)
ylabel('$\phi_{_D}$','Interpreter','latex','FontSize',18)
% xlim([2.5 22.5])
% ylim([0 0.5])
xtickformat('%.1f')
ytickformat('%.2f')
ax = gca;
ax.LineWidth = 1.1;
% ax1.XTick = 2.5:2.5:22.5;
% ax1.YTick = 0.05:.025:0.20;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.TickLength = [.015 .015];
ax.XAxisLocation = 'origin';
grid on
% ax1.XGrid = 'on';
if size(varargin) > 0
    lgd = legend('Linear','Point',...
        'Interpreter','latex','Location','northeast');
    lgd.FontSize = 12;
    legend boxoff
end
x0=1;
y0=6;
width=2.9;
height=2.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
hold off

print(gcf,'5_Dub.png','-dpng','-r300');

% Source Plot
figure(2)
hold on
set(gca,'box','off')
xlabel(xAxLab,'Interpreter','latex','FontSize',12)
ylabel('$\phi_{_S}$','Interpreter','latex','FontSize',18)
% xlim([2.5 22.5])
% ylim([0 0.5])
xtickformat('%.1f')
ytickformat('%.2f')
ax = gca;
ax.LineWidth = 1.1;
% ax1.XTick = 2.5:2.5:22.5;
% ax1.YTick = 0.05:.025:0.20;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.TickLength = [.015 .015];
ax.XAxisLocation = 'origin';
grid on
% ax1.XGrid = 'on';
if size(varargin) > 0
    lgd = legend('Constant','Point',...
        'Interpreter','latex','Location','northeast');
    lgd.FontSize = 12;
    legend boxoff
end
x0=1;
y0=1;
width=2.9;
height=2.5;
set(gcf,'units','inches','position',[x0,y0,width,height])
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
hold off

print(gcf,'5_Src.png','-dpng','-r300');

end
