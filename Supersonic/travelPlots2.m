
% Jake Davis
% 01/05/2019

function travelPlots2(vec,phiD,phiS,xAxLab)

fig = figure;
set(fig,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]])
set(gca,'box','off')
hold on

% Doublet Plot
yyaxis left
plot(vec,phiD,'-xb','LineWidth',1,'MarkerSize',7)
xlabel(xAxLab,'Interpreter','latex','FontSize',12)
ylabel('$\phi_{_D}$','Interpreter','latex','FontSize',18)
% xlim([1.5 3.5])
% ylim([-0.7 0.20])
xtickformat('%.1f')
ytickformat('%.1f')
ax1 = gca;
ax1.LineWidth = 1.1;
% ax1.XTick = -0.5:.25:2.0;
% ax1.YTick = -.4:.05:0.3;
ax1.XMinorTick = 'on';
ax1.YMinorTick = 'on';
ax1.TickLength = [.015 .015];
ax1.XGrid = 'on';

% Source Plot
yyaxis right
plot(vec,phiS,'-or','LineWidth',1,'MarkerSize',5)
ylabel('$\phi_{_S}$','Interpreter','latex','FontSize',18)
% ylim([0 0.5])
ytickformat('%.2f')
ax2 = gca;
ax2.LineWidth = 1.1;
% ax2.XTick = 2.5:2.5:22.5;
% ax2.YTick = 0.05:.025:0.20;
ax2.YMinorTick = 'on';
ax2.TickLength = [.015 .015];

lgd = legend('Doublet','Source','Interpreter','latex','Location','east');
lgd.FontSize = 12;
legend boxoff

outerpos = ax1.OuterPosition;
ti = ax1.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax1.Position = [left bottom ax_width ax_height];

hold off

% print(gcf,'figs\5_xPos.png','-dpng','-r300');
% print(gcf,'figs\5_yPos.png','-dpng','-r300');

end
