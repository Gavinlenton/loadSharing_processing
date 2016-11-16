function CreatePlotExtraVar(sd15, sd30, mean15, mean30, conditionLabels, limits)
%Create bar plot with error bars of pad data
%   Detailed explanation goes here
SD = [sd15, sd30];
MeanValue = [mean15, mean30];

errY = zeros(size(SD,1), size(SD,2), 2);
errY(:,:,1) = 0.*SD;   % 0% lower error
errY(:,:,2) = 1.*SD;   % normal upper error

h = barwitherr(errY,MeanValue);

% % Uncomment for blue and red color
% set(h(1),'BarWidth',0.7, 'FaceColor', [0 0.447058826684952 0.74117648601532]);
% set(h(2),'BarWidth',0.7, 'FaceColor', [0.8 0.2 0.4]); % The bars will nOT TOUCH EACHOTHER

% Black and white color
set(h(1),'BarWidth',0.8, 'FaceColor', [0.1, 0.6, 0.8]);
set(h(2),'BarWidth',0.8, 'FaceColor', [0.9, 0.85, 0]); % The bars will nOT TOUCH EACHOTHER

set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 14)
set(gca,'YGrid','off')
set(gca,'box', 'off')
ylim(limits)
set(gca,'XGrid','off', 'xticklabel', conditionLabels)
set(gca,'GridLineStyle','none')

% Insert horizontal line at 100%
hold on
g = plot([0 6], [100,100]);
set(g, 'LineWidth', 2.0)
set(g, 'Color', [0 0 0], 'LineStyle', '--')
end
