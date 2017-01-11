function plot_trajectories_diffTBAS(metrics)
%Plot the ensemble averages for waveform data
%   INPUT - structure containing ensemble average of some time series data

metricsNames = fieldnames(metrics);
conditions = fieldnames(metrics.ankle_angle_r_moment);
plotNumber = tight_subplot(2,3,[.05 .03],[.15 .1],[.07 .03]);
x = 0:1:100;
cmap = colormap(parula(250));

conditionLabels = {'cARM1'; 'cARM2'; 'pARM1'; 'pARM2'; 'pARM3'};

% Ankle moments
axes(plotNumber(3))
k = 1;
% Plot of 15 kg slow walking
for cName = 1:2:length(conditions)/2
	conditionName = conditions{cName};
	plotColor = cmap(round(1+40*(k-1)),:);
	a1 = plot(matfiltfilt(0.01, 6, 4, (metrics.(metricsNames{7}).(conditionName).mean)), 'LineWidth', 1, 'Color',plotColor);
	hold on
	% Plot SD
	[ph,msg]=jbfill(x,matfiltfilt(0.01, 6, 4, (metrics.(metricsNames{7}).(conditionName).mean - metrics.(metricsNames{7}).(conditionName).var))',...
		matfiltfilt(0.01, 6, 4, (metrics.(metricsNames{7}).(conditionName).mean + metrics.(metricsNames{7}).(conditionName).var))', plotColor, plotColor ,0,.1);
	% Make it so the SDs won't appear in legend
	set(get(get(ph,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
	k = k+1;
end
t1 = title('Ankle'); set(t1, 'FontSize', 16, 'Color', [1 1 1]);
set(gca, 'xlim', [0,100], 'Color', [0 0 0], 'xtick', [], 'xcolor', [0 0 0], 'ycolor', [1 1 1], 'ylim', [-0.3, 0.3], 'FontName', 'Calibri', 'FontSize', 16, 'box', 'off');
plot([60,60], [-0.3,0.3], '-.', 'Color', [1 1 1], 'LineWidth', 2);
set(gcf, 'Color', [0 0 0]);
l = legend(conditionLabels, 'Location', 'northeast', 'box', 'off', 'orientation', 'horizontal', 'fontsize', 18, 'TextColor', [1 1 1]);
% title(l, 'Armour Conditions');

axes(plotNumber(6))
k = 1;
% Plot of 30 kg slow walking
for cName = 1:2:length(conditions)/2
	conditionName = conditions{cName+10};
	plotColor = cmap(round(1+40*(k-1)),:);
	a2 = plot(matfiltfilt(0.01, 6, 4, (metrics.(metricsNames{7}).(conditionName).mean)), 'LineWidth', 1, 'Color',plotColor);
	hold on
	% Plot SD
	[ph,msg]=jbfill(x,matfiltfilt(0.01, 6, 4, (metrics.(metricsNames{7}).(conditionName).mean - metrics.(metricsNames{7}).(conditionName).var))',...
		matfiltfilt(0.01, 6, 4, (metrics.(metricsNames{7}).(conditionName).mean + metrics.(metricsNames{7}).(conditionName).var))', plotColor, plotColor ,0,.1);
	k = k+1;
end
x1 = xlabel('Gait Cycle (%)'); set(x1, 'Color', [1 1 1]);
set(gca, 'xlim', [0,100], 'Color', [0 0 0], 'xcolor', [1 1 1], 'ycolor', [1 1 1], 'ylim', [-0.5, 0.5], 'FontName', 'Calibri', 'FontSize', 16, 'box', 'off');
plot([62.4,62.4], [-0.5,0.5], '-.', 'Color', [1 1 1], 'LineWidth', 2);

% Knee moments
axes(plotNumber(2))
k = 1;
% Plot of 15 kg slow walking
for cName = 1:2:length(conditions)/2
	conditionName = conditions{cName};
	plotColor = cmap(round(1+40*(k-1)),:);
	k1 = plot(matfiltfilt(0.01, 6, 4, (metrics.(metricsNames{4}).(conditionName).mean)), 'LineWidth', 1, 'Color',plotColor);
	hold on
	% Plot SD
	[ph,msg]=jbfill(x,matfiltfilt(0.01, 6, 4, (metrics.(metricsNames{4}).(conditionName).mean - metrics.(metricsNames{4}).(conditionName).var))',...
		matfiltfilt(0.01, 6, 4, (metrics.(metricsNames{4}).(conditionName).mean + metrics.(metricsNames{4}).(conditionName).var))', plotColor, plotColor ,0,.1);
	k = k+1;
end
t2 = title('Knee'); set(t2, 'FontSize', 16, 'Color', [1 1 1]);
set(gca, 'xlim', [0,100], 'Color', [0 0 0], 'xtick', [], 'xcolor', [0 0 0], 'ycolor', [1 1 1], 'ylim', [-0.3, 0.3], 'FontName', 'Calibri', 'FontSize', 16, 'box', 'off');
plot([60,60], [-0.3,0.3], '-.', 'Color', [0.9 0.9 0.9], 'LineWidth', 2);

% Plot of 30 kg slow walking
axes(plotNumber(5))
k = 1;
for cName = 1:2:length(conditions)/2
	conditionName = conditions{cName+10};
	plotColor = cmap(round(1+40*(k-1)),:);
	k2 = plot(matfiltfilt(0.01, 6, 4, (metrics.(metricsNames{4}).(conditionName).mean)), 'LineWidth', 1, 'Color',plotColor);
	hold on
	% Plot SD
	[ph,msg]=jbfill(x,matfiltfilt(0.01, 6, 4, (metrics.(metricsNames{4}).(conditionName).mean - metrics.(metricsNames{4}).(conditionName).var))',...
		matfiltfilt(0.01, 6, 4, (metrics.(metricsNames{4}).(conditionName).mean + metrics.(metricsNames{4}).(conditionName).var))', plotColor, plotColor ,0,.1);
	k = k+1;
end
x2 = xlabel('Gait Cycle (%)'); set(x2, 'Color', [1 1 1]);
set(gca, 'xlim', [0,100], 'Color', [0 0 0], 'xcolor', [1 1 1], 'ycolor', [1 1 1], 'ylim', [-0.5, 0.5], 'FontName', 'Calibri', 'FontSize', 16, 'box', 'off');
plot([62.4,62.4], [-0.5,0.5], '-.', 'Color', [0.9 0.9 0.9], 'LineWidth', 2);

% Hip moments
axes(plotNumber(1))
k = 1;
% Plot of 15 kg slow walking
for cName = 1:2:length(conditions)/2
	conditionName = conditions{cName};
	plotColor = cmap(round(1+40*(k-1)),:);
	h1 = plot(matfiltfilt(0.01, 6, 4, (metrics.(metricsNames{1}).(conditionName).mean)), 'LineWidth', 1, 'Color',plotColor);
	hold on
	% Plot SD
	[ph,msg]=jbfill(x,matfiltfilt(0.01, 6, 4, (metrics.(metricsNames{1}).(conditionName).mean - metrics.(metricsNames{1}).(conditionName).var))',...
		matfiltfilt(0.01, 6, 4, (metrics.(metricsNames{1}).(conditionName).mean + metrics.(metricsNames{1}).(conditionName).var))', plotColor, plotColor ,0,.1);
	k = k+1;
end
t3 = title('Hip'); set(t3, 'FontSize', 16, 'Color', [1 1 1]);
y = ylabel({'15 kg'; 'Moment (N.m. kg^-^1)'}); set(y, 'Position', [-8.5, 0, 0], 'Color', [1 1 1]);
set(gca, 'xlim', [0,100], 'Color', [0 0 0], 'xtick', [], 'xcolor', [0 0 0], 'ycolor', [1 1 1], 'ylim', [-0.5, 0.5], 'FontName', 'Calibri', 'FontSize', 16, 'box', 'off');
plot([60,60], [-0.5,0.5], '-.', 'Color', [.9 .9 .9], 'LineWidth', 2);

% Plot of 30 kg slow walking
axes(plotNumber(4))
k = 1;
for cName = 1:2:length(conditions)/2
	conditionName = conditions{cName+10};
	plotColor = cmap(round(1+40*(k-1)),:);
	h2 = plot(matfiltfilt(0.01, 6, 4, (metrics.(metricsNames{1}).(conditionName).mean)), 'LineWidth', 1, 'Color',plotColor);
	hold on
	% Plot SD
	[ph,msg]=jbfill(x,matfiltfilt(0.01, 6, 4, (metrics.(metricsNames{1}).(conditionName).mean - metrics.(metricsNames{1}).(conditionName).var))',...
		matfiltfilt(0.01, 6, 4, (metrics.(metricsNames{1}).(conditionName).mean + metrics.(metricsNames{1}).(conditionName).var))', plotColor, plotColor ,0,.1);
	k = k+1;
end
y = ylabel({'30 kg'; 'Moment (N.m. kg^-^1)'}); set(y, 'Position', [-8.5, 0, 0.5], 'Color', [1 1 1]); x3 = xlabel('Gait Cycle (%)');
set(x3, 'Color', [1 1 1]); set(gca, 'xlim', [0,100], 'Color', [0 0 0], 'xcolor', [1 1 1], 'ycolor', [1 1 1], 'ylim', [-0.5, 0.5], 'FontName', 'Calibri', 'FontSize', 16, 'box', 'off');
plot([62.4,62.4], [-0.5,0.5], '-.', 'Color', [.9 .9 .9], 'LineWidth', 2);

end

%% Option to save the figures if you want


