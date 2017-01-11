function plot_trajectories(metrics)
%Plot the ensemble averages for waveform data
%   INPUT - structure containing ensemble average of some time series data

conditions = fieldnames(metrics.ankle_angle_r_moment);
plotNumber = tight_subplot(2,3,[.05 .03],[.15 .1],[.07 .03]);
x = 0:1:100;
cmap = colormap(parula(250));

conditionLabels = {'TBAS';'cARM1'; 'cARM2'; 'pARM1'; 'pARM2'; 'pARM3'};

% Ankle moments
axes(plotNumber(3))

% Plot of 15 kg slow walking
for cName = 1:2:length(conditions)
	conditionName15 = conditions{cName};
	plotColor = cmap(round(1+20*(cName-1)),:);
	a1 = plot(matfiltfilt(0.01, 6, 4, (metrics.ankle_angle_r_moment.(conditionName15).mean)), 'LineWidth', 1, 'Color',plotColor);
	hold on
	% Plot SD
	[ph,msg]=jbfill(x,matfiltfilt(0.01, 6, 4, (metrics.ankle_angle_r_moment.(conditionName15).mean - metrics.ankle_angle_r_moment.(conditionName15).var))',...
		matfiltfilt(0.01, 6, 4, (metrics.ankle_angle_r_moment.(conditionName15).mean + metrics.ankle_angle_r_moment.(conditionName15).var))', plotColor, plotColor ,0,.1);
	% Make it so the SDs won't appear in legend
	set(get(get(ph,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
t1 = title({'Ankle'; 'dorsiflexion (+)       plantarflexion (-)'}); set(t1, 'FontSize', 16, 'Color', [1 1 1]);
set(gca, 'xlim', [0,100], 'Color', [0 0 0], 'xtick', [], 'xcolor', [0 0 0], 'ycolor', [1 1 1], 'ylim', [-3, 1], 'FontName', 'Calibri', 'FontSize', 16, 'box', 'off');
plot([60,60], [-3,1], '-.', 'Color', [1 1 1], 'LineWidth', 2);
set(gcf, 'Color', [0 0 0]);
l = legend(conditionLabels, 'Location', 'northeast', 'box', 'off', 'orientation', 'horizontal', 'fontsize', 18, 'TextColor', [1 1 1]);
% title(l, 'Armour Conditions');

axes(plotNumber(6))
% Plot of 30 kg slow walking
for cName = 1:6
	tName = cName + 14;
	conditionName30 = conditions{tName};
	plotColor = cmap(round(1+40*(cName-1)),:);
	a2 = plot(matfiltfilt(0.01, 6, 4, (metrics.ankle_angle_r_moment.(conditionName30).mean)), 'LineWidth', 1, 'Color',plotColor);
	hold on
	% Plot SD
	[ph,msg]=jbfill(x,matfiltfilt(0.01, 6, 4, (metrics.ankle_angle_r_moment.(conditionName30).mean - metrics.ankle_angle_r_moment.(conditionName30).var))',...
		matfiltfilt(0.01, 6, 4, (metrics.ankle_angle_r_moment.(conditionName30).mean + metrics.ankle_angle_r_moment.(conditionName30).var))', plotColor, plotColor ,0,.1);
end
x1 = xlabel('Gait Cycle (%)'); set(x1, 'Color', [1 1 1]);
set(gca, 'xlim', [0,100], 'Color', [0 0 0], 'xcolor', [1 1 1], 'ycolor', [1 1 1], 'ylim', [-3, 1], 'FontName', 'Calibri', 'FontSize', 16, 'box', 'off');
plot([62.4,62.4], [-3,1], '-.', 'Color', [1 1 1], 'LineWidth', 2);

% Knee moments
axes(plotNumber(2))
% Plot of 15 kg slow walking
for cName = 1:6
	tName = cName + 2;
	conditionName15 = conditions{tName};
	plotColor = cmap(round(1+40*(cName-1)),:);
	k1 = plot(matfiltfilt(0.01, 6, 4, (metrics.knee_angle_r_moment.(conditionName15).mean)), 'LineWidth', 1, 'Color',plotColor);
	hold on
	% Plot SD
	[ph,msg]=jbfill(x,matfiltfilt(0.01, 6, 4, (metrics.knee_angle_r_moment.(conditionName15).mean - metrics.knee_angle_r_moment.(conditionName15).var))',...
		matfiltfilt(0.01, 6, 4, (metrics.knee_angle_r_moment.(conditionName15).mean + metrics.knee_angle_r_moment.(conditionName15).var))', plotColor, plotColor ,0,.1);
end
t2 = title({'Knee'; 'flexion (+)      extension (-)'}); set(t2, 'FontSize', 16, 'Color', [1 1 1]);
set(gca, 'xlim', [0,100], 'Color', [0 0 0], 'xtick', [], 'xcolor', [0 0 0], 'ycolor', [1 1 1], 'ylim', [-1, 2], 'FontName', 'Calibri', 'FontSize', 16, 'box', 'off');
plot([60,60], [-1,2], '-.', 'Color', [0.9 0.9 0.9], 'LineWidth', 2);

% Plot of 30 kg slow walking
axes(plotNumber(5))
for cName = 1:6
	plotColor = cmap(round(1+40*(cName-1)),:);
	k2 = plot(matfiltfilt(0.01, 6, 4, (knee30_mean(:,cName,1))), 'LineWidth', 1, 'Color',plotColor);
	hold on
	% Plot SD
	[ph,msg]=jbfill(x,matfiltfilt(0.01, 6, 4, knee30Lower(:,cName,1))', matfiltfilt(0.01, 6, 4, knee30Upper(:,cName,1))', plotColor, plotColor ,0,.1);
end
x2 = xlabel('Gait Cycle (%)'); set(x2, 'Color', [1 1 1]);
set(gca, 'xlim', [0,100], 'Color', [0 0 0], 'xcolor', [1 1 1], 'ycolor', [1 1 1], 'ylim', [-1, 2], 'FontName', 'Calibri', 'FontSize', 16, 'box', 'off');
plot([62.4,62.4], [-1,2], '-.', 'Color', [0.9 0.9 0.9], 'LineWidth', 2);

% Hip moments
axes(plotNumber(1))
% Plot of 15 kg slow walking
for cName = 1:6
	plotColor = cmap(round(1+40*(cName-1)),:);
	h1 = plot(matfiltfilt(0.01, 6, 4, (hip15_mean(:,cName,1))), 'LineWidth', 1, 'Color',plotColor);
	hold on
	% Plot SD
	[ph,msg]=jbfill(x,matfiltfilt(0.01, 6, 4, hip15Lower(:,cName,1))', matfiltfilt(0.01, 6, 4, hip15Upper(:,cName,1))', plotColor, plotColor ,0,.1);
end
t3 = title({'Hip'; 'flexion (+)      extension (-)'}); set(t3, 'FontSize', 16, 'Color', [1 1 1]);
y = ylabel({'15 kg'; 'Moment (N.m. kg^-^1)'}); set(y, 'Position', [-8.5, -0.5, 0], 'Color', [1 1 1]);
set(gca, 'xlim', [0,100], 'Color', [0 0 0], 'xtick', [], 'xcolor', [0 0 0], 'ycolor', [1 1 1], 'ylim', [-3, 2], 'FontName', 'Calibri', 'FontSize', 16, 'box', 'off');
plot([60,60], [-3,2], '-.', 'Color', [.9 .9 .9], 'LineWidth', 2);

% Plot of 30 kg slow walking
axes(plotNumber(4))
for cName = 1:6
	plotColor = cmap(round(1+40*(cName-1)),:);
	h2 = plot(matfiltfilt(0.01, 6, 4, (hip30_mean(:,cName,1))), 'LineWidth', 1, 'Color',plotColor);
	hold on
	% Plot SD
	[ph,msg]=jbfill(x,matfiltfilt(0.01, 6, 4, hip30Lower(:,cName,1))', matfiltfilt(0.01, 6, 4, hip30Upper(:,cName,1))', plotColor, plotColor ,0,.1);
end
y = ylabel({'30 kg'; 'Moment (N.m. kg^-^1)'}); set(y, 'Position', [-8.5, -0.5, 0], 'Color', [1 1 1]); x3 = xlabel('Gait Cycle (%)');
set(x3, 'Color', [1 1 1]); set(gca, 'xlim', [0,100], 'Color', [0 0 0], 'xcolor', [1 1 1], 'ycolor', [1 1 1], 'ylim', [-3, 2], 'FontName', 'Calibri', 'FontSize', 16, 'box', 'off');
plot([62.4,62.4], [-3,2], '-.', 'Color', [.9 .9 .9], 'LineWidth', 2);

end

%% Option to save the figures if you want


