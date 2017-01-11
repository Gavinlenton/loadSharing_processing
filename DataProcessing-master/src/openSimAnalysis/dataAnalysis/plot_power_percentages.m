function plot_power_percentages(power_percent, conditionLabels)
%Create bar subplots of the power percentage contribution data
%   INPUT - power_percent: structure containing the power percentage means
%           and SDs
%         - conditionLabels: cell arrau containing the names of the conditions

f = figure();
% Use tight_subplot function to bring subplots closer
hz = tight_subplot(2,2,[.03 .03],[.06 .05],[.05 .03]);

xTickLabel = {'Hip'; 'Knee'; 'Ankle'};

% --- POSITIVE POWER --- %
p1 = power_percent.light_slow_mean_pos'; p2 = power_percent.light_fast_mean_pos'; 
p3 = power_percent.heavy_slow_mean_pos'; p4 = power_percent.heavy_fast_mean_pos'; 
p5 = power_percent.light_slow_mean_neg'; p6 = power_percent.light_fast_mean_pos'; 
p7 = power_percent.heavy_slow_mean_pos'; p8 = power_percent.heavy_fast_mean_pos'; 

% 15 kg slow walking
axes(hz(1));
bar(p1(:,2:7)); ylabel('Relative contribution (%)'); t = title('Slow walking speed');
set(t, 'fontsize', 16, 'FontName', 'Calibri');
set(gca, 'xticklabel',[], 'ylim', [0,100], 'fontsize', 14, 'FontName', 'Calibri',...
	'ytick', [0, 25,50,75,100], 'box', 'off');
legend(conditionLabels, 'Location', 'southoutside', 'box', 'off', 'orientation', 'horizontal', 'FontSize', 16);

% 15 kg fast walking
axes(hz(2))
bar(p2(:,2:7)); t = title('Fast walking speed');
set(t, 'fontsize', 16, 'FontName', 'Calibri');
set(gca,'xticklabel',[], 'ylim', [0,100], 'fontsize', 14, 'FontName', 'Calibri',...
	'ytick', [0, 25,50,75,100], 'box', 'off');

% 30 kg slow walking
axes(hz(3))
bar(p3(:,:)); ylabel('Relative contribution (%)');
% set(t, 'Units', 'Normalized', 'Position', [1.03, 1.05, 0], 'fontsize', 16, 'FontName', 'Calibri');
set(gca,'xticklabel',xTickLabel, 'ylim', [0,100], 'fontsize', 14, 'FontName', 'Calibri',...
	'ytick', [0, 25,50,75,100], 'box', 'off');

% 30 kg fast walking
axes(hz(4))
bar(p4(:,:)); %ylabel('Relative contribution (%)'); 
% set(t, 'Units', 'Normalized', 'Position', [1.03, 1.05, 0], 'fontsize', 16, 'FontName', 'Calibri');
set(gca,'xticklabel',xTickLabel, 'ylim', [0,100], 'fontsize', 14, 'FontName', 'Calibri',...
	'ytick', [0, 25,50,75,100], 'box', 'off');

% % --- NEGATIVE POWER --- %
% 
% % 15 kg slow walking
% axes(hz(2))
% bar(power_percent.light_slow_mean_neg(2:7,:));
% set(gca,'xticklabel',conditionLabels, 'xcolor', [0 0 0], 'yticklabel', {}, 'ylim', [0,100], 'fontsize', 14, 'FontName', 'Calibri', 'box', 'off',...
% 	'ytick', [0, 25,50,75,100]);
% 
% % 15 kg fast walking
% axes(hz(4))
% bar(power_percent.light_fast_mean_neg(2:7,:));
% set(gca,'xticklabel',conditionLabels, 'yticklabel', {}, 'ylim', [0,100], 'fontsize', 14, 'FontName', 'Calibri', 'box', 'off',...
% 	'ytick', [0, 25,50,75,100]);
% 
% % 30 kg slow walking
% axes(hz(6))
% bar(power_percent.heavy_slow_mean_neg);
% set(gca,'xticklabel',conditionLabels, 'yticklabel', {}, 'ylim', [0,100], 'fontsize', 14, 'FontName', 'Calibri', 'box', 'off',...
% 	'ytick', [0, 25,50,75,100]);
% 
% % 30 kg fast walking
% axes(hz(8))
% bar(power_percent.heavy_fast_mean_neg);
% set(gca,'xticklabel',conditionLabels, 'yticklabel', {}, 'ylim', [0,100], 'fontsize', 14, 'FontName', 'Calibri', 'box', 'off',...
% 	'ytick', [0, 25,50,75,100]);

applyhatch_pluscolor(f,'xkw.+\/');

end

