function plot_peak_moments(peak_moments, conditionLabels)
%Create bar subplots of the power percentage contribution data
%   INPUT - peak_moments: structure containing the peak moment means
%           and SDs
%         - conditionLabels: cell array containing the names of the conditions

% Function to bring plots together
ha = tight_subplot(2,3,[.05 .03],[.1 .1],[.07 .01]);

% ANKLE
axes(ha(3))
CreatePlotExtraVar(peak_moments.ank15_SD(:,1), peak_moments.ank15_SD(:,2), peak_moments.ank15_mean(:,1), peak_moments.ank15_mean(:,2), conditionLabels, [0,3])
t = title('Ankle Plantarflexion'); set(t, 'Units', 'Normalized', 'Position', [0.5, 1, 0]);
legend({'Slow', 'Fast'}, 'Location', 'northeast', 'box', 'off', 'orientation', 'vertical', 'fontsize', 14);
set(gca, 'YTick', [0, 1, 2, 3]);

axes(ha(6))
CreatePlotExtraVar(peak_moments.ank30_SD(:,1), peak_moments.ank30_SD(:,2), peak_moments.ank30_mean(:,1), peak_moments.ank30_mean(:,2), conditionLabels, [0, 3.5])
set(gca, 'YTick', [0, 1, 2, 3]);

% KNEE
axes(ha(2))
CreatePlotExtraVar(peak_moments.kneeExt15_SD(:,1), peak_moments.kneeExt15_SD(:,2), peak_moments.kneeExt15_mean(:,1)*-1, peak_moments.kneeExt15_mean(:,2)*-1, conditionLabels, [0, 1])
t = title('Knee Extension'); set(t, 'Units', 'Normalized', 'Position', [0.5, 1, 0]);
set(gca, 'YTick', [0, 0.5, 1]);

axes(ha(5))
CreatePlotExtraVar(peak_moments.kneeExt30_SD(:,1), peak_moments.kneeExt30_SD(:,2), peak_moments.kneeExt30_mean(:,1)*-1, peak_moments.kneeExt30_mean(:,2)*-1, conditionLabels, [0, 2])
set(gca, 'xticklabel', conditionLabels, 'fontsize', 14, 'ylim', [0,1], 'YTick', [0, 0.5, 1], 'FontName', 'Calibri', 'box', 'off');

% HIP
axes(ha(1))
CreatePlotExtraVar(peak_moments.hipExt15_SD(:,1), peak_moments.hipExt15_SD(:,2), peak_moments.hipExt15_mean(:,1), peak_moments.hipExt15_mean(:,2), conditionLabels, [0, 4])
y = ylabel({'15 kg'; 'Moment (N.m. kg^-^1)'});  t = title('Hip Extension');
set(y, 'Units', 'Normalized', 'Position', [-0.07, 0.5, 0]); set(t, 'Units', 'Normalized', 'Position', [0.5, 1, 0]);

axes(ha(4))
CreatePlotExtraVar(peak_moments.hipExt30_SD(:,1), peak_moments.hipExt30_SD(:,2), peak_moments.hipExt30_mean(:,1), peak_moments.hipExt30_mean(:,2), conditionLabels, [0, 4])
y = ylabel({'30 kg'; 'Moment (N.m. kg^-^1)'}); set(y, 'Units', 'Normalized', 'Position', [-0.07, 0.5, 0]);

end

