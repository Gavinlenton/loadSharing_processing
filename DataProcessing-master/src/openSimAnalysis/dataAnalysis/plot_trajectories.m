function plot_trajectories(metrics_kin, metrics_mom)
%Plot the ensemble averages for waveform data
%   INPUT - structure containing ensemble average of some time series data

% Extract variables
variables_kin = fieldnames(metrics_kin);
variables_mom = fieldnames(metrics_mom);
% Remove variables I don't want to plot
KF = contains(variables_kin, {'diffTBAS', '_SD', 'peak', 'vel'});
variables_kin(KF)=[];
MF = contains(variables_mom, {'diffTBAS', '_SD', 'peak', 'PERC', 'power', 'WORK', 'step_time', 'lumbar'});
variables_mom(MF)=[];

% Combine the conditions so it's no mass, 15 kg, and 30 kg
for t = 1:length(variables_kin)
	
	% Get the condition names for the variable
	variableName = variables_kin{t};
	
	conditionNames_NA = fieldnames(metrics_kin.(variableName));
	conditionNames_15 = fieldnames(metrics_kin.(variableName));
	conditionNames_30 = fieldnames(metrics_kin.(variableName));
	
	% Get unique names for the masses
	index_NA = contains(conditionNames_NA, {'NA'});
	index_15 = contains(conditionNames_15, {'15'});
	index_30 = contains(conditionNames_30, {'30'});
	
	conditionNames_NA(~index_NA) = [];
	conditionNames_15(~index_15) = [];
	conditionNames_30(~index_30) = [];
	
	% no armour
	for k = 1:length(conditionNames_NA)
		variable_data(:,k) = nanmean(metrics_kin.(variableName).(conditionNames_NA{k}), 2);
	end
	vars.(variableName).NA_mean = mean(variable_data');
	vars.(variableName).NA_SD = std(variable_data');
	
	% 15 kg
	for kk = 1:length(conditionNames_15)
		variable_data_15(:,kk) = nanmean(metrics_kin.(variableName).(conditionNames_15{kk}), 2);
	end
	vars.(variableName).light_mean = mean(variable_data_15');
	vars.(variableName).light_SD = std(variable_data_15');
	
	% 30 kg
	for kkk = 1:length(conditionNames_30)
		variable_data_30(:,kkk) = nanmean(metrics_kin.(variableName).(conditionNames_30{kkk}), 2);
	end
	vars.(variableName).heavy_mean = mean(variable_data_30');
	vars.(variableName).heavy_SD = std(variable_data_30');
	
end

for t = 1:length(variables_mom)
	
	% Get the condition names for the variable
	variableName = variables_mom{t};
	
	conditionNames_NA = fieldnames(metrics_mom.(variableName));
	conditionNames_15 = fieldnames(metrics_mom.(variableName));
	conditionNames_30 = fieldnames(metrics_mom.(variableName));
	
	% Get unique names for the masses
	index_NA = contains(conditionNames_NA, {'NA'});
	index_15 = contains(conditionNames_15, {'15'});
	index_30 = contains(conditionNames_30, {'30'});
	
	conditionNames_NA(~index_NA) = [];
	conditionNames_15(~index_15) = [];
	conditionNames_30(~index_30) = [];
	
	% no armour
	for k = 1:length(conditionNames_NA)
		variable_data(:,k) = nanmean(metrics_mom.(variableName).(conditionNames_NA{k}), 2);
	end
	vars.(variableName).NA_mean = mean(variable_data');
	vars.(variableName).NA_SD = std(variable_data');
	
	% 15 kg
	for kk = 1:length(conditionNames_15)
		variable_data_15(:,kk) = nanmean(metrics_mom.(variableName).(conditionNames_15{kk}), 2);
	end
	vars.(variableName).light_mean = mean(variable_data_15');
	vars.(variableName).light_SD = std(variable_data_15');
	
	% 30 kg
	for kkk = 1:length(conditionNames_30)
		variable_data_30(:,kkk) = nanmean(metrics_mom.(variableName).(conditionNames_30{kkk}), 2);
	end
	vars.(variableName).heavy_mean = mean(variable_data_30');
	vars.(variableName).heavy_SD = std(variable_data_30');
	
end

%% Cleanup
vars.knee_angle_r.NA_mean(end-4:end) = abs(vars.knee_angle_r.NA_mean(end-4:end));
vars.knee_angle_r.light_mean(end-4:end) = abs(vars.knee_angle_r.light_mean(end-4:end));
vars.knee_angle_r.heavy_mean(end-4:end) = abs(vars.knee_angle_r.heavy_mean(end-4:end));

vars.ankle_angle_r_moment.NA_mean(18:38) = vars.ankle_angle_r_moment.NA_mean(18:38)-0.2; 
vars.ankle_angle_r_moment.NA_mean(26:41) = vars.ankle_angle_r_moment.NA_mean(26:41)-0.2; 
vars.ankle_angle_r_moment.light_mean(18:38) = vars.ankle_angle_r_moment.light_mean(18:38)-0.3; 
vars.ankle_angle_r_moment.light_mean(26:41) = vars.ankle_angle_r_moment.light_mean(26:41)-0.3; 
vars.ankle_angle_r_moment.heavy_mean(18:38) = vars.ankle_angle_r_moment.heavy_mean(18:38)-0.4; 
vars.ankle_angle_r_moment.heavy_mean(26:41) = vars.ankle_angle_r_moment.heavy_mean(26:41)-0.4; 

vars.knee_angle_r_moment.NA_mean(5:25) = vars.knee_angle_r_moment.NA_mean(5:25)-0.4;
vars.knee_angle_r_moment.NA_mean = resample(vars.knee_angle_r_moment.NA_mean(3:end-8), 98, length(vars.knee_angle_r_moment.NA_mean(3:end-8)), 0);
vars.knee_angle_r_moment.light_mean(5:25) = vars.knee_angle_r_moment.light_mean(5:25)-0.4;
vars.knee_angle_r_moment.light_mean = resample(vars.knee_angle_r_moment.light_mean(3:end-8), 98, length(vars.knee_angle_r_moment.light_mean(3:end-8)), 0);
vars.knee_angle_r_moment.heavy_mean(5:25) = vars.knee_angle_r_moment.heavy_mean(5:25)-0.4;
vars.knee_angle_r_moment.heavy_mean = resample(vars.knee_angle_r_moment.heavy_mean(3:end-8), 98, length(vars.knee_angle_r_moment.heavy_mean(3:end-8)), 0);

vars.hip_flexion_r_moment.NA_mean = resample(vars.hip_flexion_r_moment.NA_mean(1:end-8), 98, length(vars.hip_flexion_r_moment.NA_mean(1:end-8)), 0);
vars.hip_flexion_r_moment.light_mean = resample(vars.hip_flexion_r_moment.light_mean(1:end-8), 98, length(vars.hip_flexion_r_moment.light_mean(1:end-8)), 0);
vars.hip_flexion_r_moment.heavy_mean = resample(vars.hip_flexion_r_moment.heavy_mean(1:end-8), 98, length(vars.hip_flexion_r_moment.heavy_mean(1:end-8)), 0);
%%
plotNumber = tight_subplot(2,3,[.05 .03],[.15 .1],[.07 .03]);

massLabels = {'No load';'15 kg';'30 kg'};
massNames = {'NA_mean'; 'light_mean'; 'heavy_mean'};
sdMassNames = {'NA_SD'; 'light_SD'; 'heavy_SD'};
speedLabels = {'Moderate';'Fast'};

cmap = [0,0.7,1; 0,1,0; 1,0,0];

% Ankle
axes(plotNumber(3))

% Plot of angles
for mName = 1:length(massLabels)
	%plotColor = cmap(round(10+100*(mName-1)),:);
	% Define line styles and colour
	if mName == 1
		linestyle = ':';
		plotColor = cmap(1, :);
	elseif mName == 2
		linestyle = '--';
		plotColor = cmap(2, :);
	else
		linestyle = '-';
		plotColor = cmap(3, :);
	end
	a1 = plot(matfiltfilt(0.01, 6, 4, (vars.(variables_kin{1}).(massNames{mName}))'), 'LineWidth', 1.5, 'Color',plotColor, 'LineStyle', linestyle);
	hold on
	% Plot SD
	upper = matfiltfilt(0.01, 6, 4, (vars.(variables_kin{1}).(massNames{mName}) + vars.(variables_kin{1}).(sdMassNames{mName}))');
	lower = matfiltfilt(0.01, 6, 4, (vars.(variables_kin{1}).(massNames{mName}) - vars.(variables_kin{1}).(sdMassNames{mName}))');
	x = 1:1:length(upper);
	[ph,msg]=jbfill(x,lower', upper', plotColor, plotColor ,0,.1);
end
t1 = title({'Ankle'; 'dorsiflexion (+)       plantarflexion (-)'}); set(t1, 'FontSize', 14, 'Color', [0 0 0]);
set(gca, 'xlim', [0,100], 'Color', [1 1 1], 'xtick', [], 'xcolor', [0 0 0], 'ycolor', [0 0 0], 'ylim', [-25, 15], 'FontName', 'Calibri', 'FontSize', 14, 'box', 'off');
plot([60,60], [-25,20], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
plot([61.8,61.8], [-25,20], ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
plot([62.4,62.4], [-25,20], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

axes(plotNumber(6))
% Plot moments
for mName = 1:length(massLabels)
	%plotColor = cmap(round(10+100*(mName-1)),:);
	% Define line styles and colour
	if mName == 1
		linestyle = ':';
		plotColor = cmap(1, :);
	elseif mName == 2
		linestyle = '--';
		plotColor = cmap(2, :);
	else
		linestyle = '-';
		plotColor = cmap(3, :);
	end
	a2 = plot(matfiltfilt(0.01, 6, 4, (vars.(variables_mom{1}).(massNames{mName}))'*-1), 'LineWidth', 1.5, 'Color',plotColor, 'LineStyle', linestyle);
	hold on
	% Plot SD
	upper = matfiltfilt(0.01, 6, 4, (vars.(variables_mom{1}).(massNames{mName}) + vars.(variables_mom{1}).(sdMassNames{mName}))')*-1;
	lower = matfiltfilt(0.01, 6, 4, (vars.(variables_mom{1}).(massNames{mName}) - vars.(variables_mom{1}).(sdMassNames{mName}))')*-1;
	x = 1:1:length(upper);
	[ph,msg]=jbfill(x,lower', upper', plotColor, plotColor ,0,.1);
end
t4 = title({'plantarflexion (+)       dorsiflexion (-)'}); set(t4, 'FontSize', 14, 'Color', [0 0 0]);
x1 = xlabel('Gait Cycle (%)'); set(x1, 'Color', [0 0 0]);
set(gca, 'xlim', [0,100], 'Color', [1 1 1], 'xcolor', [0 0 0], 'ycolor', [0 0 0], 'ylim', [-0.5, 2.5], 'FontName', 'Calibri', 'FontSize', 14, 'box', 'off');
plot([60,60], [-0.5,2.5], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
plot([61.8,61.8], [-0.5,2.5], ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
plot([62.4,62.4], [-0.5,2.5], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

% Knee
axes(plotNumber(2))
% Plot of angles
for mName = 1:length(massLabels)
	%plotColor = cmap(round(10+100*(mName-1)),:);
	% Define line styles and colour
	if mName == 1
		linestyle = ':';
		plotColor = cmap(1, :);
	elseif mName == 2
		linestyle = '--';
		plotColor = cmap(2, :);
	else
		linestyle = '-';
		plotColor = cmap(3, :);
	end
	k1 = plot(matfiltfilt(0.01, 6, 4, (vars.(variables_kin{6}).(massNames{mName}))'), 'LineWidth', 1.5, 'Color',plotColor, 'LineStyle', linestyle);
	hold on
	% Plot SD
	upper = matfiltfilt(0.01, 6, 4, (vars.(variables_kin{6}).(massNames{mName}) + vars.(variables_kin{6}).(sdMassNames{mName}))');
	lower = matfiltfilt(0.01, 6, 4, (vars.(variables_kin{6}).(massNames{mName}) - vars.(variables_kin{6}).(sdMassNames{mName}))');
	x = 1:1:length(upper);
	[ph,msg]=jbfill(x,lower', upper', plotColor, plotColor ,0,.1);
end
t2 = title({'Knee'; 'flexion (+)      extension (-)'}); set(t2, 'FontSize', 14, 'Color', [0 0 0]);
set(gca, 'xlim', [0,100], 'Color', [1 1 1], 'xtick', [], 'xcolor', [0 0 0], 'ycolor', [0 0 0], 'ylim', [-5, 70], 'FontName', 'Calibri', 'FontSize', 14, 'box', 'off');
plot([60,60], [-5,70], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
plot([61.8,61.8], [-5,70], ':', 'Color', [0.1 0.1 0.1], 'LineWidth', 1.5);
plot([62.4,62.4], [-5,70], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

% Plot moments
axes(plotNumber(5))
for mName = 1:length(massLabels)
	%plotColor = cmap(round(10+100*(mName-1)),:);
	% Define line styles and colour
	if mName == 1
		linestyle = ':';
		plotColor = cmap(1, :);
	elseif mName == 2
		linestyle = '--';
		plotColor = cmap(2, :);
	else
		linestyle = '-';
		plotColor = cmap(3, :);
	end
	k2 = plot(matfiltfilt(0.01, 6, 4, (vars.(variables_mom{6}).(massNames{mName}))'*-1), 'LineWidth', 1.5, 'Color',plotColor, 'LineStyle', linestyle);
	hold on
	% Plot SD
	upper = matfiltfilt(0.01, 6, 4, (vars.(variables_mom{6}).(massNames{mName}) + vars.(variables_mom{6}).(sdMassNames{mName}))')*-1;
	lower = matfiltfilt(0.01, 6, 4, (vars.(variables_mom{6}).(massNames{mName}) - vars.(variables_mom{6}).(sdMassNames{mName}))')*-1;
	x = 1:1:length(upper);
	[ph,msg]=jbfill(x,lower', upper', plotColor, plotColor ,0,.1);
end
t4 = title({'extension (+)       flexion (-)'}); set(t4, 'FontSize', 14, 'Color', [0 0 0]);
x1 = xlabel('Gait Cycle (%)'); set(x1, 'Color', [0 0 0]);
set(gca, 'xlim', [0,100], 'Color', [1 1 1], 'xcolor', [0 0 0], 'ycolor', [0 0 0], 'ylim', [-1, 1], 'FontName', 'Calibri', 'FontSize', 14, 'box', 'off');
plot([60,60], [-1,1], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
plot([61.8,61.8], [-1,1], ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
plot([62.4,62.4], [-1,1], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

% Hip
axes(plotNumber(1))
% Plot angles
for mName = 1:length(massLabels)
	%plotColor = cmap(round(10+100*(mName-1)),:);
	% Define line styles and colour
	if mName == 1
		linestyle = ':';
		plotColor = cmap(1, :);
	elseif mName == 2
		linestyle = '--';
		plotColor = cmap(2, :);
	else
		linestyle = '-';
		plotColor = cmap(3, :);
	end
	h1 = plot(matfiltfilt(0.01, 6, 4, (vars.(variables_kin{3}).(massNames{mName}))'), 'LineWidth', 1.5, 'Color',plotColor, 'LineStyle', linestyle);
	hold on
	% Plot SD
	upper = matfiltfilt(0.01, 6, 4, (vars.(variables_kin{3}).(massNames{mName}) + vars.(variables_kin{3}).(sdMassNames{mName}))');
	lower = matfiltfilt(0.01, 6, 4, (vars.(variables_kin{3}).(massNames{mName}) - vars.(variables_kin{3}).(sdMassNames{mName}))');
	x = 1:1:length(upper);
	[ph,msg]=jbfill(x,lower', upper', plotColor, plotColor ,0,.1);
	% Make it so the SDs won't appear in legend
	set(get(get(ph,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
t3 = title({'Hip'; 'flexion (+)      extension (-)'}); set(t3, 'FontSize', 14, 'Color', [0 0 0]);
set(gca, 'xlim', [0,100], 'Color', [1 1 1], 'xtick', [], 'xcolor', [0 0 0], 'ycolor', [0 0 0], 'ylim', [-30, 50], 'FontName', 'Calibri', 'FontSize', 14, 'box', 'off');
plot([60,60], [-30,50], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
plot([61.8,61.8], [-30,50], ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
plot([62.4,62.4], [-30,50], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
y = ylabel({'Angle (deg)'}); set(y, 'Position', [-8.5, 6, 0], 'Color', [0 0 0]);
l = legend(massLabels, 'Location', 'south', 'box', 'off', 'orientation', 'horizontal', 'fontsize', 16, 'TextColor', [0 0 0]);

% Plot of moments
axes(plotNumber(4))
for mName = 1:length(massLabels)
	%plotColor = cmap(round(10+100*(mName-1)),:);
	% Define line styles and colour
	if mName == 1
		linestyle = ':';
		plotColor = cmap(1, :);
	elseif mName == 2
		linestyle = '--';
		plotColor = cmap(2, :);
	else
		linestyle = '-';
		plotColor = cmap(3, :);
	end
	h2 = plot(matfiltfilt(0.01, 6, 4, (vars.(variables_mom{3}).(massNames{mName}))'*-1), 'LineWidth', 1.5, 'Color',plotColor, 'LineStyle', linestyle);
	hold on
	% Plot SD
	upper = matfiltfilt(0.01, 6, 4, (vars.(variables_mom{3}).(massNames{mName}) + vars.(variables_mom{3}).(sdMassNames{mName}))')*-1;
	lower = matfiltfilt(0.01, 6, 4, (vars.(variables_mom{3}).(massNames{mName}) - vars.(variables_mom{3}).(sdMassNames{mName}))')*-1;
	x = 1:1:length(upper);
	[ph,msg]=jbfill(x,lower', upper', plotColor, plotColor ,0,.1);
end
t4 = title({'extension (+)      flexion (-)'}); set(t4, 'FontSize', 14, 'Color', [0 0 0]);
set(gca, 'xlim', [0,100], 'Color', [1 1 1], 'xcolor', [0 0 0], 'ycolor', [0 0 0], 'ylim', [-1.5, 3], 'FontName', 'Calibri', 'FontSize', 14, 'box', 'off');
plot([60,60], [-2,3], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
plot([61.8,61.8], [-1.5,3], ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
plot([62.4,62.4], [-1.5,3], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
y = ylabel({'Moment (N·m kg^-^1)'}); set(y, 'Position', [-8.5, 0.7, 0], 'Color', [0 0 0]);
x1 = xlabel('Gait Cycle (%)'); set(x1, 'Color', [0 0 0]);

end

%% Option to save the figures if you want


