%% Script to Process data from Load sharing trials

% IK, ID, muscle analysis, and SO taken from Batch OpenSim Processing Scripts (BOPS)
% Copyright (C) 2015 Alice Mantoan, Monica Reggiani
% <https://simtk.org/home/bops>
% Please acknowledge the authors of BOPS when using this script.
clc; clear;

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath(genpath('GeneratePlots'));

%% Load the data
BasePath=uigetdir('../../../', 'Select Elaborated Data Folder');
cd(BasePath);
load('IDMetrics_all.mat');

% Sort IDMetrics so it's in correct order
subjects = fieldnames(ID_metrics);
subjectNumber = zeros(length(subjects), 1);
newOrder = {'NA_slow', 'NA_fast', 'TBAS15_slow', 'TBAS15_fast', 'CRYE15_slow', 'CRYE15_fast',...
	'TYR15_slow', 'TYR15_fast', 'USMC15_slow', 'USMC15_fast', 'CORE15_slow', 'CORE15_fast',...
	'SORD15_slow', 'SORD15_fast', 'TBAS30_slow', 'TBAS30_fast', 'CRYE30_slow', 'CRYE30_fast',...
	'TYR30_slow', 'TYR30_fast', 'USMC30_slow', 'USMC30_fast', 'CORE30_slow', 'CORE30_fast',...
	'SORD30_slow', 'SORD30_fast'};

for t = 1:length(subjects)
	
	% Determine conditions for that subject
	subjectFields = fieldnames(ID_metrics.(subjects{t}));
	subjectName = subjects{t};
	Key   = '_';
	Index = strfind(subjectName, Key);
	Value = sscanf(subjectName(Index(1) + length(Key):end), '%g', 1);
	subjectNumber(t, 1) = Value;
	
	if length(subjectFields) ~= 26
		% Find which conditions are missing
		diffFields = setdiff(newOrder, subjectFields);
		
		% Remove those from new order variable
		for k = 1:length(diffFields)
			if k == 1
			subject_new_order = newOrder(cellfun('isempty', strfind(newOrder, diffFields(k))));
			else
				subject_new_order = subject_new_order(cellfun('isempty', strfind(subject_new_order, diffFields(k))));
			end
		end
		
		% Reorder conditions
		ID_metrics_ordered.(subjects{t}) = orderfields(ID_metrics.(subjects{t}), subject_new_order);
	else
		ID_metrics_ordered.(subjects{t}) = orderfields(ID_metrics.(subjects{t}), newOrder);
	end
end

% Determine how many condition fields exist
conditions = fieldnames(ID_metrics_ordered.(subjects{1}))';

% Create structure to store outputs
all_variables = fieldnames(ID_metrics_ordered.(subjects{1}).(conditions{1}));
all_variables(ismember(all_variables,{'step_time', 'POWER_POS_TOTAL', 'POWER_NEG_TOTAL', 'lumbar_extension_moment',...
	'lumbar_bending_moment', 'lumbar_rotation_moment'}))=[];
metrics_for_plot = struct();

% Loop through all conditions
for cName = 1:length(conditions)
	
	conditionName = conditions{cName};
	
	% Loop through subjects
	for sName = 1:length(subjects)
		
		% Determine if subject has that condition name
		if any(strcmp(fieldnames(ID_metrics_ordered.(subjects{sName})), conditionName))
			
			% Loop through the variables of interest
			for k = 1:length(all_variables)
				variableName = all_variables{k};
				
				switch variableName
					
					% Integer variables are stored as subjects (rows) * conditions
					% (columns) and vectors stored as variable * subject (columns) *
					% conditions (z direction)
					case 'hip_flexion_r_moment'
						metrics_for_plot.(variableName)(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						metrics_for_plot.hip_flexion_p(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
						metrics_for_plot.hip_ext_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						metrics_for_plot.hip_flex_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
					case 'hip_adduction_r_moment'
						metrics_for_plot.(variableName)(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						metrics_for_plot.hip_flexion_p(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
						metrics_for_plot.hip_add_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						metrics_for_plot.hip_abd_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
					case 'hip_rotation_r_moment'
						metrics_for_plot.(variableName)(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						metrics_for_plot.hip_flexion_p(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
						metrics_for_plot.hip_int_rot_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						metrics_for_plot.hip_ext_rot_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
					case 'knee_angle_r_moment'
						metrics_for_plot.(variableName)(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						metrics_for_plot.knee_flexion_p(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
						metrics_for_plot.knee_ext_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						metrics_for_plot.knee_flex_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
					case 'knee_adduction_r_moment'
						metrics_for_plot.(variableName)(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						metrics_for_plot.knee_add_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						
					case 'knee_rotation_r_moment'
						metrics_for_plot.(variableName)(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						metrics_for_plot.knee_int_rot_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						metrics_for_plot.knee_ext_rot_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
					case 'ankle_angle_r_moment'
						metrics_for_plot.(variableName)(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						metrics_for_plot.ankle_flexion_p(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
						metrics_for_plot.ankle_plant_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						metrics_for_plot.ankle_dors_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
						% Joint work metrics
					case 'WORK_POS_TOTAL'
						metrics_for_plot.work_pos(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
					case 'WORK_NEG_TOTAL'
						metrics_for_plot.work_neg(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
						
						% Joint positive power metrics
					case 'PERC_POS_HIP'
						metrics_for_plot.perc_pos_hip(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
					case 'PERC_POS_KNEE'
						metrics_for_plot.perc_pos_knee(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
					case 'PERC_POS_ANKLE'
						metrics_for_plot.perc_pos_ankle(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
						
						% Joint negative power metrics
					case 'PERC_NEG_HIP'
						metrics_for_plot.perc_neg_hip(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
					case 'PERC_NEG_KNEE'
						metrics_for_plot.perc_neg_knee(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
					case 'PERC_NEG_ANKLE'
						metrics_for_plot.perc_neg_ankle(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
				end
			end
		end
	end
end

% Quick loop through variables to convert zeros to NaNs
variable_names = fieldnames(metrics_for_plot);
for vName = 1:length(variable_names)
	metrics_for_plot.(variable_names{vName})(~metrics_for_plot.(variable_names{vName})) = nan;
end

%% Plotting power percentage

conditionLabels = {'TBAS'; 'ARM1'; 'ARM2'; 'ARM3'; 'ARM4'; 'ARM5'};

% Obtain the mean and SD for the hip, knee, and ankle percent power
% contributions. Data output is each column represent a different joint and
% each row represents a different armour condition
% for t = 1:3
% 	if t == 1
% 		variable = 'hip';
% 	elseif t == 2
% 		variable = 'knee';
% 	else
% 		variable = 'ankle';
% 	end	
% 	
% 	% Positive
% 	% Light and slow
% 	[ls_mean(:,t), ls_SD(:,t)] = Mean_SD(metrics_for_plot.(['perc_pos_', variable])(:, 1:2:13));
% 	% Light and fast
% 	[lf_mean(:,t), lf_SD(:,t)] = Mean_SD(metrics_for_plot.(['perc_pos_', variable])(:, 2:2:14));
% 	% Heavy and slow
% 	[hs_mean(:,t), hs_SD(:,t)] = Mean_SD(metrics_for_plot.(['perc_pos_', variable])(:, 15:2:25));
% 	% Heavy and fast
% 	[hf_mean(:,t), lf_SD(:,t)] = Mean_SD(metrics_for_plot.(['perc_pos_', variable])(:, 16:2:26));
% 	
% 	% Negative
% 	% Light and slow
% 	[ls_mean_neg(:,t), ls_SD_neg(:,t)] = Mean_SD(metrics_for_plot.(['perc_neg_', variable])(:, 1:2:13));
% 	% Light and fast
% 	[lf_mean_neg(:,t), lf_SD_neg(:,t)] = Mean_SD(metrics_for_plot.(['perc_neg_', variable])(:, 2:2:14));
% 	% Heavy and slow
% 	[hs_mean_neg(:,t), hs_SD_neg(:,t)] = Mean_SD(metrics_for_plot.(['perc_neg_', variable])(:, 15:2:25));
% 	% Heavy and fast
% 	[hf_mean_neg(:,t), lf_SD_neg(:,t)] = Mean_SD(metrics_for_plot.(['perc_neg_', variable])(:, 16:2:26));
% 	
% end


%% subplots for the joint power percentage analysis

% hz = tight_subplot(4,2,[.07 .03],[.07 .05],[.05 .01]);
% 
% % Pos power
% % 15 kg slow walking
% axes(hz(1))
% bar(light_slow); ylabel('Relative contribution (%)'); t = title('15 kg slow speed');
% set(t, 'Units', 'Normalized', 'Position', [1.03, 1.05, 0], 'fontsize', 16, 'FontName', 'Calibri');
% set(gca,'xticklabel',{'Hip', 'Knee', 'Ankle'}, 'ylim', [0,100], 'fontsize', 14, 'FontName', 'Calibri',...
% 	'ytick', [0, 25,50,75,100], 'box', 'off');
% legend(conditionLabels, 'Location', 'southoutside', 'box', 'off', 'orientation', 'horizontal', 'FontSize', 16);
% 
% % 15 kg fast walking
% axes(hz(3))
% bar(light_fast); ylabel('Relative contribution (%)'); t = title('15 kg fast speed');
% set(t, 'Units', 'Normalized', 'Position', [1.03, 1.05, 0], 'fontsize', 16, 'FontName', 'Calibri');
% set(gca,'xticklabel',{'Hip', 'Knee', 'Ankle'}, 'ylim', [0,100], 'fontsize', 14, 'FontName', 'Calibri',...
% 	'ytick', [0, 25,50,75,100], 'box', 'off');
% 
% % 30 kg slow walking
% axes(hz(5))
% bar(heavy_slow); ylabel('Relative contribution (%)'); t = title('30 kg slow speed');
% set(t, 'Units', 'Normalized', 'Position', [1.03, 1.05, 0], 'fontsize', 16, 'FontName', 'Calibri');
% set(gca,'xticklabel',{'Hip', 'Knee', 'Ankle'}, 'ylim', [0,100], 'fontsize', 14, 'FontName', 'Calibri',...
% 	'ytick', [0, 25,50,75,100], 'box', 'off');
% 
% % 30 kg fast walking
% axes(hz(7))
% bar(heavy_fast); ylabel('Relative contribution (%)'); t = title('30 kg fast speed');
% set(t, 'Units', 'Normalized', 'Position', [1.03, 1.05, 0], 'fontsize', 16, 'FontName', 'Calibri');
% set(gca,'xticklabel',{'Hip', 'Knee', 'Ankle'}, 'ylim', [0,100], 'fontsize', 14, 'FontName', 'Calibri',...
% 	'ytick', [0, 25,50,75,100], 'box', 'off');
% 
% % Neg power
% % 15 kg slow walking
% axes(hz(2))
% bar(light_slow_n);
% set(gca,'xticklabel',{'Hip', 'Knee', 'Ankle'}, 'xcolor', [0 0 0], 'yticklabel', {}, 'ylim', [0,100], 'fontsize', 14, 'FontName', 'Calibri', 'box', 'off',...
% 	'ytick', [0, 25,50,75,100]);
% 
% % 15 kg fast walking
% axes(hz(4))
% bar(light_fast_n);
% set(gca,'xticklabel',{'Hip', 'Knee', 'Ankle'}, 'yticklabel', {}, 'ylim', [0,100], 'fontsize', 14, 'FontName', 'Calibri', 'box', 'off',...
% 	'ytick', [0, 25,50,75,100]);
% 
% % 30 kg slow walking
% axes(hz(6))
% bar(heavy_slow_n);
% set(gca,'xticklabel',{'Hip', 'Knee', 'Ankle'}, 'yticklabel', {}, 'ylim', [0,100], 'fontsize', 14, 'FontName', 'Calibri', 'box', 'off',...
% 	'ytick', [0, 25,50,75,100]);
% 
% % 30 kg fast walking
% axes(hz(8))
% bar(heavy_fast_n);
% set(gca,'xticklabel',{'Hip', 'Knee', 'Ankle'}, 'yticklabel', {}, 'ylim', [0,100], 'fontsize', 14, 'FontName', 'Calibri', 'box', 'off',...
% 	'ytick', [0, 25,50,75,100]);

%% Moment plotting

% Peak values
% Obtain the mean and SD for the hip, knee, and ankle peak moments of
% interest. Data output is each column represent a different joint and
% each row represents a different armour condition

% Initialize vars
ankP15_mean = zeros(6, 2); ankP15_SD = zeros(6,2); ankP30_mean = zeros(6, 2); ankP30_SD =  zeros(6,2);

kneeF15_mean = zeros(6, 2);  kneeF15_SD = zeros(6, 2); kneeF30_mean = zeros(6,2); kneeF30_SD =  zeros(6,2);
kneeE15_mean = zeros(6, 2);  kneeE15_SD = zeros(6, 2); kneeE30_mean = zeros(6,2); kneeE30_SD =  zeros(6,2);

hipE15_mean = zeros(6, 2);  hipE15_SD = zeros(6, 2); hipE30_mean = zeros(6,2); hipE30_SD =  zeros(6,2);

for tt = 1:2
	increment = 2;
	tt_30 = tt + 14;
	
	% Ankle
	[ankP15_mean(:,tt), ankP15_SD(:, tt)] = Mean_SD(abs(metrics_for_plot.ankle_plant_peak(:, (tt+2):increment:(increment*6)+tt)));
	[ankP30_mean(:,tt), ankP30_SD(:, tt)] = Mean_SD(abs(metrics_for_plot.ankle_plant_peak(:, tt_30:increment:(increment*5)+tt_30)));
	
	% Knee
	[kneeF15_mean(:,tt), kneeF15_SD(:, tt)] = Mean_SD(metrics_for_plot.knee_flex_peak(:, tt+2:increment:(increment*6)+tt));
	[kneeF30_mean(:,tt), kneeF30_SD(:, tt)] = Mean_SD(metrics_for_plot.knee_flex_peak(:, tt_30:increment:(increment*5)+tt_30));
	
	% Knee extension
	[kneeE15_mean(:,tt), kneeE15_SD(:, tt)] = Mean_SD(metrics_for_plot.knee_ext_peak(:, tt+2:increment:(increment*6)+tt));
	[kneeE30_mean(:,tt), kneeE30_SD(:, tt)] = Mean_SD(metrics_for_plot.knee_ext_peak(:, tt_30:increment:(increment*5)+tt_30));
	
	% Hip
	[hipE15_mean(:,tt), hipE15_SD(:, tt)] = Mean_SD(abs(metrics_for_plot.hip_ext_peak(:, tt+2:increment:(increment*6)+tt)));
	[hipE30_mean(:,tt), hipE30_SD(:, tt)] = Mean_SD(abs(metrics_for_plot.hip_ext_peak(:, tt_30:increment:(increment*5)+tt_30)));
	
end

% Make dir if it's not there
if ~isdir([BasePath, filesep, 'results_moments'])
	mkdir(BasePath, 'results_moments');
end

%% Save peaks in .mat files and .csv files
cd(fullfile(BasePath, 'results_moments'));
save('Moments.mat', '-struct', 'metrics_for_plot');

createTableCSV(metrics_for_plot, subjectNumber, conditions);

%% Subplots
% Function to bring plots together
ha = tight_subplot(2,3,[.05 .03],[.1 .1],[.07 .01]);

% ANKLE
axes(ha(3))
CreatePlotExtraVar(ankP15_SD(:,1), ankP15_SD(:,2), ankP15_mean(:,1), ankP15_mean(:,2), conditionLabels, [0,3])
t = title('Ankle Plantarflexion'); set(t, 'Units', 'Normalized', 'Position', [0.5, 1, 0]);
legend({'Slow', 'Fast'}, 'Location', 'northeast', 'box', 'off', 'orientation', 'vertical', 'fontsize', 14);
set(gca, 'YTick', [0, 1, 2, 3]);

axes(ha(6))
CreatePlotExtraVar(ankP30_SD(:,1), ankP30_SD(:,2), ankP30_mean(:,1), ankP30_mean(:,2), conditionLabels, [0, 3.5])
set(gca, 'YTick', [0, 1, 2, 3]);

% KNEE
axes(ha(2))
CreatePlotExtraVar(kneeE15_SD(:,1), kneeE15_SD(:,2), kneeE15_mean(:,1)*-1, kneeE15_mean(:,2)*-1, conditionLabels, [0, 1])
t = title('Knee Extension'); set(t, 'Units', 'Normalized', 'Position', [0.5, 1, 0]);
set(gca, 'YTick', [0, 0.5, 1]);

axes(ha(5))
CreatePlotExtraVar(kneeE30_SD(:,1), kneeE30_SD(:,2), kneeE30_mean(:,1)*-1, kneeE30_mean(:,2)*-1, conditionLabels, [0, 2])
set(gca, 'xticklabel', conditionLabels, 'fontsize', 14, 'ylim', [0,1], 'YTick', [0, 0.5, 1], 'FontName', 'Calibri', 'box', 'off');

% HIP
axes(ha(1))
CreatePlotExtraVar(hipE15_SD(:,1), hipE15_SD(:,2), hipE15_mean(:,1), hipE15_mean(:,2), conditionLabels, [0, 4])
y = ylabel({'15 kg'; 'Moment (N.m. kg^-^1)'});  t = title('Hip Extension'); 
set(y, 'Units', 'Normalized', 'Position', [-0.07, 0.5, 0]); set(t, 'Units', 'Normalized', 'Position', [0.5, 1, 0]);

axes(ha(4))
CreatePlotExtraVar(hipE30_SD(:,1), hipE30_SD(:,2), hipE30_mean(:,1), hipE30_mean(:,2), conditionLabels, [0, 4])
y = ylabel({'30 kg'; 'Moment (N.m. kg^-^1)'}); set(y, 'Units', 'Normalized', 'Position', [-0.07, 0.5, 0]);


%% Moment plotting

% Trajectories
% Obtain the mean and SD for the hip, knee, and ankle moments of
% interest. Data output is each column represent a different joint and
% each row represents a different armour condition

% Initialize vars
ank15_mean = zeros(101, 6, 2); ank15_SD = zeros(101, 6,2); ank30_mean = zeros(101, 6, 2); ank30_SD =  zeros(101, 6,2);
ank15Upper = zeros(101, 6, 2); ank15Lower = zeros(101, 6, 2); ank30Upper = zeros(101, 6, 2); ank30Lower = zeros(101, 6, 2);

knee15_mean = zeros(101, 6, 2);  knee15_SD = zeros(101, 6, 2); knee30_mean = zeros(101, 6,2); knee30_SD =  zeros(101, 6,2);
knee15Upper = zeros(101, 6, 2); knee15Lower = zeros(101, 6, 2); knee30Upper = zeros(101, 6, 2); knee30Lower = zeros(101, 6, 2);

hip15_mean = zeros(101, 6, 2);  hip15_SD = zeros(101, 6, 2); hip30_mean = zeros(101, 6,2); hip30_SD =  zeros(101, 6,2);
hip15Upper = zeros(101, 6, 2); hip15Lower = zeros(101, 6, 2); hip30Upper = zeros(101, 6, 2); hip30Lower = zeros(101, 6, 2);

for m = 1:2
	increment = 2;
	m_30 = m + 14;
	
	% Ankle
	[ank15_mean(:, :, m), ank15_SD(:, :, m), ank15Upper(:, :, m), ank15Lower(:, :, m)]...
		= Mean_SD_3D(metrics_for_plot.ankle_angle_r_moment(:, :, (m+2):increment:(increment*6)+m));
	[ank30_mean(:, :, m), ank30_SD(:, :, m), ank30Upper(:, :, m), ank30Lower(:, :, m)]...
		= Mean_SD_3D(metrics_for_plot.ankle_angle_r_moment(:, :, m_30:increment:(increment*5)+m_30));
	
	% Knee
	[knee15_mean(:, :, m), knee15_SD(:, :, m), knee15Upper(:,:,m), knee15Lower(:,:,m)]...
		= Mean_SD_3D(metrics_for_plot.knee_angle_r_moment(:, :,  m+2:increment:(increment*6)+m));
	[knee30_mean(:, :, m), knee30_SD(:, :, m), knee30Upper(:,:,m), knee30Lower(:,:,m)]...
		= Mean_SD_3D(metrics_for_plot.knee_angle_r_moment(:, :, m_30:increment:(increment*5)+m_30));
	
	% Hip
	[hip15_mean(:, :, m), hip15_SD(:, :, m), hip15Upper(:,:,m), hip15Lower(:,:,m)]...
		= Mean_SD_3D(metrics_for_plot.hip_flexion_r_moment(:, :, m+2:increment:(increment*6)+m));
	[hip30_mean(:, :, m), hip30_SD(:, :, m), hip30Upper(:,:,m), hip30Lower(:,:,m)]...
		= Mean_SD_3D(metrics_for_plot.hip_flexion_r_moment(:, :, m_30:increment:(increment*5)+m_30));
end

%% plots
moments = tight_subplot(2,3,[.05 .03],[.15 .1],[.07 .03]);
x = 0:1:100;
cmap = colormap(parula(250));

% Ankle moments
axes(moments(3))
% Plot of 15 kg slow walking
for cName = 1:6
	plotColor = cmap(round(1+40*(cName-1)),:);
	a1 = plot(matfiltfilt(0.01, 6, 4, (ank15_mean(:,cName,1))), 'LineWidth', 1, 'Color',plotColor);
	hold on
	% Plot SD
	[ph,msg]=jbfill(x,matfiltfilt(0.01, 6, 4, ank15Lower(:,cName,1))', matfiltfilt(0.01, 6, 4, ank15Upper(:,cName,1))', plotColor, plotColor ,0,.1);
	% Make it so the SDs won't appear in legend
	set(get(get(ph,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
t1 = title({'Ankle'; 'dorsiflexion (+)       plantarflexion (-)'}); set(t1, 'FontSize', 16, 'Color', [1 1 1]);
set(gca, 'xlim', [0,100], 'Color', [0 0 0], 'xtick', [], 'xcolor', [0 0 0], 'ycolor', [1 1 1], 'ylim', [-3, 1], 'FontName', 'Calibri', 'FontSize', 16, 'box', 'off');
plot([60,60], [-3,1], '-.', 'Color', [1 1 1], 'LineWidth', 2);
set(gcf, 'Color', [0 0 0]);
l = legend(conditionLabels, 'Location', 'northeast', 'box', 'off', 'orientation', 'horizontal', 'fontsize', 18, 'TextColor', [1 1 1]);
% title(l, 'Armour Conditions');

axes(moments(6))
% Plot of 30 kg slow walking
for cName = 1:6
	plotColor = cmap(round(1+40*(cName-1)),:);
	a2 = plot(matfiltfilt(0.01, 6, 4, (ank30_mean(:,cName,1))), 'LineWidth', 1, 'Color',plotColor);
	hold on
	% Plot SD
	[ph,msg]=jbfill(x,matfiltfilt(0.01, 6, 4, ank30Lower(:,cName,1))', matfiltfilt(0.01, 6, 4, ank30Upper(:,cName,1))', plotColor, plotColor ,0,.1);
end
x1 = xlabel('Gait Cycle (%)'); set(x1, 'Color', [1 1 1]);
set(gca, 'xlim', [0,100], 'Color', [0 0 0], 'xcolor', [1 1 1], 'ycolor', [1 1 1], 'ylim', [-3, 1], 'FontName', 'Calibri', 'FontSize', 16, 'box', 'off');
plot([62.4,62.4], [-3,1], '-.', 'Color', [1 1 1], 'LineWidth', 2);

% Knee moments
axes(moments(2))
% Plot of 15 kg slow walking
for cName = 1:6
	plotColor = cmap(round(1+40*(cName-1)),:);
	k1 = plot(matfiltfilt(0.01, 6, 4, (knee15_mean(:,cName,1))), 'LineWidth', 1, 'Color',plotColor);
	hold on
	% Plot SD
	[ph,msg]=jbfill(x,matfiltfilt(0.01, 6, 4, knee15Lower(:,cName,1))', matfiltfilt(0.01, 6, 4, knee15Upper(:,cName,1))', plotColor, plotColor ,0,.1);
end
t2 = title({'Knee'; 'flexion (+)      extension (-)'}); set(t2, 'FontSize', 16, 'Color', [1 1 1]);
set(gca, 'xlim', [0,100], 'Color', [0 0 0], 'xtick', [], 'xcolor', [0 0 0], 'ycolor', [1 1 1], 'ylim', [-1, 2], 'FontName', 'Calibri', 'FontSize', 16, 'box', 'off');
plot([60,60], [-1,2], '-.', 'Color', [0.9 0.9 0.9], 'LineWidth', 2);

% Plot of 30 kg slow walking
axes(moments(5))
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
axes(moments(1))
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
axes(moments(4))
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