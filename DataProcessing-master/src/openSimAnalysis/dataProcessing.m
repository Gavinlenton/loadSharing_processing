%% Script to Process data from Load sharing trials

% IK, ID, muscle analysis, and SO taken from Batch OpenSim Processing Scripts (BOPS)
% Copyright (C) 2015 Alice Mantoan, Monica Reggiani
% <https://simtk.org/home/bops>
% Please acknowledge the authors of BOPS when using this script.
clc; clear;

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%% Load the data
BasePath=uigetdir('../../../', 'Select Elaborated Data Folder');
cd(BasePath);
load('IDMetrics_all.mat');

% Sort IDMetrics so it's in correct order
subjects = fieldnames(ID_metrics);
newOrder = {'NA_slow', 'NA_fast', 'TBAS15_slow', 'TBAS15_fast', 'CRYE15_slow', 'CRYE15_fast',...
	'TYR15_slow', 'TYR15_fast', 'USMC15_slow', 'USMC15_fast', 'CORE15_slow', 'CORE15_fast',...
	'SORD15_slow', 'SORD15_fast', 'TBAS30_slow', 'TBAS30_fast', 'CRYE30_slow', 'CRYE30_fast',...
	'TYR30_slow', 'TYR30_fast', 'USMC30_slow', 'USMC30_fast', 'CORE30_slow', 'CORE30_fast',...
	'SORD30_slow', 'SORD30_fast'};

for t = 1:length(subjects)
	ID_metrics_ordered.(subjects{t}) = orderfields(ID_metrics.(subjects{t}), newOrder);
end

% Determine how many condition fields exist
conditions = fieldnames(ID_metrics_ordered.(subjects{1}))';

% Loop through all conditions
for i = 1:length(conditions)
	
	conditionName = conditions{i};

	% Loop through armour and speed conditions
	for ii = 1:length(subjects)
		
		% Determine if subject has that condition name
		if any(strcmp(fieldnames(ID_metrics_ordered.(subjects{ii})), conditionName))
			
			% Determine the names of the conditions which have ROM data for the
			% subject
			all_variables = fieldnames(ID_metrics_ordered.(subjects{ii}).(conditionName));
			all_variables(ismember(all_variables,{'step_time',...
				'WORK_POS_TOTAL', 'WORK_NEG_TOTAL', 'lumbar_extension_moment',...
				'lumbar_bending_moment', 'lumbar_rotation_moment'}))=[];
			
			% Loop through the variables of interest
			for k = 1:length(all_variables)
				variableName = all_variables{k};
				
				switch variableName
					
					% Integer variables are stored as subjects (rows) * conditions
					% (columns) and vectors stored as variable * subject (columns) *
					% conditions (z direction)
					case 'hip_flexion_r_moment'
						hip_flexion_m(:, ii, i) = ID_metrics_ordered.(subjects{ii}).(conditionName).(variableName).mean;
						hip_flexion_p(:, ii, i) = ID_metrics_ordered.(subjects{ii}).(conditionName).(variableName).JOINT_POWER;
						hip_ext_peak(ii, i) = ID_metrics_ordered.(subjects{ii}).(conditionName).(variableName).max_min_range(2);
						
					case 'knee_angle_r_moment'
						knee_flexion_m(:, ii, i) = ID_metrics_ordered.(subjects{ii}).(conditionName).(variableName).mean;
						knee_flexion_p(:, ii, i) = ID_metrics_ordered.(subjects{ii}).(conditionName).(variableName).JOINT_POWER;
						knee_ext_peak(ii, i) = ID_metrics_ordered.(subjects{ii}).(conditionName).(variableName).max_min_range(2);
						knee_flex_peak(ii, i) = ID_metrics_ordered.(subjects{ii}).(conditionName).(variableName).max_min_range(1);
						
					case 'knee_adduction_r_moment'
						knee_flexion_m(:, ii, i) = ID_metrics_ordered.(subjects{ii}).(conditionName).(variableName).mean;
						knee_add_peak(ii, i) = ID_metrics_ordered.(subjects{ii}).(conditionName).(variableName).max_min_range(2);
						
					case 'ankle_angle_r_moment'
						ankle_flexion_m(:, ii, i) = ID_metrics_ordered.(subjects{ii}).(conditionName).(variableName).mean;
						ankle_flexion_p(:, ii, i) = ID_metrics_ordered.(subjects{ii}).(conditionName).(variableName).JOINT_POWER;
						ankle_plant_peak(ii, i) = ID_metrics_ordered.(subjects{ii}).(conditionName).(variableName).max_min_range(2);
						ankle_dors_peak(ii, i) = ID_metrics_ordered.(subjects{ii}).(conditionName).(variableName).max_min_range(1);
						
					% Joint work metrics
					case 'WORK_POS_TOTAL'
						work_pos(ii, i) = ID_metrics_ordered.(subjects{ii}).(conditionName).(variableName);
					case 'WORK_NEG_TOTAL'
						work_neg(ii, i) = ID_metrics_ordered.(subjects{ii}).(conditionName).(variableName);
					
					% Joint positive power metrics
					case 'PERC_POS_HIP'
						perc_pos_hip(ii, i) = ID_metrics_ordered.(subjects{ii}).(conditionName).(variableName);
					case 'PERC_POS_KNEE'
						perc_pos_knee(ii, i) = ID_metrics_ordered.(subjects{ii}).(conditionName).(variableName);
					case 'PERC_POS_ANKLE'
						perc_pos_ankle(ii, i) = ID_metrics_ordered.(subjects{ii}).(conditionName).(variableName);
						
					% Joint negative power metrics
					case 'PERC_NEG_HIP'
						perc_neg_hip(ii, i) = ID_metrics_ordered.(subjects{ii}).(conditionName).(variableName);
					case 'PERC_NEG_KNEE'
						perc_neg_knee(ii, i) = ID_metrics_ordered.(subjects{ii}).(conditionName).(variableName);
					case 'PERC_NEG_ANKLE'
						perc_neg_ankle(ii, i) = ID_metrics_ordered.(subjects{ii}).(conditionName).(variableName);
				end
			end
		end		
	end
end

%% Plotting power percentage

conditionLabels = {'TBAS'; 'CRYE'; 'TYR'; 'USMC'; 'CORE'; 'SORD'};

% Manipulate the data so it's in a format I can plot with
light_slow = [perc_pos_hip(1, 3:2:13); perc_pos_knee(1, 3:2:13); perc_pos_ankle(1, 3:2:13)];
light_fast = [perc_pos_hip(1, 4:2:14); perc_pos_knee(1, 4:2:14); perc_pos_ankle(1, 4:2:14)];
heavy_slow = [perc_pos_hip(1, 15:2:25); perc_pos_knee(1, 15:2:25); perc_pos_ankle(1, 15:2:25)];
heavy_fast = [perc_pos_hip(1, 16:2:26); perc_pos_knee(1, 16:2:26); perc_pos_ankle(1, 16:2:26)];

light_slow_n = [perc_neg_hip(1, 3:2:13); perc_neg_knee(1, 3:2:13); perc_neg_ankle(1, 3:2:13)];
light_fast_n = [perc_neg_hip(1, 4:2:14); perc_neg_knee(1, 4:2:14); perc_neg_ankle(1, 4:2:14)];
heavy_slow_n = [perc_neg_hip(1, 15:2:25); perc_neg_knee(1, 15:2:25); perc_neg_ankle(1, 15:2:25)];
heavy_fast_n = [perc_neg_hip(1, 16:2:26); perc_neg_knee(1, 16:2:26); perc_neg_ankle(1, 16:2:26)];

%% subplots for the joint power percentage analysis

% Pos power
% 15 kg slow walking
subplot(4,2,1)
bar(light_slow); ylabel('Relative contribution (%)'); title('15 kg slow speed');
set(gca,'xticklabel',{'Hip', 'Knee', 'Ankle'}, 'ylim', [0,100], 'fontsize', 14);
legend(conditionLabels, 'Location', 'southoutside', 'box', 'off', 'orientation', 'horizontal');

% 15 kg fast walking
subplot(4,2,3)
bar(light_fast); ylabel('Relative contribution (%)'); title('15 kg fast speed');
set(gca,'xticklabel',{'Hip', 'Knee', 'Ankle'}, 'ylim', [0,100], 'fontsize', 14);

% 30 kg slow walking
subplot(4,2,5)
bar(heavy_slow); ylabel('Relative contribution (%)'); title('30 kg slow speed');
set(gca,'xticklabel',{'Hip', 'Knee', 'Ankle'}, 'ylim', [0,100], 'fontsize', 14);

% 30 kg fast walking
subplot(4,2,7)
bar(heavy_fast); ylabel('Relative contribution (%)'); title('30 kg fast speed');
set(gca,'xticklabel',{'Hip', 'Knee', 'Ankle'}, 'ylim', [0,100], 'fontsize', 14);

% Neg power
% 15 kg slow walking
subplot(4,2,2)
bar(light_slow_n); ylabel('Relative contribution (%)'); title('15 kg slow speed');
set(gca,'xticklabel',{'Hip', 'Knee', 'Ankle'}, 'ylim', [0,100], 'fontsize', 14);

% 15 kg fast walking
subplot(4,2,4)
bar(light_fast_n); ylabel('Relative contribution (%)'); title('15 kg fast speed');
set(gca,'xticklabel',{'Hip', 'Knee', 'Ankle'}, 'ylim', [0,100], 'fontsize', 14);

% 30 kg slow walking
subplot(4,2,6)
bar(heavy_slow_n); ylabel('Relative contribution (%)'); title('30 kg slow speed');
set(gca,'xticklabel',{'Hip', 'Knee', 'Ankle'}, 'ylim', [0,100], 'fontsize', 14);

% 30 kg fast walking
subplot(4,2,8)
bar(heavy_fast_n); ylabel('Relative contribution (%)'); title('30 kg fast speed');
set(gca,'xticklabel',{'Hip', 'Knee', 'Ankle'}, 'ylim', [0,100], 'fontsize', 14);

%% Moment plotting

% Manipulate the data so it's in a format I can plot with
ank_plant_slow = [ankle_plant_peak(1, 3:2:13); ankle_plant_peak(1, 15:2:25)];
ank_plant_fast = [ankle_plant_peak(1, 4:2:14); ankle_plant_peak(1, 16:2:26)];

knee_flex_slow = [knee_flex_peak(1, 3:2:13); knee_flex_peak(1, 15:2:25)];
knee_flex_fast = [knee_flex_peak(1, 4:2:14); knee_flex_peak(1, 16:2:26)];

hip_ext_slow = [hip_ext_peak(1, 3:2:13); hip_ext_peak(1, 15:2:25)];
hip_ext_fast = [hip_ext_peak(1, 4:2:14); hip_ext_peak(1, 16:2:26)];

test = [0.5:1:5.5];

%% Subplots
subplot(3,2,1)
bar(ank_plant_slow'*-1); ylabel('Plantar flexion moment (N kg^-^1)'); title('Slow speed');
set(gca, 'xticklabel', {}, 'fontsize', 12, 'ylim', [0,2.5]);
legend({'15 kg', '30 kg'}, 'Location', 'southoutside', 'box', 'off', 'orientation', 'horizontal');

subplot(3,2,2)
bar(ank_plant_fast'*-1); ylabel('Plantar flexion moment (N kg^-^1)'); title('Fast speed');
set(gca, 'xticklabel', {}, 'fontsize', 12, 'ylim', [0,2.5]);

subplot(3,2,3)
bar(knee_flex_slow'); ylabel('Knee flexion moment (N kg^-^1)');
set(gca, 'xticklabel', {}, 'fontsize', 12, 'ylim', [0,1.5]);

subplot(3,2,4)
bar(knee_flex_fast'); ylabel('Knee flexion moment (N kg^-^1)');
set(gca, 'xticklabel', {}, 'fontsize', 12, 'ylim', [0,1.5]);

subplot(3,2,5)
bar(hip_ext_slow'*-1); ylabel('Hip extension moment (N kg^-^1)');
set(gca,'xticklabel', conditionLabels, 'fontsize', 12, 'ylim', [0,3]);

subplot(3,2,6)
bar(hip_ext_fast'*-1); ylabel('Hip extension moment (N kg^-^1)');
set(gca,'xticklabel', conditionLabels, 'fontsize', 12, 'ylim', [0,3]);

