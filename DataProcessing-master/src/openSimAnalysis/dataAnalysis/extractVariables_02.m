function [ID_metrics_for_analysis, IK_metrics_for_analysis]  = extractVariables_02(subjects, metricsID, metricsIK, conditions, ID_metrics_ordered, IK_metrics_ordered)
% Extracts the metrics of interest from the load sharing data for all participants
% and puts them into structures for further analysis

%   INPUT: subjects - names of subjects for analysis
%          metrics - names of the metrics of interest
%          conditions - names of the conditions to be analysed
%          ID_metrics - structure containing the metrics from ID
%          IK_metrics - structure containing the metrics from IK

%   OUTPUT - structure containing the metrics for further analysis

% Initialize structure
ID_metrics_for_analysis = struct();
IK_metrics_for_analysis = struct();

% Loop through all conditions
for cName = 1:length(conditions)
	
	conditionName = conditions{cName};
	
	% Loop through subjects
	for sName = 1:length(subjects)
		
		% Get condition names for the participant
		conditionNames_participant = fieldnames(ID_metrics_ordered.(subjects{sName}));
		
		% Check to see if condition is not TBAS and NA
		if ~strcmp(conditionName(1:4), 'TBAS') && ~strcmp(conditionName(1:2), 'NA')
			% Define TBAS condition name
			tbasName = regexprep(conditionName, conditionName(1:end-7), 'TBAS');
			
			% Check to make sure that the participant has the TBAS condition
			% as well -  if not then don't subtract TBAS value
			if any(strcmp(conditionNames_participant, tbasName))
				notTBAS = 1;
			else
				notTBAS = 0;
			end
		else
			notTBAS = 0;
		end
		
		% Determine if subject has that condition name
		if any(strcmp(fieldnames(ID_metrics_ordered.(subjects{sName})), conditionName))
			
			% Loop through the variables of interest
			for k = 1:length(metricsID)
				variableName = metricsID{k};
				
				switch variableName
					
					% Integer variables are stored as subjects (rows) * conditions
					% (columns) and vectors stored as variable * subject (columns) *
					% conditions (z direction)
					case 'lumbar_extension_moment'
						ID_metrics_for_analysis.(variableName).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end);
						ID_metrics_for_analysis.([variableName, '_SD']).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd(3:end);
						ID_metrics_for_analysis.lumb_ext_peak.(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						ID_metrics_for_analysis.lumb_flex_peak.(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
						% If it's not TBAS then subtract TBAS from the mean
						% values to get the diff
						if notTBAS == 1
							% Subtract values from TBAS
							ID_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(:, sName) =...
								ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end)...
								- ID_metrics_ordered.(subjects{sName}).(tbasName).(variableName).mean(3:end);
						end
						
					case 'lumbar_bending_moment'
						ID_metrics_for_analysis.(variableName).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end);
						ID_metrics_for_analysis.([variableName, '_SD']).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd(3:end);
						ID_metrics_for_analysis.lumb_rightBend_peak.(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						ID_metrics_for_analysis.lumb_leftBend_peak.(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
						% If it's not TBAS then subtract TBAS from the mean
						% values to get the diff
						if notTBAS == 1
							% Subtract values from TBAS
							ID_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(:, sName) =...
								ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end)...
								- ID_metrics_ordered.(subjects{sName}).(tbasName).(variableName).mean(3:end);
						end
						
					case 'lumbar_rotation_moment'
						ID_metrics_for_analysis.(variableName).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end);
						ID_metrics_for_analysis.([variableName, '_SD']).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd(3:end);
						ID_metrics_for_analysis.lumb_intRot_peak.(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						ID_metrics_for_analysis.lumb_extRot_peak.(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
						% If it's not TBAS then subtract TBAS from the mean
						% values to get the diff
						if notTBAS == 1
							% Subtract values from TBAS
							ID_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(:, sName) =...
								ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end)...
								- ID_metrics_ordered.(subjects{sName}).(tbasName).(variableName).mean(3:end);
						end
						
					case 'hip_flexion_r_moment'
						ID_metrics_for_analysis.(variableName).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end);
						ID_metrics_for_analysis.([variableName, '_SD']).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd(3:end);
						ID_metrics_for_analysis.hip_flexion_power.(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).JOINT_POWER(3:end);
						ID_metrics_for_analysis.hip_ext_peak.(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						ID_metrics_for_analysis.hip_flex_peak.(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
						% If it's not TBAS then subtract TBAS from the mean
						% values to get the diff
						if notTBAS == 1
							% Subtract values from TBAS
							ID_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(:, sName) =...
								ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end)...
								- ID_metrics_ordered.(subjects{sName}).(tbasName).(variableName).mean(3:end);
						end
						
					case 'hip_adduction_r_moment'
						ID_metrics_for_analysis.(variableName).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end);
						ID_metrics_for_analysis.([variableName, '_SD']).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd(3:end);
						ID_metrics_for_analysis.hip_adduction_power.(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).JOINT_POWER(3:end);
						ID_metrics_for_analysis.hip_add_peak.(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						ID_metrics_for_analysis.hip_abd_peak.(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
						if notTBAS == 1
							% Subtract values from TBAS
							ID_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(:, sName) =...
								ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end)...
								- ID_metrics_ordered.(subjects{sName}).(tbasName).(variableName).mean(3:end);
						end
						
					case 'hip_rotation_r_moment'
						ID_metrics_for_analysis.(variableName).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end);
						ID_metrics_for_analysis.([variableName, '_SD']).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd(3:end);
						ID_metrics_for_analysis.hip_rotation_power.(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).JOINT_POWER(3:end);
						ID_metrics_for_analysis.hip_int_rot_peak.(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						ID_metrics_for_analysis.hip_ext_rot_peak.(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
						if notTBAS == 1
							% Subtract values from TBAS
							ID_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(:, sName) =...
								ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end)...
								- ID_metrics_ordered.(subjects{sName}).(tbasName).(variableName).mean(3:end);
						end
						
					case 'knee_angle_r_moment'
						ID_metrics_for_analysis.(variableName).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end);
						ID_metrics_for_analysis.([variableName, '_SD']).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd(3:end);
						ID_metrics_for_analysis.knee_flexion_power.(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).JOINT_POWER(3:end);
						ID_metrics_for_analysis.knee_ext_peak.(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						[pks, locs] = findpeaks(ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(1:20),'MinPeakDistance', 10);
						ID_metrics_for_analysis.knee_flex_peak.(conditionName)(sName) = pks(1);
						
						if notTBAS == 1
							% Subtract values from TBAS
							ID_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(:, sName) =...
								ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end)...
								- ID_metrics_ordered.(subjects{sName}).(tbasName).(variableName).mean(3:end);
						end
						
					case 'knee_adduction_r_moment'
						ID_metrics_for_analysis.(variableName).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end);
						ID_metrics_for_analysis.([variableName, '_SD']).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd(3:end);
						ID_metrics_for_analysis.knee_adduction_power.(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).JOINT_POWER(3:end);
						ID_metrics_for_analysis.knee_add_peak.(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						ID_metrics_for_analysis.knee_abd_peak.(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
						if notTBAS == 1
							% Subtract values from TBAS
							ID_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(:, sName) =...
								ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end)...
								- ID_metrics_ordered.(subjects{sName}).(tbasName).(variableName).mean(3:end);
						end
						
					case 'knee_rotation_r_moment'
						ID_metrics_for_analysis.(variableName).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end);
						ID_metrics_for_analysis.([variableName, '_SD']).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd(3:end);
						ID_metrics_for_analysis.knee_rotation_power.(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).JOINT_POWER(3:end);
						ID_metrics_for_analysis.knee_int_rot_peak.(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						ID_metrics_for_analysis.knee_ext_rot_peak.(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
						if notTBAS == 1
							% Subtract values from TBAS
							ID_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(:, sName) =...
								ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end)...
								- ID_metrics_ordered.(subjects{sName}).(tbasName).(variableName).mean(3:end);
						end
						
					case 'ankle_angle_r_moment'
						ID_metrics_for_analysis.(variableName).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end);
						ID_metrics_for_analysis.([variableName, '_SD']).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd(3:end);
						ID_metrics_for_analysis.ankle_flexion_power.(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).JOINT_POWER(3:end);
						ID_metrics_for_analysis.ankle_plant_peak.(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						ID_metrics_for_analysis.ankle_dors_peak.(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
						if notTBAS == 1
							% Subtract values from TBAS
							ID_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(:, sName) =...
								ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end)...
								- ID_metrics_ordered.(subjects{sName}).(tbasName).(variableName).mean(3:end);
						end
						
					case 'step_time'
						ID_metrics_for_analysis.(variableName).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
						
						% Joint work metrics
					case 'WORK_POS_TOTAL'
						ID_metrics_for_analysis.(variableName).(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
						
						if notTBAS == 1
							% Subtract values from TBAS
							ID_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(sName) =...
								ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName)...
								- ID_metrics_ordered.(subjects{sName}).(tbasName).(variableName);
						end
						
					case 'WORK_NEG_TOTAL'
						ID_metrics_for_analysis.(variableName).(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
						
						if notTBAS == 1
							% Subtract values from TBAS
							ID_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(sName) =...
								ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName)...
								- ID_metrics_ordered.(subjects{sName}).(tbasName).(variableName);
						end
						
						% Joint positive power metrics
					case 'PERC_POS_HIP'
						ID_metrics_for_analysis.(variableName).(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
						
						if notTBAS == 1
							% Subtract values from TBAS
							ID_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(sName) =...
								ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName)...
								- ID_metrics_ordered.(subjects{sName}).(tbasName).(variableName);
						end
						
					case 'PERC_POS_KNEE'
						ID_metrics_for_analysis.(variableName).(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
						
						if notTBAS == 1
							% Subtract values from TBAS
							ID_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(sName) =...
								ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName)...
								- ID_metrics_ordered.(subjects{sName}).(tbasName).(variableName);
						end
						
					case 'PERC_POS_ANKLE'
						ID_metrics_for_analysis.(variableName).(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
						
						if notTBAS == 1
							% Subtract values from TBAS
							ID_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(sName) =...
								ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName)...
								- ID_metrics_ordered.(subjects{sName}).(tbasName).(variableName);
						end
						
						% Joint negative power metrics
					case 'PERC_NEG_HIP'
						ID_metrics_for_analysis.(variableName).(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
						
						if notTBAS == 1
							% Subtract values from TBAS
							ID_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(sName) =...
								ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName)...
								- ID_metrics_ordered.(subjects{sName}).(tbasName).(variableName);
						end
						
					case 'PERC_NEG_KNEE'
						ID_metrics_for_analysis.(variableName).(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
						
						if notTBAS == 1
							% Subtract values from TBAS
							ID_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(sName) =...
								ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName)...
								- ID_metrics_ordered.(subjects{sName}).(tbasName).(variableName);
						end
						
					case 'PERC_NEG_ANKLE'
						ID_metrics_for_analysis.(variableName).(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
						
						if notTBAS == 1
							% Subtract values from TBAS
							ID_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(sName) =...
								ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName)...
								- ID_metrics_ordered.(subjects{sName}).(tbasName).(variableName);
						end
				end
			end
			
			% IK data extraction
			for k = 1:length(metricsIK)
				variableName = metricsIK{k};
				
				switch variableName
					
					% Integer variables are stored as subjects (rows) * conditions
					% (columns) and vectors stored as variable * subject (columns) *
					% conditions (z direction)
					case 'hip_flexion_r'
						IK_metrics_for_analysis.(variableName).(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end);
						IK_metrics_for_analysis.([variableName, '_SD']).(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd(3:end);
						IK_metrics_for_analysis.hip_flexion_vel.(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).ANGULAR_VEL(3:end);
						IK_metrics_for_analysis.hip_ext_peak.(conditionName)(sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						IK_metrics_for_analysis.hip_flex_peak.(conditionName)(sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
						if notTBAS == 1
							% Subtract values from TBAS
							IK_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(:, sName) =...
								IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end)...
								- IK_metrics_ordered.(subjects{sName}).(tbasName).(variableName).mean(3:end);
						end
						
					case 'hip_adduction_r'
						IK_metrics_for_analysis.(variableName).(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end);
						IK_metrics_for_analysis.([variableName, '_SD']).(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd(3:end);
						IK_metrics_for_analysis.hip_adduction_vel.(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).ANGULAR_VEL(3:end);
						IK_metrics_for_analysis.hip_add_peak.(conditionName)(sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						IK_metrics_for_analysis.hip_abd_peak.(conditionName)(sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
						if notTBAS == 1
							% Subtract values from TBAS
							IK_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(:, sName) =...
								IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end)...
								- IK_metrics_ordered.(subjects{sName}).(tbasName).(variableName).mean(3:end);
						end
						
					case 'hip_rotation_r'
						IK_metrics_for_analysis.(variableName).(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end);
						IK_metrics_for_analysis.([variableName, '_SD']).(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd(3:end);
						IK_metrics_for_analysis.hip_rotation_vel.(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).ANGULAR_VEL(3:end);
						IK_metrics_for_analysis.hip_int_rot_peak.(conditionName)(sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						IK_metrics_for_analysis.hip_ext_rot_peak.(conditionName)(sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
						if notTBAS == 1
							% Subtract values from TBAS
							IK_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(:, sName) =...
								IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end)...
								- IK_metrics_ordered.(subjects{sName}).(tbasName).(variableName).mean(3:end);
						end
						
					case 'knee_angle_r'
						IK_metrics_for_analysis.(variableName).(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end);
						IK_metrics_for_analysis.([variableName, '_SD']).(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd(3:end);
						IK_metrics_for_analysis.knee_flexion_vel.(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).ANGULAR_VEL(3:end);
						IK_metrics_for_analysis.knee_ext_peak.(conditionName)(sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						
						% Find a new peak knee flexion moment in early
						% stance
						IK_metrics_for_analysis.knee_flex_peak.(conditionName)(sName) = max(IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(1:40));
						
						if notTBAS == 1
							% Subtract values from TBAS
							IK_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(:, sName) =...
								IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end)...
								- IK_metrics_ordered.(subjects{sName}).(tbasName).(variableName).mean(3:end);
						end
						
					case 'knee_adduction_r'
						IK_metrics_for_analysis.(variableName).(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end);
						IK_metrics_for_analysis.([variableName, '_SD']).(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd(3:end);
						IK_metrics_for_analysis.knee_adduction_vel.(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).ANGULAR_VEL(3:end);
						IK_metrics_for_analysis.knee_add_peak.(conditionName)(sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						IK_metrics_for_analysis.knee_abd_peak.(conditionName)(sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
						if notTBAS == 1
							% Subtract values from TBAS
							IK_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(:, sName) =...
								IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end)...
								- IK_metrics_ordered.(subjects{sName}).(tbasName).(variableName).mean(3:end);
						end
						
					case 'knee_rotation_r'
						IK_metrics_for_analysis.(variableName).(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end);
						IK_metrics_for_analysis.([variableName, '_SD']).(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd(3:end);
						IK_metrics_for_analysis.knee_rotation_vel.(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).ANGULAR_VEL(3:end);
						IK_metrics_for_analysis.knee_int_rot_peak.(conditionName)(sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						IK_metrics_for_analysis.knee_ext_rot_peak.(conditionName)(sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
						if notTBAS == 1
							% Subtract values from TBAS
							IK_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(:, sName) =...
								IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end)...
								- IK_metrics_ordered.(subjects{sName}).(tbasName).(variableName).mean(3:end);
						end
						
					case 'ankle_angle_r'
						IK_metrics_for_analysis.(variableName).(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end);
						IK_metrics_for_analysis.([variableName, '_SD']).(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd(3:end);
						IK_metrics_for_analysis.ankle_flexion_vel.(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).ANGULAR_VEL(3:end);
						IK_metrics_for_analysis.ankle_plant_peak.(conditionName)(sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						IK_metrics_for_analysis.ankle_dors_peak.(conditionName)(sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
						if notTBAS == 1
							% Subtract values from TBAS
							IK_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(:, sName) =...
								IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end)...
								- IK_metrics_ordered.(subjects{sName}).(tbasName).(variableName).mean(3:end);
						end
						
					case 'lumbar_extension'
						IK_metrics_for_analysis.(variableName).(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end);
						IK_metrics_for_analysis.([variableName, '_SD']).(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd(3:end);
						IK_metrics_for_analysis.lumbar_ext_peak.(conditionName)(sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						IK_metrics_for_analysis.lumbar_flex_peak.(conditionName)(sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						IK_metrics_for_analysis.lumbar_flex_mean.(conditionName)(sName) = mean(IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end));
						
						if notTBAS == 1
							% Subtract values from TBAS
							IK_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(:, sName) =...
								IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end)...
								- IK_metrics_ordered.(subjects{sName}).(tbasName).(variableName).mean(3:end);
						end
						
					case 'lumbar_bending'
						IK_metrics_for_analysis.(variableName).(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end);
						IK_metrics_for_analysis.([variableName, '_SD']).(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd(3:end);
						IK_metrics_for_analysis.lumbar_rightBend_peak.(conditionName)(sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						IK_metrics_for_analysis.lumbar_leftBend_peak.(conditionName)(sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						
						if notTBAS == 1
							% Subtract values from TBAS
							IK_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(:, sName) =...
								IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end)...
								- IK_metrics_ordered.(subjects{sName}).(tbasName).(variableName).mean(3:end);
						end
						
					case 'lumbar_rotation'
						IK_metrics_for_analysis.(variableName).(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end);
						IK_metrics_for_analysis.([variableName, '_SD']).(conditionName)(:, sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd(3:end);
						IK_metrics_for_analysis.lumbar_intRot_peak.(conditionName)(sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						IK_metrics_for_analysis.lumbar_extRot_peak.(conditionName)(sName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						
						if notTBAS == 1
							% Subtract values from TBAS
							IK_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(:, sName) =...
								IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean(3:end)...
								- IK_metrics_ordered.(subjects{sName}).(tbasName).(variableName).mean(3:end);
						end
						
				end
			end
		end
	end
end

end

