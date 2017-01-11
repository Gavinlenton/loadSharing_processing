function [ID_metrics_for_analysis, IK_metrics_for_analysis]  = extractVariables(subjects, metricsID, metricsIK, conditions, ID_metrics_ordered, IK_metrics_ordered)
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
	
	% Check to see if condition is not TBAS
	if ~strcmp(conditionName(1:4), 'TBAS')
		% Define TBAS condition name
		tbasName = regexprep(conditionName, conditionName(1:end-7), 'TBAS');
		notTBAS = 1;
	end
		
	% Loop through subjects
	for sName = 1:length(subjects)
		
		% Determine if subject has that condition name
		if any(strcmp(fieldnames(ID_metrics_ordered.(subjects{sName})), conditionName))
			
			% Loop through the variables of interest
			for k = 1:length(metricsID)
				variableName = metricsID{k};
				
				switch variableName
					
					% Integer variables are stored as subjects (rows) * conditions
					% (columns) and vectors stored as variable * subject (columns) *
					% conditions (z direction)
					case 'hip_flexion_r_moment'
						ID_metrics_for_analysis.(variableName).(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						ID_metrics_for_analysis.([variableName, '_SD'])(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd;
						ID_metrics_for_analysis.hip_flexion_power.(conditionName)(:, sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
						ID_metrics_for_analysis.hip_ext_peak.(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						ID_metrics_for_analysis.hip_flex_peak.(conditionName)(sName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
						% If it's not TBAS then subtract TBAS from the mean
						% values to get the diff
						if notTBAS == 1
							% Subtract values from TBAS
							ID_metrics_for_analysis.([variableName, '_diffTBAS']).(conditionName)(:, sName) =...
								ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean...
								- ID_metrics_ordered.(subjects{sName}).(tbasName).(variableName).mean;
						end
						
					case 'hip_adduction_r_moment'
						ID_metrics_for_analysis.(variableName)(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						ID_metrics_for_analysis.([variableName, '_SD'])(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd;
						ID_metrics_for_analysis.hip_adduction_power(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
						ID_metrics_for_analysis.hip_add_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						ID_metrics_for_analysis.hip_abd_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
					case 'hip_rotation_r_moment'
						ID_metrics_for_analysis.(variableName)(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						ID_metrics_for_analysis.([variableName, '_SD'])(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd;
						ID_metrics_for_analysis.hip_rotation_power(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
						ID_metrics_for_analysis.hip_int_rot_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						ID_metrics_for_analysis.hip_ext_rot_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
					case 'knee_angle_r_moment'
						ID_metrics_for_analysis.(variableName)(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						ID_metrics_for_analysis.([variableName, '_SD'])(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd;
						ID_metrics_for_analysis.knee_flexion_power(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
						ID_metrics_for_analysis.knee_ext_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						ID_metrics_for_analysis.knee_flex_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
					case 'knee_adduction_r_moment'
						ID_metrics_for_analysis.(variableName)(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						ID_metrics_for_analysis.([variableName, '_SD'])(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd;
						ID_metrics_for_analysis.knee_adduction_power(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
						ID_metrics_for_analysis.knee_add_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						ID_metrics_for_analysis.knee_abd_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
					case 'knee_rotation_r_moment'
						ID_metrics_for_analysis.(variableName)(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						ID_metrics_for_analysis.([variableName, '_SD'])(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd;
						ID_metrics_for_analysis.knee_rotation_power(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
						ID_metrics_for_analysis.knee_int_rot_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						ID_metrics_for_analysis.knee_ext_rot_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
					case 'ankle_angle_r_moment'
						ID_metrics_for_analysis.(variableName)(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						ID_metrics_for_analysis.([variableName, '_SD'])(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd;
						ID_metrics_for_analysis.ankle_flexion_power(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
						ID_metrics_for_analysis.ankle_plant_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						ID_metrics_for_analysis.ankle_dors_peak(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
					case 'step_time'
						ID_metrics_for_analysis.(variableName)(:, sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
						
						% Joint work metrics
					case 'WORK_POS_TOTAL'
						ID_metrics_for_analysis.work_pos(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
					case 'WORK_NEG_TOTAL'
						ID_metrics_for_analysis.work_neg(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
						
						% Joint positive power metrics
					case 'PERC_POS_HIP'
						ID_metrics_for_analysis.perc_pos_hip(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
					case 'PERC_POS_KNEE'
						ID_metrics_for_analysis.perc_pos_knee(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
					case 'PERC_POS_ANKLE'
						ID_metrics_for_analysis.perc_pos_ankle(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
						
						% Joint negative power metrics
					case 'PERC_NEG_HIP'
						ID_metrics_for_analysis.perc_neg_hip(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
					case 'PERC_NEG_KNEE'
						ID_metrics_for_analysis.perc_neg_knee(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
					case 'PERC_NEG_ANKLE'
						ID_metrics_for_analysis.perc_neg_ankle(sName, cName) = ID_metrics_ordered.(subjects{sName}).(conditionName).(variableName);
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
						IK_metrics_for_analysis.(variableName)(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						IK_metrics_for_analysis.([variableName, '_SD'])(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd;
						IK_metrics_for_analysis.hip_flexion_vel(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).ANGULAR_VEL;
						IK_metrics_for_analysis.hip_ext_peak(sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						IK_metrics_for_analysis.hip_flex_peak(sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
					case 'hip_adduction_r'
						IK_metrics_for_analysis.(variableName)(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						IK_metrics_for_analysis.([variableName, '_SD'])(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd;
						IK_metrics_for_analysis.hip_adduction_vel(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).ANGULAR_VEL;
						IK_metrics_for_analysis.hip_add_peak(sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						IK_metrics_for_analysis.hip_abd_peak(sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
					case 'hip_rotation_r'
						IK_metrics_for_analysis.(variableName)(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						IK_metrics_for_analysis.([variableName, '_SD'])(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd;
						IK_metrics_for_analysis.hip_rotation_vel(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).ANGULAR_VEL;
						IK_metrics_for_analysis.hip_int_rot_peak(sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						IK_metrics_for_analysis.hip_ext_rot_peak(sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
					case 'knee_angle_r'
						IK_metrics_for_analysis.(variableName)(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						IK_metrics_for_analysis.([variableName, '_SD'])(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd;
						IK_metrics_for_analysis.knee_flexion_vel(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).ANGULAR_VEL;
						IK_metrics_for_analysis.knee_ext_peak(sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						IK_metrics_for_analysis.knee_flex_peak(sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
					case 'knee_adduction_r'
						IK_metrics_for_analysis.(variableName)(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						IK_metrics_for_analysis.([variableName, '_SD'])(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd;
						IK_metrics_for_analysis.knee_adduction_vel(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).ANGULAR_VEL;
						IK_metrics_for_analysis.knee_add_peak(sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						IK_metrics_for_analysis.knee_abd_peak(sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
					case 'knee_rotation_r'
						IK_metrics_for_analysis.(variableName)(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						IK_metrics_for_analysis.([variableName, '_SD'])(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd;
						IK_metrics_for_analysis.knee_rotation_vel(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).ANGULAR_VEL;
						IK_metrics_for_analysis.knee_int_rot_peak(sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						IK_metrics_for_analysis.knee_ext_rot_peak(sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
					case 'ankle_angle_r'
						IK_metrics_for_analysis.(variableName)(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						IK_metrics_for_analysis.([variableName, '_SD'])(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd;
						IK_metrics_for_analysis.ankle_flexion_vel(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).ANGULAR_VEL;
						IK_metrics_for_analysis.ankle_plant_peak(sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						IK_metrics_for_analysis.ankle_dors_peak(sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						
					case 'lumbar_extension'
						IK_metrics_for_analysis.(variableName)(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						IK_metrics_for_analysis.([variableName, '_SD'])(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd;
						IK_metrics_for_analysis.lumbar_ext_peak(sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						IK_metrics_for_analysis.lumbar_flex_peak(sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						
					case 'lumbar_bending'
						IK_metrics_for_analysis.(variableName)(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						IK_metrics_for_analysis.([variableName, '_SD'])(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd;
						IK_metrics_for_analysis.lumbar_rightBend_peak(sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						IK_metrics_for_analysis.lumbar_leftBend_peak(sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
						
					case 'lumbar_rotation'
						IK_metrics_for_analysis.(variableName)(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).mean;
						IK_metrics_for_analysis.([variableName, '_SD'])(:, sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).sd;
						IK_metrics_for_analysis.lumbar_intRot_peak(sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
						IK_metrics_for_analysis.lumbar_extRot_peak(sName, cName) = IK_metrics_ordered.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
				end
			end
		end
	end
end

end

