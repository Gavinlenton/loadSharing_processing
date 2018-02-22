function metrics_for_plot = extractVarFromStrct(metrics, subjects, fileID)
%Function to select variables depending on whether analysis is for
%kinematic, joint moment, point kinematics, or reaction force data
%   Input structure containing metrics of interest, and cell array of subjects and conditions for analysis. Function then extracts
%   the fieldnames of all variables and depending on the identification
%   label, will extract pertinent information into an output structure

conditions = fieldnames(metrics.(subjects{1}))';
all_variables = fieldnames(metrics.(subjects{1}).(conditions{1}));
metrics_for_plot = struct();

% Loop through all conditions
for cName = 1:length(conditions)
	
	conditionName = conditions{cName};
	
	% Loop through subjects
	for sName = 1:length(subjects)
		
		% Determine if subject has that condition name
		if any(strcmp(fieldnames(metrics.(subjects{sName})), conditionName))
			
			switch fileID
				
				% Inverse dynamics
				case 'ID'
					
					% Loop through the variables of interest
					for k = 1:length(all_variables)
						variableName = all_variables{k};
						
						switch variableName
							
							% Integer variables are stored as subjects (rows) * conditions
							% (columns) and vectors stored as variable * subject (columns) *
							% conditions (z direction)
							
							case 'lumbar_extension_moment'
								metrics_for_plot.(variableName)(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).mean;
								metrics_for_plot.lumbar_flexion_p(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
								metrics_for_plot.lumbar_ext_peak(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
								metrics_for_plot.lumbar_flex_peak(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
								
							case 'lumbar_bending_moment'
								metrics_for_plot.(variableName)(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).mean;
								metrics_for_plot.lumbar_bending_p(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
								metrics_for_plot.hip_add_peak(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
								metrics_for_plot.hip_abd_peak(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
								
							case 'lumbar_rotation_moment'
								metrics_for_plot.(variableName)(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).mean;
								metrics_for_plot.lumbar_rot_p(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
								metrics_for_plot.hip_int_rot_peak(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
								metrics_for_plot.hip_ext_rot_peak(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
							
							case 'hip_flexion_r_moment'
								metrics_for_plot.(variableName)(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).mean;
								metrics_for_plot.hip_flexion_p(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
								metrics_for_plot.hip_ext_peak(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
								metrics_for_plot.hip_flex_peak(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
								
							case 'hip_adduction_r_moment'
								metrics_for_plot.(variableName)(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).mean;
								metrics_for_plot.hip_add_p(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
								metrics_for_plot.hip_add_peak(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
								metrics_for_plot.hip_abd_peak(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
								
							case 'hip_rotation_r_moment'
								metrics_for_plot.(variableName)(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).mean;
								metrics_for_plot.hip_rot_p(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
								metrics_for_plot.hip_int_rot_peak(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
								metrics_for_plot.hip_ext_rot_peak(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
								
							case 'knee_angle_r_moment'
								metrics_for_plot.(variableName)(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).mean;
								metrics_for_plot.knee_flexion_p(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
								metrics_for_plot.knee_ext_peak(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
								metrics_for_plot.knee_flex_peak(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
								
							case 'knee_adduction_r_moment'
								metrics_for_plot.(variableName)(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).mean;
								metrics_for_plot.knee_add_p(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
								metrics_for_plot.knee_add_peak(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
								
							case 'knee_rotation_r_moment'
								metrics_for_plot.(variableName)(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).mean;
								metrics_for_plot.knee_int_rot_peak(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
								metrics_for_plot.knee_ext_rot_peak(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
								metrics_for_plot.knee_rot_p(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
								
							case 'ankle_angle_r_moment'
								metrics_for_plot.(variableName)(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).mean;
								metrics_for_plot.ankle_flexion_p(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).JOINT_POWER;
								metrics_for_plot.ankle_plant_peak(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).max_min_range(2);
								metrics_for_plot.ankle_dors_peak(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).max_min_range(1);
								
								% Joint work metrics
							case 'WORK_POS_TOTAL'
								metrics_for_plot.work_pos(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName);
							case 'WORK_NEG_TOTAL'
								metrics_for_plot.work_neg(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName);
								
								% Joint positive power metrics
							case 'PERC_POS_HIP'
								metrics_for_plot.perc_pos_hip(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName);
							case 'PERC_POS_KNEE'
								metrics_for_plot.perc_pos_knee(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName);
							case 'PERC_POS_ANKLE'
								metrics_for_plot.perc_pos_ankle(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName);
								
								% Joint negative power metrics
							case 'PERC_NEG_HIP'
								metrics_for_plot.perc_neg_hip(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName);
							case 'PERC_NEG_KNEE'
								metrics_for_plot.perc_neg_knee(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName);
							case 'PERC_NEG_ANKLE'
								metrics_for_plot.perc_neg_ankle(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName);
						end
					end
					
					% Point kinematics
				case 'PK'
					
					% Loop through the variables of interest
					for k = 1:length(all_variables)
						variableName = all_variables{k};
						
						% Change what I extract based on variable name
						switch variableName
							
							% Vertical position mean waveform, mean value,
							% and horizontal displacement ROM
							case 'pos'
								metrics_for_plot.(variableName)(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).('Y').mean;
								metrics_for_plot.('vPosMean')(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).('Y').('vPosMean');
								metrics_for_plot.('horzDispROM')(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).('X').max_min_range(3);
								% Range of COM displacement vertical
							case 'range_COM_disp_vert'
								metrics_for_plot.(variableName)(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName);
							case 'jointStiffnessMean'
								metrics_for_plot.(variableName)(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName);
								
						end
					end
					
					% Reaction force
				case 'GRF'
					
					% Loop through the variables of interest
					for k = 1:length(all_variables)
						variableName = all_variables{k};
						
						switch variableName
							
							case 'ground_force1_pz' % Take range of M/L GRF
								metrics_for_plot.(variableName)(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName).('max_min_range')(3);
								
							case 'powerCOM'
								metrics_for_plot.(variableName)(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName);
								
							case 'powerID'
								metrics_for_plot.(variableName)(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName);
								
							case 'powerSoftTissue'
								metrics_for_plot.(variableName)(:, sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName);
								
							case 'workCollisionID'
								metrics_for_plot.(variableName)(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName);
								metrics_for_plot.('workCollisionCOM')(sName, cName) = metrics.(subjects{sName}).(conditionName).('workCollisionCOM');
								metrics_for_plot.('workCollisionST')(sName, cName) = metrics.(subjects{sName}).(conditionName).('workCollisionST');
								
							case 'workReboundID'
								metrics_for_plot.(variableName)(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName);
								metrics_for_plot.('workReboundCOM')(sName, cName) = metrics.(subjects{sName}).(conditionName).('workReboundCOM');
								metrics_for_plot.('workReboundST')(sName, cName) = metrics.(subjects{sName}).(conditionName).('workReboundST');
								
							case 'workPreloadID'
								metrics_for_plot.(variableName)(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName);
								metrics_for_plot.('workPreloadCOM')(sName, cName) = metrics.(subjects{sName}).(conditionName).('workPreloadCOM');
								metrics_for_plot.('workPreloadST')(sName, cName) = metrics.(subjects{sName}).(conditionName).('workPreloadST');
											
							case 'workToeoffID'
								metrics_for_plot.(variableName)(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName);
								metrics_for_plot.('workToeoffCOM')(sName, cName) = metrics.(subjects{sName}).(conditionName).('workToeoffCOM');
								metrics_for_plot.('workToeoffST')(sName, cName) = metrics.(subjects{sName}).(conditionName).('workToeoffST');
								
							case 'posCOMwork'
								metrics_for_plot.(variableName)(sName, cName) = metrics.(subjects{sName}).(conditionName).(variableName);
								metrics_for_plot.('negCOMwork')(sName, cName) = metrics.(subjects{sName}).(conditionName).('negCOMwork');
								metrics_for_plot.('netCOMwork')(sName, cName) = metrics.(subjects{sName}).(conditionName).('netCOMwork');
						end
					end
			end
		end
	end
end

