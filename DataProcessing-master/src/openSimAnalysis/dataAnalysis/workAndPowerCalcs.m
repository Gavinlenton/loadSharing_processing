function [ID_metrics, ID_metrics_powers] = workAndPowerCalcs(IK_metrics, ID_metrics, conditions, step_time)
%Compute the work and power performed at each joint
%   INPUT - structures containing the IK and ID waveforms for each
%   participant and condition
%   OUTPUT - structure containing the work and power performed at the
%   ankle, knee, and hip joints.

nameMetrics = fieldnames(ID_metrics);
nameMetricsIK = fieldnames(IK_metrics);
ID_metrics_powers = struct();

g = 2;
for i = 2:length(nameMetrics)
	
	for t = 1:length(conditions)
		
		for k = 1:size(ID_metrics.(nameMetrics{i}).(conditions{t}),2)
			
			% Average step time for the condition
			step_time_subject = step_time.(conditions{t})(k);
			stance_time = ID_metrics.toe_off_time.(conditions{t})(k);
			toe_off_frame = floor(stance_time / step_time_subject * 100);
			stride_time = step_time_subject-stance_time;
			
			% Compute joint powers as product of moments and joint angular
			% velocity
			ID_metrics.([nameMetrics{i}(1:end-6),'JOINT_POWER_TOTAL']).(conditions{t})(:,k) = (ID_metrics.(nameMetrics{i}).(conditions{t})(:,k)...
				.* IK_metrics.(nameMetricsIK{g}).(conditions{t})(:,k));
			
			if toe_off_frame > 10
				% Stride power
				ID_metrics.([nameMetrics{i}(1:end-6),'JOINT_POWER_STRIDE']).(conditions{t})(:,k)...
					= (ID_metrics.(nameMetrics{i}).(conditions{t})(toe_off_frame:end,k)...
					.* IK_metrics.(nameMetricsIK{g}).(conditions{t})(toe_off_frame:end, k));
				
				% Stance power
				ID_metrics.([nameMetrics{i}(1:end-6),'JOINT_POWER_STANCE']).(conditions{t})(:,k)...
					= (ID_metrics.(nameMetrics{i}).(conditions{t})(1:toe_off_frame,k)...
					.* IK_metrics.(nameMetricsIK{g}).(conditions{t})(1:toe_off_frame, k));
				
				for kk = 1:3
					% Normal power analysis
					if kk == 1
						
						time = 'JOINT_POWER_TOTAL';
						ID_metrics_powers = power_work_calculations(ID_metrics, ID_metrics_powers, nameMetrics, conditions, i, t, k, time, step_time_subject);
						
						% Just stride time
					elseif kk == 2
						time = 'JOINT_POWER_STRIDE';
						ID_metrics_powers = power_work_calculations(ID_metrics, ID_metrics_powers, nameMetrics, conditions, i, t, k, time, stride_time);
						
						% Just stance time
					else
						time = 'JOINT_POWER_STANCE';
						ID_metrics_powers = power_work_calculations(ID_metrics, ID_metrics_powers, nameMetrics, conditions, i, t, k, time, stance_time);
						
					end
				end
			end
		end
	end
	g = g + 2;
end
end

