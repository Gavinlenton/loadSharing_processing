function ID_metrics = workAndPowerCalcs(IK_metrics, ID_metrics, conditions, step_time)
%Compute the work and power performed at each joint
%   INPUT - structures containing the IK and ID waveforms for each
%   participant and condition
%   OUTPUT - structure containing the work and power performed at the
%   ankle, knee, and hip joints.

nameMetrics = fieldnames(ID_metrics);
nameMetricsIK = fieldnames(IK_metrics);
dt = 1/100;

g = 2;
for i = 1:length(nameMetrics)
	
	for t = 1:length(conditions)
		
		% Average step time for the condition
		step_time_subject = nanmean(step_time(1,:,t));
		
		% Compute joint powers as product of moments and joint angular
		% velocity
		ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER') = (ID_metrics.(nameMetrics{i}).(conditions{t}).mean...
			.* IK_metrics.(nameMetricsIK{g}).(conditions{t}).mean);
		
		abovePositiveThreshold = ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER') >= 0;
		spanLocs = bwlabel(abovePositiveThreshold);   %identify contiguous ones
		spanLength = regionprops(abovePositiveThreshold, 'area');  %length of each span
		spanLength = [ spanLength.Area];
		goodSpans = find(spanLength>=5);   %get only spans of 5+ frames
		positiveInSpans = ismember(spanLocs, goodSpans);  %indices of these spans
		ID_metrics.(nameMetrics{i}).(conditions{t}).('positive_power_index') = positiveInSpans;
		
		% Negative periods of interest
		belowPositiveThreshold = ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER') < 0;
		spanLocsBelow = bwlabel(belowPositiveThreshold);   %identify contiguous ones
		spanLengthBelow = regionprops(belowPositiveThreshold, 'area');  %length of each span
		spanLengthBelow = [ spanLengthBelow.Area];
		goodSpansBelow = find(spanLengthBelow>=5);   %get only spans of 10+ frames
		negativeInSpans = ismember(spanLocsBelow, goodSpansBelow);  %indices of these spans
		ID_metrics.(nameMetrics{i}).(conditions{t}).('negative_power_index') = negativeInSpans;
		
		% POWER BURSTS - define burst periods as described in Winter
		% 1987 text.
		
		% HIP
		if strcmp(nameMetrics{i}, 'hip_flexion_r_moment')
			
			% Positive work
			H1 = ismember(spanLocs, goodSpans(1));
			work_H1 = trapz(ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER')(H1,1));
			H3 = ismember(spanLocs, goodSpans(2));
			work_H3 = trapz(ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER')(H3,1));
			
			% Sum positive work to get total positive work
			work_pos_hip =  work_H1 + work_H3;
			% Divide by sampling freq to get normalised work
			work_pos_hip_norm = work_pos_hip*dt;
			
			% Multiply by 2 to get both limbs (assumed symmetry)
			work_pos_hip_bothLimbs = work_pos_hip_norm *2;
			power_pos_hip = work_pos_hip_bothLimbs / step_time_subject;
			
			% Negative work
			if numel(goodSpansBelow) == 1
				H2 = ismember(spanLocsBelow, goodSpansBelow(1));
			else
				H2 = ismember(spanLocsBelow, goodSpansBelow(2));
			end
			work_neg_hip = trapz(ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER')(H2,1));
			work_neg_hip_norm = work_neg_hip * dt;
			
			% Multiply by 2 to get both limbs (assumed symmetry)
			work_neg_hip_bothLimbs = work_neg_hip_norm *2;
			power_neg_hip = work_neg_hip_bothLimbs / step_time_subject;
			
			% Net work
			net_work_hip = work_pos_hip_norm + work_neg_hip_norm;
			
			% Add to structure
			ID_metrics.(nameMetrics{i}).(conditions{t}).('H1') = work_H1 * dt;
			ID_metrics.(nameMetrics{i}).(conditions{t}).('H2') = work_neg_hip_norm;
			ID_metrics.(nameMetrics{i}).(conditions{t}).('H3') = work_H3 * dt;
			ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_WORK_NET') = net_work_hip;
			ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_WORK_POS') = work_pos_hip_norm;
			ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_WORK_NEG') = work_neg_hip_norm;
			ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER_POS') = power_pos_hip;
			ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER_NEG') = power_neg_hip;
			
			% KNEE
		elseif strcmp(nameMetrics{i}, 'knee_angle_r_moment')
			
			% Negative work
			K1 = ismember(spanLocsBelow, goodSpansBelow(1));
			work_K1 = trapz(ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER')(K1,1));
			
			% Check to see how many spans were detected and assign
			% appropriately.
			if length(goodSpansBelow) > 4
				K3 = ismember(spanLocsBelow, goodSpansBelow(3));
				work_K3 = trapz(ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER')(K3,1));
				K4 = ismember(spanLocsBelow, goodSpansBelow(end-1));
				K5 = ismember(spanLocsBelow, goodSpansBelow(end));
				work_K4 = trapz(ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER')(K4,1)) +...
					trapz(ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER')(K5,1));
				
			elseif length(goodSpansBelow) == 4
				K3 = ismember(spanLocsBelow, goodSpansBelow(3));
				work_K3 = trapz(ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER')(K3,1));
				K4 = ismember(spanLocsBelow, goodSpansBelow(4));
				work_K4 = trapz(ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER')(K4,1));
				
			elseif length(goodSpansBelow) <= 3
				K3 = ismember(spanLocsBelow, goodSpansBelow(end-1));
				work_K3 = trapz(ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER')(K3,1));
				K4 = ismember(spanLocsBelow, goodSpansBelow(end));
				work_K4 = trapz(ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER')(K4,1));
			end
			
			% Sum positive work to get total positive work
			work_neg_knee =  work_K1 + work_K3 + work_K4;
			work_neg_knee_norm = work_neg_knee*dt;
			
			% Multiply by 2 to get both limbs (assumed symmetry)
			work_neg_knee_bothLimbs = work_neg_knee_norm *2;
			power_neg_knee = work_neg_knee_bothLimbs / step_time_subject;
			
			% Positive work
			A = find(goodSpans == 2, 1);
			% If positive span was not detected (i.e., was less
			% than 5 frames)
			if isempty(A)
				goodSpans = find(spanLength>=1);   %get spans of 1 frame or greater
				positiveInSpans = ismember(spanLocs, goodSpans);  %indices of these spans
				ID_metrics.(nameMetrics{i}).(conditions{t}).('positive_power_index') = positiveInSpans;
			end
			
			K2 = ismember(spanLocs, goodSpans(2));
			K2_half = ismember(spanLocs, goodSpans(3));
			
			% Determine if there is more than one period above zero
			b = find(K2 > 0);
			bb = find(K2_half>0);
			bbb = find(K3>0);

			% If the assignment is in the range I'm looking for
			% (e.g., 15-40% of stride) then perform calcs
			if b(1) < 70
				
				% If there is then add them together
				if bb(1) < bbb(1)
					work_K2_half = trapz(ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER')(K2_half,1));
					work_K2_half2 =  trapz(ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER')(K2,1));
					work_K2 = work_K2_half + work_K2_half2;
					
				else
					work_K2 = trapz(ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER')(K2_half,1));
				end
				
				work_pos_knee_norm = work_K2 * dt;
			else
				work_pos_knee_norm = 0;
			end
			
			% Multiply by 2 to get both limbs (assumed symmetry)
			work_pos_knee_bothLimbs = work_pos_knee_norm  *2;
			power_pos_knee = work_pos_knee_bothLimbs / step_time_subject;
			
			% Net work
			net_work_knee = work_pos_knee_norm + work_neg_knee_norm;
			
			% Add to structure
			ID_metrics.(nameMetrics{i}).(conditions{t}).('K1') = work_K1 * dt;
			ID_metrics.(nameMetrics{i}).(conditions{t}).('K2') = work_K2 * dt;
			ID_metrics.(nameMetrics{i}).(conditions{t}).('K3') = work_K3 * dt;
			ID_metrics.(nameMetrics{i}).(conditions{t}).('K4') = work_K4 * dt;
			ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_WORK_NET') = net_work_knee;
			ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_WORK_POS') = work_pos_knee_norm;
			ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_WORK_NEG') = work_neg_knee_norm;
			ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER_POS') = power_pos_knee;
			ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER_NEG') = power_neg_knee;
			
			% ANKLE
		elseif strcmp(nameMetrics{i}, 'ankle_angle_r_moment')
			
			% Positive work
			A2 = ismember(spanLocs, goodSpans(end-1));
			c = find(A2 > 0);
			
			% Make sure it's the correct time stamp
			if c(1) < 40
				A2 = ismember(spanLocs, goodSpans(end));
			elseif c(end) > 70
				A2 = ismember(spanLocs, goodSpans(end-2));
			end
			
			work_A2 = trapz(ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER')(A2,1));
			
			% Sum positive work to get total positive work
			work_pos_ankle =  work_A2;
			work_pos_ankle_norm = work_pos_ankle*dt;
			
			% Multiply by 2 to get both limbs (assumed symmetry)
			work_pos_ankle_bothLimbs = work_pos_ankle_norm  *2;
			power_pos_ankle = work_pos_ankle_bothLimbs / step_time_subject;
			
			% Negative work - from 5-40% of gait cycle
			work_neg_ankle = trapz(ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER')(5:40,1));

			work_neg_ankle_norm = work_neg_ankle * dt;
			
			% Multiply by 2 to get both limbs (assumed symmetry)
			work_neg_ankle_bothLimbs = work_neg_ankle_norm  *2;
			power_neg_ankle = work_neg_ankle_bothLimbs / step_time_subject;
			
			% Net work
			net_work_ankle = work_pos_ankle_norm + work_neg_ankle_norm;
			
			% Add to structure
			ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_WORK_NET') = net_work_ankle;
			ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_WORK_POS') = work_pos_ankle_norm;
			ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_WORK_NEG') = work_neg_ankle_norm;
			ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER_POS') = power_pos_ankle;
			ID_metrics.(nameMetrics{i}).(conditions{t}).('JOINT_POWER_NEG') = power_neg_ankle;
			
		else
		end
	end
	g = g + 2;
end

for t = 1:length(conditions)
	% Positive work at all joints
	work_pos_total = ID_metrics.('hip_flexion_r_moment').(conditions{t}).('JOINT_WORK_POS') + ...
		ID_metrics.('knee_angle_r_moment').(conditions{t}).('JOINT_WORK_POS') + ...
		ID_metrics.('ankle_angle_r_moment').(conditions{t}).('JOINT_WORK_POS');
	% Negative work at all joints
	work_neg_total = ID_metrics.('hip_flexion_r_moment').(conditions{t}).('JOINT_WORK_NEG') + ...
		ID_metrics.('knee_angle_r_moment').(conditions{t}).('JOINT_WORK_NEG') + ...
		ID_metrics.('ankle_angle_r_moment').(conditions{t}).('JOINT_WORK_NEG');
	% Positive power at all joints
	power_pos_total = ID_metrics.('hip_flexion_r_moment').(conditions{t}).('JOINT_POWER_POS') + ...
		ID_metrics.('knee_angle_r_moment').(conditions{t}).('JOINT_POWER_POS') + ...
		ID_metrics.('ankle_angle_r_moment').(conditions{t}).('JOINT_POWER_POS');
	% Negative power at all joints
	power_neg_total = ID_metrics.('hip_flexion_r_moment').(conditions{t}).('JOINT_POWER_NEG') + ...
		ID_metrics.('knee_angle_r_moment').(conditions{t}).('JOINT_POWER_NEG') + ...
		ID_metrics.('ankle_angle_r_moment').(conditions{t}).('JOINT_POWER_NEG');
	
	% Determine percent contribution of each joint to total power
	powerP_perc_hip = (ID_metrics.('hip_flexion_r_moment').(conditions{t}).('JOINT_POWER_POS') / power_pos_total) * 100;
	powerP_perc_knee = (ID_metrics.('knee_angle_r_moment').(conditions{t}).('JOINT_POWER_POS') / power_pos_total) * 100;
	powerP_perc_ankle = (ID_metrics.('ankle_angle_r_moment').(conditions{t}).('JOINT_POWER_POS') / power_pos_total) * 100;
	
	powerN_perc_hip = (ID_metrics.('hip_flexion_r_moment').(conditions{t}).('JOINT_POWER_NEG') / power_neg_total) * 100;
	powerN_perc_knee = (ID_metrics.('knee_angle_r_moment').(conditions{t}).('JOINT_POWER_NEG') / power_neg_total) * 100;
	powerN_perc_ankle = (ID_metrics.('ankle_angle_r_moment').(conditions{t}).('JOINT_POWER_NEG') / power_neg_total) * 100;
	
	% Add to structure
	ID_metrics.('POWER_POS_TOTAL').(conditions{t}) = power_pos_total;
	ID_metrics.('POWER_NEG_TOTAL').(conditions{t}) = power_neg_total;
	ID_metrics.('WORK_POS_TOTAL').(conditions{t}) = work_pos_total;
	ID_metrics.('WORK_NEG_TOTAL').(conditions{t}) = work_neg_total;
	ID_metrics.('PERC_POS_HIP').(conditions{t}) = powerP_perc_hip;
	ID_metrics.('PERC_POS_KNEE').(conditions{t}) = powerP_perc_knee;
	ID_metrics.('PERC_POS_ANKLE').(conditions{t}) = powerP_perc_ankle;
	ID_metrics.('PERC_NEG_HIP').(conditions{t}) = powerN_perc_hip;
	ID_metrics.('PERC_NEG_KNEE').(conditions{t}) = powerN_perc_knee;
	ID_metrics.('PERC_NEG_ANKLE').(conditions{t}) = powerN_perc_ankle;
end
end

