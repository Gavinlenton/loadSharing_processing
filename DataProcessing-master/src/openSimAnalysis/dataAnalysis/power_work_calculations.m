function [metrics_with_powers] = power_work_calculations(ID_metrics, metrics_with_powers, nameMetrics, conditions, i, t, k, time, step_time_subject)
%Given an input structure with joint powers, compute the work and power
%parameters of interest

% Sampling time
dt = 1/100;
metric_name = [nameMetrics{i}(1:end-6), time];

abovePositiveThreshold = ID_metrics.(metric_name).(conditions{t})(:,k)>= 0;
spanLocs = bwlabel(abovePositiveThreshold);   %identify contiguous ones
spanLength = regionprops(abovePositiveThreshold, 'area');  %length of each span
spanLength = [ spanLength.Area];
goodSpans = find(spanLength>=5);   %get only spans of 5+ frames
positiveInSpans = ismember(spanLocs, goodSpans);  %indices of these spans
metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'positive_power_index']).(conditions{t})(:,k) = positiveInSpans;

% Negative periods of interest
belowPositiveThreshold = ID_metrics.(metric_name).(conditions{t})(:,k) < 0;
spanLocsBelow = bwlabel(belowPositiveThreshold);   %identify contiguous ones
spanLengthBelow = regionprops(belowPositiveThreshold, 'area');  %length of each span
spanLengthBelow = [ spanLengthBelow.Area];
goodSpansBelow = find(spanLengthBelow>=5);   %get only spans of 10+ frames
negativeInSpans = ismember(spanLocsBelow, goodSpansBelow);  %indices of these spans
metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'negative_power_index']).(conditions{t})(:,k) = negativeInSpans;

% POWER BURSTS - define burst periods as described in Winter
% 1987 text

if strcmp(time, 'JOINT_POWER_TOTAL')
	
	% HIP
	if strcmp(nameMetrics{i}, 'hip_flexion_r_moment')
		
		% Positive work
		if numel(goodSpans) == 3
			H1 = ismember(spanLocs, goodSpans(1));
			work_H1 = trapz(ID_metrics.(metric_name).(conditions{t})(H1,k));
			H3_half1 = ismember(spanLocs, goodSpans(2));
			work_H3_half1 = trapz(ID_metrics.(metric_name).(conditions{t})(H3_half1,k));
			H3_half2 = ismember(spanLocs, goodSpans(3));
			work_H3_half2 = trapz(ID_metrics.(metric_name).(conditions{t})(H3_half2,k));
			
			work_H3 = work_H3_half1 + work_H3_half2;
			% Sum positive work to get total positive work
			work_pos_hip =  work_H1 + work_H3;
			
		elseif numel(goodSpans) > 3
			
			H1 = ismember(spanLocs, goodSpans(1));
			work_H1 = trapz(ID_metrics.(metric_name).(conditions{t})(H1,k));
			H3_third1 = ismember(spanLocs, goodSpans(2));
			work_H3_third1 = trapz(ID_metrics.(metric_name).(conditions{t})(H3_third1,k));
			H3_third2 = ismember(spanLocs, goodSpans(3));
			work_H3_third2 = trapz(ID_metrics.(metric_name).(conditions{t})(H3_third2,k));
			H3_third3 = ismember(spanLocs, goodSpans(4));
			work_H3_third3 = trapz(ID_metrics.(metric_name).(conditions{t})(H3_third3,k));
			
			work_H3 = work_H3_third1 + work_H3_third2 + work_H3_third3;
			% Sum positive work to get total positive work
			work_pos_hip =  work_H1 + work_H3;
			
		else
			
			H1 = ismember(spanLocs, goodSpans(1));
			work_H1 = trapz(ID_metrics.(metric_name).(conditions{t})(H1,k));
			H3 = ismember(spanLocs, goodSpans(2));
			work_H3 = trapz(ID_metrics.(metric_name).(conditions{t})(H3,k));
			% Sum positive work to get total positive work
			work_pos_hip =  work_H1 + work_H3;
		end
		
		% Divide by sampling freq to get normalised work
		work_pos_hip_norm = work_pos_hip*dt;
		
		% Multiply by 2 to get both limbs (assumed symmetry)
		work_pos_hip_bothLimbs = work_pos_hip_norm *2;
		power_pos_hip = work_pos_hip_bothLimbs / step_time_subject;
		
		% Negative work
		if numel(goodSpansBelow) == 1
			H2 = ismember(spanLocsBelow, goodSpansBelow(1));
			work_neg_hip = trapz(ID_metrics.(metric_name).(conditions{t})(H2,k));
			
		elseif numel(goodSpansBelow) == 2
			H2_half1 = ismember(spanLocsBelow, goodSpansBelow(1));
			h2_half2 = ismember(spanLocsBelow, goodSpansBelow(2));
			
			work_h2half1 = trapz(ID_metrics.(metric_name).(conditions{t})(H2_half1,k));
			work_h2half2 = trapz(ID_metrics.(metric_name).(conditions{t})(h2_half2,k));
			work_neg_hip = work_h2half1 + work_h2half2;
			
		elseif numel(goodSpansBelow) > 2
			work_neg_hip = 0;
			for tt = 1:numel(goodSpansBelow)
				H2 = ismember(spanLocsBelow, goodSpansBelow(tt));
				work_H2 = trapz(ID_metrics.(metric_name).(conditions{t})(H2,k));
				work_neg_hip = work_neg_hip + work_H2;
			end
			
			
		end
		
		work_neg_hip_norm = work_neg_hip * dt;
		
		% Multiply by 2 to get both limbs (assumed symmetry)
		work_neg_hip_bothLimbs = work_neg_hip_norm *2;
		power_neg_hip = work_neg_hip_bothLimbs / step_time_subject;
		
		% Net work
		net_work_hip = work_pos_hip_norm + work_neg_hip_norm;
		
		% Add to structure
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'H1']).(conditions{t})(:,k) = work_H1 * dt;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'H2']).(conditions{t})(:,k)  = work_neg_hip_norm;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'H3']).(conditions{t})(:,k)  = work_H3 * dt;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_NET']).(conditions{t})(:,k)  = net_work_hip;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_POS']).(conditions{t})(:,k)  = work_pos_hip_norm;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_NEG']).(conditions{t})(:,k)  = work_neg_hip_norm;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k)  = power_pos_hip;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k)  = power_neg_hip;
		
		% KNEE
	elseif strcmp(nameMetrics{i}, 'knee_angle_r_moment')
		
		% Negative work
		work_neg_knee = 0;
		for tt = 1:numel(goodSpansBelow)
			K1 = ismember(spanLocsBelow, goodSpansBelow(tt));
			work_K1 = trapz(ID_metrics.(metric_name).(conditions{t})(K1,k));
			work_neg_knee = work_neg_knee + work_K1;
		end
		
		% Sum negative work to get total positive work
		work_neg_knee_norm = work_neg_knee*dt;
		
		% Multiply by 2 to get both limbs (assumed symmetry)
		work_neg_knee_bothLimbs = work_neg_knee_norm *2;
		power_neg_knee = work_neg_knee_bothLimbs / step_time_subject;
		
		% Positive work
		goodSpans = find(spanLength>=1);   %get spans of 1 frame or greater
		positiveInSpans = ismember(spanLocs, goodSpans);  %indices of these spans
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'positive_power_index']).(conditions{t})(:,k) = positiveInSpans;
		
		work_pos_knee = 0;
		for tt = 1:numel(goodSpans)
			K2 = ismember(spanLocs, goodSpans(tt));
			work_K2 = trapz(ID_metrics.(metric_name).(conditions{t})(K2,k));
			work_pos_knee = work_pos_knee + work_K2;
		end
		
		work_pos_knee_norm = work_pos_knee * dt;
		
		% Multiply by 2 to get both limbs (assumed symmetry)
		work_pos_knee_bothLimbs = work_pos_knee_norm  *2;
		power_pos_knee = work_pos_knee_bothLimbs / step_time_subject;
		
		% Net work
		net_work_knee = work_pos_knee_norm + work_neg_knee_norm;
		
		% Add to structure
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_NET']).(conditions{t})(:,k)  = net_work_knee;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_POS']).(conditions{t})(:,k)  = work_pos_knee_norm;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_NEG']).(conditions{t})(:,k)  = work_neg_knee_norm;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k)  = power_pos_knee;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k)  = power_neg_knee;
		
		% ANKLE
	elseif strcmp(nameMetrics{i}, 'ankle_angle_r_moment')
		
		% Positive work
		work_pos_ankle = 0;
		for tt = 1:numel(goodSpans)
			A2 = ismember(spanLocs, goodSpans(tt));
			work_A2 = trapz(ID_metrics.(metric_name).(conditions{t})(A2,k));
			% Sum positive work to get total positive work
			work_pos_ankle = work_pos_ankle + work_A2;
		end
		
		% Normalise to collection frequency
		work_pos_ankle_norm = work_pos_ankle*dt;
		
		% Multiply by 2 to get both limbs (assumed symmetry)
		work_pos_ankle_bothLimbs = work_pos_ankle_norm  *2;
		power_pos_ankle = work_pos_ankle_bothLimbs / step_time_subject;
		
		% Negative work - from 5-40% of gait cycle
		work_neg_ankle = trapz(ID_metrics.(metric_name).(conditions{t})(5:40,k));
		
		work_neg_ankle_norm = work_neg_ankle * dt;
		
		% Multiply by 2 to get both limbs (assumed symmetry)
		work_neg_ankle_bothLimbs = work_neg_ankle_norm  *2;
		power_neg_ankle = work_neg_ankle_bothLimbs / step_time_subject;
		
		% Net work
		net_work_ankle = work_pos_ankle_norm + work_neg_ankle_norm;
		
		% Add to structure
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_NET']).(conditions{t})(:,k)  = net_work_ankle;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_POS']).(conditions{t})(:,k)  = work_pos_ankle_norm;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_NEG']).(conditions{t})(:,k)  = work_neg_ankle_norm;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k)  = power_pos_ankle;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k)  = power_neg_ankle;
		
		%% Get power percentages
		
		% Positive work at all joints
		work_pos_total(:,k) = metrics_with_powers.(time).([nameMetrics{2}(1:end-6), 'JOINT_WORK_POS']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{5}(1:end-6), 'JOINT_WORK_POS']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{8}(1:end-6), 'JOINT_WORK_POS']).(conditions{t})(:,k);
		% Negative work at all joints
		work_neg_total(:,k) = metrics_with_powers.(time).([nameMetrics{2}(1:end-6), 'JOINT_WORK_NEG']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{5}(1:end-6), 'JOINT_WORK_NEG']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{8}(1:end-6), 'JOINT_WORK_NEG']).(conditions{t})(:,k);
		% Positive power at all joints
		power_pos_total(:,k) = metrics_with_powers.(time).([nameMetrics{2}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{5}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{8}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k);
		% Negative power at all joints
		power_neg_total(:,k) = metrics_with_powers.(time).([nameMetrics{2}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{5}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{8}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k);
		
		% Determine percent contribution of each joint to total power
		powerP_perc_hip(:,k) = (metrics_with_powers.(time).([nameMetrics{2}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k) / power_pos_total(:,k)) * 100;
		powerP_perc_knee(:,k) = (metrics_with_powers.(time).([nameMetrics{5}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k) / power_pos_total(:,k)) * 100;
		powerP_perc_ankle(:,k) = (metrics_with_powers.(time).([nameMetrics{8}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k) / power_pos_total(:,k)) * 100;
		
		powerN_perc_hip(:,k) = (metrics_with_powers.(time).([nameMetrics{2}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k) / power_neg_total(:,k)) * 100;
		powerN_perc_knee(:,k) = (metrics_with_powers.(time).([nameMetrics{5}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k) / power_neg_total(:,k)) * 100;
		powerN_perc_ankle(:,k) = (metrics_with_powers.(time).([nameMetrics{8}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k) / power_neg_total(:,k)) * 100;
		
		% Add to structure
		metrics_with_powers.(time).('WORK_POS_TOTAL').(conditions{t})(:,k) = work_pos_total(:,k);
		metrics_with_powers.(time).('WORK_NEG_TOTAL').(conditions{t})(:,k)= work_neg_total(:,k);
		metrics_with_powers.(time).('POWER_POS_TOTAL').(conditions{t})(:,k) = power_pos_total(:,k);
		metrics_with_powers.(time).('POWER_NEG_TOTAL').(conditions{t})(:,k) = power_neg_total(:,k);
		metrics_with_powers.(time).('PERC_POS_HIP').(conditions{t})(:,k) = powerP_perc_hip(:,k);
		metrics_with_powers.(time).('PERC_POS_KNEE').(conditions{t})(:,k) = powerP_perc_knee(:,k);
		metrics_with_powers.(time).('PERC_POS_ANKLE').(conditions{t})(:,k) = powerP_perc_ankle(:,k);
		metrics_with_powers.(time).('PERC_NEG_HIP').(conditions{t})(:,k) = powerN_perc_hip(:,k);
		metrics_with_powers.(time).('PERC_NEG_KNEE').(conditions{t})(:,k) = powerN_perc_knee(:,k);
		metrics_with_powers.(time).('PERC_NEG_ANKLE').(conditions{t})(:,k) = powerN_perc_ankle(:,k);
		
	else
	end
	
elseif strcmp(time, 'JOINT_POWER_STRIDE')
	% HIP
	if strcmp(nameMetrics{i}, 'hip_flexion_r_moment')
		
		% Positive work
		work_pos_hip = 0;
		for tt = 1:numel(goodSpans)
			H1 = ismember(spanLocs, goodSpans(tt));
			work_H1 = trapz(ID_metrics.(metric_name).(conditions{t})(H1,k));
			work_pos_hip = work_pos_hip + work_H1;
		end
		
		% Divide by sampling freq to get normalised work
		work_pos_hip_norm = work_pos_hip*dt;
		
		% Multiply by 2 to get both limbs (assumed symmetry)
		work_pos_hip_bothLimbs = work_pos_hip_norm *2;
		power_pos_hip = work_pos_hip_bothLimbs / step_time_subject;
		
		% Negative work
		work_neg_hip = 0;
		for tt = 1:numel(goodSpansBelow)
			H2 = ismember(spanLocsBelow, goodSpansBelow(tt));
			work_H2 = trapz(ID_metrics.(metric_name).(conditions{t})(H2,k));
			work_neg_hip = work_neg_hip + work_H2;
		end
		work_neg_hip_norm = work_neg_hip * dt;
		
		% Multiply by 2 to get both limbs (assumed symmetry)
		work_neg_hip_bothLimbs = work_neg_hip_norm *2;
		power_neg_hip = work_neg_hip_bothLimbs / step_time_subject;
		
		% Net work
		net_work_hip = work_pos_hip_norm + work_neg_hip_norm;
		
		% Add to structure
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_NET']).(conditions{t})(:,k)  = net_work_hip;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_POS']).(conditions{t})(:,k)  = work_pos_hip_norm;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_NEG']).(conditions{t})(:,k)  = work_neg_hip_norm;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k)  = power_pos_hip;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k)  = power_neg_hip;
		
	elseif strcmp(nameMetrics{i}, 'knee_angle_r_moment')
		% KNEE
		% Negative work
		work_neg_knee = 0;
		for tt = 1:numel(goodSpansBelow)
			K1 = ismember(spanLocsBelow, goodSpansBelow(tt));
			work_K1 = trapz(ID_metrics.(metric_name).(conditions{t})(K1,k));
			work_neg_knee = work_neg_knee + work_K1;
		end
		
		% Sum negative work to get total positive work
		work_neg_knee_norm = work_neg_knee*dt;
		
		% Multiply by 2 to get both limbs (assumed symmetry)
		work_neg_knee_bothLimbs = work_neg_knee_norm *2;
		power_neg_knee = work_neg_knee_bothLimbs / step_time_subject;
		
		% Positive work
		goodSpans = find(spanLength>=1);   %get spans of 1 frame or greater
		positiveInSpans = ismember(spanLocs, goodSpans);  %indices of these spans
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'positive_power_index']).(conditions{t})(:,k) = positiveInSpans;
		
		work_pos_knee = 0;
		for tt = 1:numel(goodSpans)
			K2 = ismember(spanLocs, goodSpans(tt));
			work_K2 = trapz(ID_metrics.(metric_name).(conditions{t})(K2,k));
			work_pos_knee = work_pos_knee + work_K2;
		end
		
		work_pos_knee_norm = work_pos_knee * dt;
		
		% Multiply by 2 to get both limbs (assumed symmetry)
		work_pos_knee_bothLimbs = work_pos_knee_norm  *2;
		power_pos_knee = work_pos_knee_bothLimbs / step_time_subject;
		
		% Net work
		net_work_knee = work_pos_knee_norm + work_neg_knee_norm;
		
		% Add to structure
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_NET']).(conditions{t})(:,k)  = net_work_knee;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_POS']).(conditions{t})(:,k)  = work_pos_knee_norm;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_NEG']).(conditions{t})(:,k)  = work_neg_knee_norm;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k)  = power_pos_knee;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k)  = power_neg_knee;
		
	elseif strcmp(nameMetrics{i}, 'ankle_angle_r_moment')
		% ANKLE
		work_pos_ankle = 0;
		for tt = 1:numel(goodSpans)
			A2 = ismember(spanLocs, goodSpans(tt));
			work_A2 = trapz(ID_metrics.(metric_name).(conditions{t})(A2,k));
			% Sum positive work to get total positive work
			work_pos_ankle = work_pos_ankle + work_A2;
		end
		
		% Normalise to collection frequency
		work_pos_ankle_norm = work_pos_ankle*dt;
		
		% Multiply by 2 to get both limbs (assumed symmetry)
		work_pos_ankle_bothLimbs = work_pos_ankle_norm  *2;
		power_pos_ankle = work_pos_ankle_bothLimbs / step_time_subject;
		
		% Negative work 
		work_neg_ankle = 0;
		for tt = 1:numel(goodSpansBelow)
			A1 = ismember(spanLocsBelow, goodSpansBelow(tt));
			work_A1 = trapz(ID_metrics.(metric_name).(conditions{t})(A1,k));
			% Sum positive work to get total positive work
			work_neg_ankle = work_neg_ankle + work_A1;
		end
		
		work_neg_ankle_norm = work_neg_ankle * dt;
		
		% Multiply by 2 to get both limbs (assumed symmetry)
		work_neg_ankle_bothLimbs = work_neg_ankle_norm  *2;
		power_neg_ankle = work_neg_ankle_bothLimbs / step_time_subject;
		
		% Net work
		net_work_ankle = work_pos_ankle_norm + work_neg_ankle_norm;
		
		% Add to structure
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_NET']).(conditions{t})(:,k)  = net_work_ankle;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_POS']).(conditions{t})(:,k)  = work_pos_ankle_norm;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_NEG']).(conditions{t})(:,k)  = work_neg_ankle_norm;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k)  = power_pos_ankle;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k)  = power_neg_ankle;
		
		%% Get power percentages
		
		% Positive work at all joints
		work_pos_total(:,k) = metrics_with_powers.(time).([nameMetrics{2}(1:end-6), 'JOINT_WORK_POS']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{5}(1:end-6), 'JOINT_WORK_POS']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{8}(1:end-6), 'JOINT_WORK_POS']).(conditions{t})(:,k);
		% Negative work at all joints
		work_neg_total(:,k) = metrics_with_powers.(time).([nameMetrics{2}(1:end-6), 'JOINT_WORK_NEG']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{5}(1:end-6), 'JOINT_WORK_NEG']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{8}(1:end-6), 'JOINT_WORK_NEG']).(conditions{t})(:,k);
		% Positive power at all joints
		power_pos_total(:,k) = metrics_with_powers.(time).([nameMetrics{2}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{5}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{8}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k);
		% Negative power at all joints
		power_neg_total(:,k) = metrics_with_powers.(time).([nameMetrics{2}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{5}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{8}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k);
		
		% Determine percent contribution of each joint to total power
		powerP_perc_hip(:,k) = (metrics_with_powers.(time).([nameMetrics{2}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k) / power_pos_total(:,k)) * 100;
		powerP_perc_knee(:,k) = (metrics_with_powers.(time).([nameMetrics{5}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k) / power_pos_total(:,k)) * 100;
		powerP_perc_ankle(:,k) = (metrics_with_powers.(time).([nameMetrics{8}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k) / power_pos_total(:,k)) * 100;
		
		powerN_perc_hip(:,k) = (metrics_with_powers.(time).([nameMetrics{2}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k) / power_neg_total(:,k)) * 100;
		powerN_perc_knee(:,k) = (metrics_with_powers.(time).([nameMetrics{5}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k) / power_neg_total(:,k)) * 100;
		powerN_perc_ankle(:,k) = (metrics_with_powers.(time).([nameMetrics{8}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k) / power_neg_total(:,k)) * 100;
		
		% Add to structure
		metrics_with_powers.(time).('WORK_POS_TOTAL').(conditions{t})(:,k) = work_pos_total(:,k);
		metrics_with_powers.(time).('WORK_NEG_TOTAL').(conditions{t})(:,k)= work_neg_total(:,k);
		metrics_with_powers.(time).('POWER_POS_TOTAL').(conditions{t})(:,k) = power_pos_total(:,k);
		metrics_with_powers.(time).('POWER_NEG_TOTAL').(conditions{t})(:,k) = power_neg_total(:,k);
		metrics_with_powers.(time).('PERC_POS_HIP').(conditions{t})(:,k) = powerP_perc_hip(:,k);
		metrics_with_powers.(time).('PERC_POS_KNEE').(conditions{t})(:,k) = powerP_perc_knee(:,k);
		metrics_with_powers.(time).('PERC_POS_ANKLE').(conditions{t})(:,k) = powerP_perc_ankle(:,k);
		metrics_with_powers.(time).('PERC_NEG_HIP').(conditions{t})(:,k) = powerN_perc_hip(:,k);
		metrics_with_powers.(time).('PERC_NEG_KNEE').(conditions{t})(:,k) = powerN_perc_knee(:,k);
		metrics_with_powers.(time).('PERC_NEG_ANKLE').(conditions{t})(:,k) = powerN_perc_ankle(:,k);
	end
	
else
	
	% HIP
	if strcmp(nameMetrics{i}, 'hip_flexion_r_moment')
		
		% Positive work
		work_pos_hip = 0;
		for tt = 1:numel(goodSpans)
			H1 = ismember(spanLocs, goodSpans(tt));
			work_H1 = trapz(ID_metrics.(metric_name).(conditions{t})(H1,k));
			work_pos_hip = work_pos_hip + work_H1;
		end
		
		% Divide by sampling freq to get normalised work
		work_pos_hip_norm = work_pos_hip*dt;
		
		% Multiply by 2 to get both limbs (assumed symmetry)
		work_pos_hip_bothLimbs = work_pos_hip_norm *2;
		power_pos_hip = work_pos_hip_bothLimbs / step_time_subject;
		
		% Negative work
		work_neg_hip = 0;
		for tt = 1:numel(goodSpansBelow)
			H2 = ismember(spanLocsBelow, goodSpansBelow(tt));
			work_H2 = trapz(ID_metrics.(metric_name).(conditions{t})(H2,k));
			work_neg_hip = work_neg_hip + work_H2;
		end
		work_neg_hip_norm = work_neg_hip * dt;
		
		% Multiply by 2 to get both limbs (assumed symmetry)
		work_neg_hip_bothLimbs = work_neg_hip_norm *2;
		power_neg_hip = work_neg_hip_bothLimbs / step_time_subject;
		
		% Net work
		net_work_hip = work_pos_hip_norm + work_neg_hip_norm;
		
		% Add to structure
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_NET']).(conditions{t})(:,k)  = net_work_hip;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_POS']).(conditions{t})(:,k)  = work_pos_hip_norm;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_NEG']).(conditions{t})(:,k)  = work_neg_hip_norm;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k)  = power_pos_hip;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k)  = power_neg_hip;
		
	elseif strcmp(nameMetrics{i}, 'knee_angle_r_moment')
		% KNEE
		% Negative work
		work_neg_knee = 0;
		for tt = 1:numel(goodSpansBelow)
			K1 = ismember(spanLocsBelow, goodSpansBelow(tt));
			work_K1 = trapz(ID_metrics.(metric_name).(conditions{t})(K1,k));
			if tt == 1
				K1_real = work_K1;
				length_span = find(K1 >0);
			end
			if length(length_span) < 10 && tt == 2
				K1_real = K1_real+work_K1;
			end
			work_neg_knee = work_neg_knee + work_K1;
		end
		
		% Sum negative work to get total positive work
		work_neg_knee_norm = work_neg_knee*dt;
		K1_norm = K1_real*dt;
		
		% Multiply by 2 to get both limbs (assumed symmetry)
		work_neg_knee_bothLimbs = work_neg_knee_norm *2;
		K1_norm_bothLimbs = K1_norm *2;
		power_neg_knee = work_neg_knee_bothLimbs / step_time_subject;
		power_K1 = K1_norm_bothLimbs / step_time_subject;

		% Positive work
		goodSpans = find(spanLength>=1);   %get spans of 1 frame or greater
		positiveInSpans = ismember(spanLocs, goodSpans);  %indices of these spans
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'positive_power_index']).(conditions{t})(:,k) = positiveInSpans;
		
		work_pos_knee = 0;
		for tt = 1:numel(goodSpans)
			K2 = ismember(spanLocs, goodSpans(tt));
			work_K2 = trapz(ID_metrics.(metric_name).(conditions{t})(K2,k));
			work_pos_knee = work_pos_knee + work_K2;
		end
		
		work_pos_knee_norm = work_pos_knee * dt;
		
		% Multiply by 2 to get both limbs (assumed symmetry)
		work_pos_knee_bothLimbs = work_pos_knee_norm  *2;
		power_pos_knee = work_pos_knee_bothLimbs / step_time_subject;
		
		% Net work
		net_work_knee = work_pos_knee_norm + work_neg_knee_norm;
		
		% Add to structure
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_NET']).(conditions{t})(:,k)  = net_work_knee;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'WORK_K1']).(conditions{t})(:,k)  = power_K1;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_POS']).(conditions{t})(:,k)  = work_pos_knee_norm;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_NEG']).(conditions{t})(:,k)  = work_neg_knee_norm;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k)  = power_pos_knee;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k)  = power_neg_knee;
		
	elseif strcmp(nameMetrics{i}, 'ankle_angle_r_moment')
		% ANKLE
		work_pos_ankle = 0;
		for tt = 1:numel(goodSpans)
			A2 = ismember(spanLocs, goodSpans(tt));
			work_A2 = trapz(ID_metrics.(metric_name).(conditions{t})(A2,k));
			% Sum positive work to get total positive work
			work_pos_ankle = work_pos_ankle + work_A2;
		end
		
		% Normalise to collection frequency
		work_pos_ankle_norm = work_pos_ankle*dt;
		
		% Multiply by 2 to get both limbs (assumed symmetry)
		work_pos_ankle_bothLimbs = work_pos_ankle_norm  *2;
		power_pos_ankle = work_pos_ankle_bothLimbs / step_time_subject;
		
		% Negative work - from 5-40% of gait cycle
		work_neg_ankle = 0;
		for tt = 1:numel(goodSpansBelow)
			A1 = ismember(spanLocsBelow, goodSpansBelow(tt));
			work_A1 = trapz(ID_metrics.(metric_name).(conditions{t})(A1,k));
			% Sum positive work to get total positive work
			work_neg_ankle = work_neg_ankle + work_A1;
		end
		
		work_neg_ankle_norm = work_neg_ankle * dt;
		
		% Multiply by 2 to get both limbs (assumed symmetry)
		work_neg_ankle_bothLimbs = work_neg_ankle_norm  *2;
		power_neg_ankle = work_neg_ankle_bothLimbs / step_time_subject;
		
		% Net work
		net_work_ankle = work_pos_ankle_norm + work_neg_ankle_norm;
		
		% Add to structure
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_NET']).(conditions{t})(:,k)  = net_work_ankle;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_POS']).(conditions{t})(:,k)  = work_pos_ankle_norm;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_WORK_NEG']).(conditions{t})(:,k)  = work_neg_ankle_norm;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k)  = power_pos_ankle;
		metrics_with_powers.(time).([nameMetrics{i}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k)  = power_neg_ankle;
		
		%% Get power percentages
		
		% Positive work at all joints
		work_pos_total(:,k) = metrics_with_powers.(time).([nameMetrics{2}(1:end-6), 'JOINT_WORK_POS']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{5}(1:end-6), 'JOINT_WORK_POS']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{8}(1:end-6), 'JOINT_WORK_POS']).(conditions{t})(:,k);
		% Negative work at all joints
		work_neg_total(:,k) = metrics_with_powers.(time).([nameMetrics{2}(1:end-6), 'JOINT_WORK_NEG']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{5}(1:end-6), 'JOINT_WORK_NEG']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{8}(1:end-6), 'JOINT_WORK_NEG']).(conditions{t})(:,k);
		% Positive power at all joints
		power_pos_total(:,k) = metrics_with_powers.(time).([nameMetrics{2}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{5}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{8}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k);
		% Negative power at all joints
		power_neg_total(:,k) = metrics_with_powers.(time).([nameMetrics{2}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{5}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k) + ...
			metrics_with_powers.(time).([nameMetrics{8}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k);
		
		% Determine percent contribution of each joint to total power
		powerP_perc_hip(:,k) = (metrics_with_powers.(time).([nameMetrics{2}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k) / power_pos_total(:,k)) * 100;
		powerP_perc_knee(:,k) = (metrics_with_powers.(time).([nameMetrics{5}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k) / power_pos_total(:,k)) * 100;
		powerP_perc_ankle(:,k) = (metrics_with_powers.(time).([nameMetrics{8}(1:end-6), 'JOINT_POWER_POS']).(conditions{t})(:,k) / power_pos_total(:,k)) * 100;
		
		powerN_perc_hip(:,k) = (metrics_with_powers.(time).([nameMetrics{2}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k) / power_neg_total(:,k)) * 100;
		powerN_perc_knee(:,k) = (metrics_with_powers.(time).([nameMetrics{5}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k) / power_neg_total(:,k)) * 100;
		powerN_perc_ankle(:,k) = (metrics_with_powers.(time).([nameMetrics{8}(1:end-6), 'JOINT_POWER_NEG']).(conditions{t})(:,k) / power_neg_total(:,k)) * 100;
		
		% Add to structure
		metrics_with_powers.(time).('WORK_POS_TOTAL').(conditions{t})(:,k) = work_pos_total(:,k);
		metrics_with_powers.(time).('WORK_NEG_TOTAL').(conditions{t})(:,k)= work_neg_total(:,k);
		metrics_with_powers.(time).('POWER_POS_TOTAL').(conditions{t})(:,k) = power_pos_total(:,k);
		metrics_with_powers.(time).('POWER_NEG_TOTAL').(conditions{t})(:,k) = power_neg_total(:,k);
		metrics_with_powers.(time).('PERC_POS_HIP').(conditions{t})(:,k) = powerP_perc_hip(:,k);
		metrics_with_powers.(time).('PERC_POS_KNEE').(conditions{t})(:,k) = powerP_perc_knee(:,k);
		metrics_with_powers.(time).('PERC_POS_ANKLE').(conditions{t})(:,k) = powerP_perc_ankle(:,k);
		metrics_with_powers.(time).('PERC_NEG_HIP').(conditions{t})(:,k) = powerN_perc_hip(:,k);
		metrics_with_powers.(time).('PERC_NEG_KNEE').(conditions{t})(:,k) = powerN_perc_knee(:,k);
		metrics_with_powers.(time).('PERC_NEG_ANKLE').(conditions{t})(:,k) = powerN_perc_ankle(:,k);
		
	end
end
