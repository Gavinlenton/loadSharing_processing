function  combinedMetrics = combineEnsembles(metrics, subjects, conditions)
%Combine two ensemble averages into a mean and SD
%   INPUT - structure containing multiple ensembles (e.g., data from a gait
%           cycle).
%   OUTPUT - structure containing the ensemble mean and SD


% Determine which names of metrics to combine
metricsToCombine = fieldnames(metrics);

if any(strcmp(metricsToCombine, 'hip_flexion_r_moment'))
	
	% Selection of variables to keep for joint moments and powers
	metricsToCombine(~ismember(metricsToCombine,{'hip_flexion_r_moment', 'hip_adduction_r_moment', 'hip_rotation_r_moment',...
		'knee_angle_r_moment', 'knee_adduction_r_moment', 'knee_rotation_r_moment', 'ankle_angle_r_moment'}))=[];
	
elseif any(strcmp(metricsToCombine, 'hip_flexion_r'))
	
	% Selection of variables to keep for joint velocities and angles
	metricsToCombine(~ismember(metricsToCombine,{'hip_flexion_r', 'hip_adduction_r', 'hip_rotation_r', 'hip_flexion_vel', 'hip_adduction_vel',...
		'hip_rotation_vel', 'knee_angle_r', 'knee_flexion_vel', 'knee_adduction_r', 'knee_adduction_vel', 'knee_rotation_r',...
		'knee_rotation_vel', 'ankle_angle_r', 'ankle_flexion_vel', 'lumbar_extension'}))=[];
end

% Number of values in trial
lengthOfTrials = 101;
numParticipants = length(subjects);


% Combine data from all participants - this assumes that all trials are the same
% length
participants_number = lengthOfTrials*numParticipants;

% Loop through metrics of interest
for k = 1:length(metricsToCombine)
	
	% Loop through conditions
	for c = 1:length(conditions)
		
		if any(strcmp(metricsToCombine, 'hip_flexion_r_moment'))
			combinedMetrics.('toe_off_time').(conditions{c}) = metrics.step_time.(conditions{c})*0.6;
		end
		
		combinedMetrics.(metricsToCombine{k}).(conditions{c}) = metrics.(metricsToCombine{k}).(conditions{c})(2:end,:);
		
% 		% Take from row 4 to end because first few frames are dodgy
% 		averages = nansum(metrics.(metricsToCombine{k}).(conditions{c})(2:end,:), 2); % mean
% 		combined_averages = resample(averages, 101, length(averages), 2);
% 		
% 		if ~strcmp(metricsToCombine{k}(end-2:end), 'vel') % If it's velocity then there is no SD
% 			
% 			STDs = nansum(metrics.([metricsToCombine{k}, '_SD']).(conditions{c})(2:end,:), 2); % stdevs
% 			combined_STDs = resample(STDs, 101, length(STDs), 2);
% 			
% 			combined_variances = combined_STDs .^2; % variance
% 			
% 			totalsub = sum(participants_number ); % Number of ensembles to combine
% 			
% 			% Initialize variables
% 			i=1;
% 			combVar = zeros(size(combined_averages, 1), 1);
% 			combAvg = zeros(size(combined_averages, 1), 1);
% 			
% 			% Only combine when i is less than the length of the vector
% 			while i<=size(combined_averages ,1)
% 				% Combine to a single mean
% 				combAvg(i,1) = (lengthOfTrials * combined_averages(i,:)) / totalsub;
% 				
% 				if strcmp(metricsToCombine{k}, 'knee_angle_r_moment')
% 					combAvg(i,1) = combAvg(i,1) - 0.2;
% 				end
% 				
% 				% Combine to a single stdev
% 				combVar(i,1) = (totalsub*(combined_variances (i,:)*(participants_number -1)'+...
% 					combined_averages(i,:)^2 * participants_number')-(combined_averages (i,:)*...
% 					participants_number')^2)/(totalsub * totalsub-1)/numParticipants;
% 				i=i+1;
% 				
% 			end
% 			
% 			% Put into new structure
% 			combinedMetrics.(metricsToCombine{k}).(conditions{c}).mean = combAvg;
% 			combinedMetrics.(metricsToCombine{k}).(conditions{c}).var = combVar;
% 			combinedMetrics.(metricsToCombine{k}).(conditions{c}).sd = sqrt(combVar);
% 			
% 		else
% 			
% 			averages_participants = zeros(100,size(metrics.(metricsToCombine{k}).(conditions{c}), 2));
% 			% Loop through participants and get the SD for velocities
% 			for t = 1:size(averages_participants,2)
% 				
% 				averages_participants(:,t) = metrics.(metricsToCombine{k}).(conditions{c})(2:end,t);
% 			end
% 			
% 			STDs_ave = nanstd(averages_participants')';
% 			STDs2 = resample(STDs_ave, 101, length(STDs_ave), 2);
% 			
% 			% Initialize variables
% 			i=1;
% 			combAvg = zeros(size(combined_averages, 1), 1);
% 			
% 			% Only combine when i is less than the length of the vector
% 			while i<=size(combined_averages ,1)
% 				% Combine to a single mean
% 				combAvg(i,1) = (lengthOfTrials * combined_averages(i,:)) / totalsub;
% 				i=i+1;
% 				
% 			end
% 			
% 			% Put into new structure
% 			combinedMetrics.(metricsToCombine{k}).(conditions{c}).mean = combAvg;
% 			combinedMetrics.(metricsToCombine{k}).(conditions{c}).sd = STDs2;
			
		end
	end
end

