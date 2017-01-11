function power_percent = jointPowerPercentageAnalysis(ID_metrics_for_analysis)
% Obtain the mean and SD for the hip, knee, and ankle percent power
% contributions.

%   INPUT - ID_metrics: structure containing the metrics from ID
 
%   OUTPUT - structure containing the joint power percentage metrics

% Obtain the mean and SD for the hip, knee, and ankle percent power
% contributions. Data output is each column represent a different joint and
% each row represents a different armour condition

power_percent = struct();

for t = 1:3
	if t == 1
		variable = 'HIP';
	elseif t == 2
		variable = 'KNEE';
	else
		variable = 'ANKLE';
	end	
	
	% Light = 15 kg, Heavy = 30kg
	% Slow = 5.5 km/h, Fast = 6.5 km/h
	
	% Positive power
	
	% Light and slow
	[power_percent.light_slow_mean_pos(:,t), power_percent.light_slow_stdev_pos(:,t)] = Mean_SD(ID_metrics_for_analysis.(['PERC_POS_', variable])(:, 1:2:13));
	% Light and fast
	[power_percent.light_fast_mean_pos(:,t), power_percent.light_fast_stdev_pos(:,t)] = Mean_SD(ID_metrics_for_analysis.(['PERC_POS_', variable])(:, 2:2:14));
	% Heavy and slow
	[power_percent.heavy_slow_mean_pos(:,t), power_percent.heavy_slow_stdev_pos(:,t)] = Mean_SD(ID_metrics_for_analysis.(['PERC_POS_', variable])(:, 15:2:25));
	% Heavy and fast
	[power_percent.heavy_fast_mean_pos(:,t), power_percent.heavy_fast_stdev_pos(:,t)] = Mean_SD(ID_metrics_for_analysis.(['PERC_POS_', variable])(:, 16:2:26));
	
	% Negative power
	
	% Light and slow
	[power_percent.light_slow_mean_neg(:,t), power_percent.light_slow_stdev_neg(:,t)] = Mean_SD(ID_metrics_for_analysis.(['PERC_NEG_', variable])(:, 1:2:13));
	% Light and fast
	[power_percent.light_fast_mean_neg(:,t), power_percent.light_fast_stdev_neg(:,t)] = Mean_SD(ID_metrics_for_analysis.(['PERC_NEG_', variable])(:, 2:2:14));
	% Heavy and slow
	[power_percent.heavy_slow_mean_neg(:,t), power_percent.heavy_slow_stdev_neg(:,t)] = Mean_SD(ID_metrics_for_analysis.(['PERC_NEG_', variable])(:, 15:2:25));
	% Heavy and fast
	[power_percent.heavy_fast_mean_neg(:,t), power_percent.heavy_fast_stdev_neg(:,t)] = Mean_SD(ID_metrics_for_analysis.(['PERC_NEG_', variable])(:, 16:2:26));
	
end

