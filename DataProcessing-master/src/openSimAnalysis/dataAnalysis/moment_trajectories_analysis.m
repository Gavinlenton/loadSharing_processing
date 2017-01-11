function moments_trajectories = moment_trajectories_analysis(ID_metrics, subjects)
% Extracts the metrics of interest from the load sharing data for all participants
% and puts them into structures for further analysis

%   INPUT: ID_metrics - structure containing the metrics from ID

%   OUTPUT - structure containing the metrics for further analysis

moments_trajectories = struct();

% Combine the ensemble averages into a means and SDs
combineEnsembles(ID_metrics, subjects);


% Loop through speeds
for m = 1:2
	increment = 2;
	m_30 = m + 14;
	
	% Ankle
	[ank15_mean(:, :, m), ank15_SD(:, :, m), ank15Upper(:, :, m), ank15Lower(:, :, m)]...
		= Mean_SD_3D(ID_metrics.ankle_angle_r_moment(:, :, (m+2):increment:(increment*6)+m));
	[ank30_mean(:, :, m), ank30_SD(:, :, m), ank30Upper(:, :, m), ank30Lower(:, :, m)]...
		= Mean_SD_3D(ID_metrics.ankle_angle_r_moment(:, :, m_30:increment:(increment*5)+m_30));
	
	% Knee
	[knee15_mean(:, :, m), knee15_SD(:, :, m), knee15Upper(:,:,m), knee15Lower(:,:,m)]...
		= Mean_SD_3D(ID_metrics.knee_angle_r_moment(:, :,  m+2:increment:(increment*6)+m));
	[knee30_mean(:, :, m), knee30_SD(:, :, m), knee30Upper(:,:,m), knee30Lower(:,:,m)]...
		= Mean_SD_3D(ID_metrics.knee_angle_r_moment(:, :, m_30:increment:(increment*5)+m_30));
	
	% Hip
	[hip15_mean(:, :, m), hip15_SD(:, :, m), hip15Upper(:,:,m), hip15Lower(:,:,m)]...
		= Mean_SD_3D(ID_metrics.hip_flexion_r_moment(:, :, m+2:increment:(increment*6)+m));
	[hip30_mean(:, :, m), hip30_SD(:, :, m), hip30Upper(:,:,m), hip30Lower(:,:,m)]...
		= Mean_SD_3D(ID_metrics.hip_flexion_r_moment(:, :, m_30:increment:(increment*5)+m_30));
end

end

