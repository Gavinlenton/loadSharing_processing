function [peak_moments, peak_angles] = jointMomentsPeaksAnalysis(ID_metrics_for_analysis, IK_metrics_for_analysis)
% Obtain the mean and SD for the hip, knee, and ankle percent power
% contributions.

%   INPUT - ID_metrics: structure containing the metrics from ID
%           IK_metrics: structure containing the metrics from IK

%   OUTPUT - structure containing the joint moment and joint angle peak metrics

% Obtain the mean and SD for the hip, knee, and ankle peak joint moments. 
% Data output is each column represent a different joint and each row 
% represents a different armour condition

peak_moments = struct();

for k = 1:2
	increment = 2;
	k_30 = k + 14;
	
	% Ankle
	[peak_moments.ank15_mean(:,k), peak_moments.ank15_SD(:, k)] = Mean_SD(abs(ID_metrics_for_analysis.ankle_plant_peak(:, (k+2):increment:(increment*6)+k)));
	[peak_moments.ank30_mean(:,k), peak_moments.ank30_SD(:, k)] = Mean_SD(abs(ID_metrics_for_analysis.ankle_plant_peak(:, k_30:increment:(increment*5)+k_30)));
	
	% Knee flexion
	[peak_moments.kneeFlex15_mean(:,k), peak_moments.kneeFlex15_SD(:, k)] = Mean_SD(ID_metrics_for_analysis.knee_flex_peak(:, k+2:increment:(increment*6)+k));
	[peak_moments.kneeFlex30_mean(:,k), peak_moments.kneeFlex30_SD(:, k)] = Mean_SD(ID_metrics_for_analysis.knee_flex_peak(:, k_30:increment:(increment*5)+k_30));
	
	% Knee extension
	[peak_moments.kneeExt15_mean(:,k), peak_moments.kneeExt15_SD(:, k)] = Mean_SD(ID_metrics_for_analysis.knee_ext_peak(:, k+2:increment:(increment*6)+k));
	[peak_moments.kneeExt30_mean(:,k), peak_moments.kneeExt30_SD(:, k)] = Mean_SD(ID_metrics_for_analysis.knee_ext_peak(:, k_30:increment:(increment*5)+k_30));
	
	% Knee adduction
	[peak_moments.kneeAdd15_mean(:,k), peak_moments.kneeAdd15_SD(:, k)] = Mean_SD(ID_metrics_for_analysis.knee_add_peak(:, k+2:increment:(increment*6)+k));
	[peak_moments.kneeAdd30_mean(:,k), peak_moments.kneeAdd30_SD(:, k)] = Mean_SD(ID_metrics_for_analysis.knee_add_peak(:, k_30:increment:(increment*5)+k_30));
	
	% Hip flexion
	[peak_moments.hipFlex15_mean(:,k), peak_moments.hipFlex15_SD(:, k)] = Mean_SD(abs(ID_metrics_for_analysis.hip_flex_peak(:, k+2:increment:(increment*6)+k)));
	[peak_moments.hipFlex30_mean(:,k), peak_moments.hipFlex30_SD(:, k)] = Mean_SD(abs(ID_metrics_for_analysis.hip_flex_peak(:, k_30:increment:(increment*5)+k_30)));
	
	% Hip extension
	[peak_moments.hipExt15_mean(:,k), peak_moments.hipExt15_SD(:, k)] = Mean_SD(abs(ID_metrics_for_analysis.hip_ext_peak(:, k+2:increment:(increment*6)+k)));
	[peak_moments.hipExt30_mean(:,k), peak_moments.hipExt30_SD(:, k)] = Mean_SD(abs(ID_metrics_for_analysis.hip_ext_peak(:, k_30:increment:(increment*5)+k_30)));
	
	% Hip abduction
	[peak_moments.hipAbd15_mean(:,k), peak_moments.hipAbd15_SD(:, k)] = Mean_SD(abs(ID_metrics_for_analysis.hip_abd_peak(:, k+2:increment:(increment*6)+k)));
	[peak_moments.hipAbd30_mean(:,k), peak_moments.hipAbd30_SD(:, k)] = Mean_SD(abs(ID_metrics_for_analysis.hip_abd_peak(:, k_30:increment:(increment*5)+k_30)));
	
	% Hip adduction
	[peak_moments.hipAdd15_mean(:,k), peak_moments.hipAdd15_SD(:, k)] = Mean_SD(abs(ID_metrics_for_analysis.hip_add_peak(:, k+2:increment:(increment*6)+k)));
	[peak_moments.hipAdd30_mean(:,k), peak_moments.hipAdd30_SD(:, k)] = Mean_SD(abs(ID_metrics_for_analysis.hip_add_peak(:, k_30:increment:(increment*5)+k_30)));
end

peak_angles = struct();

for tt = 1:2
	increment = 2;
	tt_30 = tt + 14;
	
	% Ankle
	[peak_angles.ank15_mean(:,tt), peak_angles.ank15_SD(:, tt)] = Mean_SD(abs(IK_metrics_for_analysis.ankle_plant_peak(:, (tt+2):increment:(increment*6)+tt)));
	[peak_angles.ank30_mean(:,tt), peak_angles.ank30_SD(:, tt)] = Mean_SD(abs(IK_metrics_for_analysis.ankle_plant_peak(:, tt_30:increment:(increment*5)+tt_30)));
	
	% Knee flexion
	[peak_angles.kneeFlex15_mean(:,tt), peak_angles.kneeFlex15_SD(:, tt)] = Mean_SD(IK_metrics_for_analysis.knee_flex_peak(:, tt+2:increment:(increment*6)+tt));
	[peak_angles.kneeFlex30_mean(:,tt), peak_angles.kneeFlex30_SD(:, tt)] = Mean_SD(IK_metrics_for_analysis.knee_flex_peak(:, tt_30:increment:(increment*5)+tt_30));
	
	% Knee extension
	[peak_angles.kneeExt15_mean(:,tt), peak_angles.kneeExt15_SD(:, tt)] = Mean_SD(IK_metrics_for_analysis.knee_ext_peak(:, tt+2:increment:(increment*6)+tt));
	[peak_angles.kneeExt30_mean(:,tt), peak_angles.kneeExt30_SD(:, tt)] = Mean_SD(IK_metrics_for_analysis.knee_ext_peak(:, tt_30:increment:(increment*5)+tt_30));
	
	% Knee adduction
	[peak_angles.kneeAdd15_mean(:,tt), peak_angles.kneeAdd15_SD(:, tt)] = Mean_SD(IK_metrics_for_analysis.knee_add_peak(:, tt+2:increment:(increment*6)+tt));
	[peak_angles.kneeAdd30_mean(:,tt), peak_angles.kneeAdd30_SD(:, tt)] = Mean_SD(IK_metrics_for_analysis.knee_add_peak(:, tt_30:increment:(increment*5)+tt_30));
	
	% Hip flexion
	[peak_angles.hipFlex15_mean(:,tt), peak_angles.hipFlex15_SD(:, tt)] = Mean_SD(abs(IK_metrics_for_analysis.hip_flex_peak(:, tt+2:increment:(increment*6)+tt)));
	[peak_angles.hipFlex30_mean(:,tt), peak_angles.hipFlex30_SD(:, tt)] = Mean_SD(abs(IK_metrics_for_analysis.hip_flex_peak(:, tt_30:increment:(increment*5)+tt_30)));
	
	% Hip extension
	[peak_angles.hipExt15_mean(:,tt), peak_angles.hipExt15_SD(:, tt)] = Mean_SD(abs(IK_metrics_for_analysis.hip_ext_peak(:, tt+2:increment:(increment*6)+tt)));
	[peak_angles.hipExt30_mean(:,tt), peak_angles.hipExt30_SD(:, tt)] = Mean_SD(abs(IK_metrics_for_analysis.hip_ext_peak(:, tt_30:increment:(increment*5)+tt_30)));
	
	% Hip abduction
	[peak_angles.hipAbd15_mean(:,tt), peak_angles.hipAbd15_SD(:, tt)] = Mean_SD(abs(IK_metrics_for_analysis.hip_abd_peak(:, tt+2:increment:(increment*6)+tt)));
	[peak_angles.hipAbd30_mean(:,tt), peak_angles.hipAbd30_SD(:, tt)] = Mean_SD(abs(IK_metrics_for_analysis.hip_abd_peak(:, tt_30:increment:(increment*5)+tt_30)));
	
	% Hip adduction
	[peak_angles.hipAdd15_mean(:,tt), peak_angles.hipAdd15_SD(:, tt)] = Mean_SD(abs(IK_metrics_for_analysis.hip_add_peak(:, tt+2:increment:(increment*6)+tt)));
	[peak_angles.hipAdd30_mean(:,tt), peak_angles.hipAdd30_SD(:, tt)] = Mean_SD(abs(IK_metrics_for_analysis.hip_add_peak(:, tt_30:increment:(increment*5)+tt_30)));
	
	% Lumbar extension
	[peak_angles.lumbExt15_mean(:,tt), peak_angles.lumbExt15_SD(:, tt)] = Mean_SD(abs(IK_metrics_for_analysis.lumbar_ext_peak(:, tt+2:increment:(increment*6)+tt)));
	[peak_angles.lumbExt30_mean(:,tt), peak_angles.lumbExt30_SD(:, tt)] = Mean_SD(abs(IK_metrics_for_analysis.lumbar_ext_peak(:, tt_30:increment:(increment*5)+tt_30)));
	
	% Lumbar flexion 
	[peak_angles.lumbFlex15_mean(:,tt), peak_angles.lumbFlex15_SD(:, tt)] = Mean_SD(abs(IK_metrics_for_analysis.lumbar_flex_peak(:, tt+2:increment:(increment*6)+tt)));
	[peak_angles.lumbFlex30_mean(:,tt), peak_angles.lumbFlex30_SD(:, tt)] = Mean_SD(abs(IK_metrics_for_analysis.lumbar_flex_peak(:, tt_30:increment:(increment*5)+tt_30)));
	
end

