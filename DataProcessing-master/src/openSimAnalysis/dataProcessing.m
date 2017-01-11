%% Script to Process data from Load sharing trials

% IK, ID, muscle analysis, and SO taken from Batch OpenSim Processing Scripts (BOPS)
% Copyright (C) 2015 Alice Mantoan, Monica Reggiani
% <https://simtk.org/home/bops>
% Please acknowledge the authors of BOPS when using this script.
clc; clear;

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath(genpath('GeneratePlots'));
addpath(genpath('dataAnalysis'));

%% Load the data
BasePath=uigetdir('../../../', 'Select Elaborated Data Folder');
cd(BasePath);
load('IDMetrics_all.mat'); load('IKMetrics_all.mat');

% Sort metrics so they are in the correct order
subjects = fieldnames(ID_metrics);
subjectNumber = zeros(length(subjects), 1);
newOrder = {'NA_slow', 'NA_fast', 'TBAS15_slow', 'TBAS15_fast', 'CRYE15_slow', 'CRYE15_fast',...
	'TYR15_slow', 'TYR15_fast', 'USMC15_slow', 'USMC15_fast', 'CORE15_slow', 'CORE15_fast',...
	'SORD15_slow', 'SORD15_fast', 'TBAS30_slow', 'TBAS30_fast', 'CRYE30_slow', 'CRYE30_fast',...
	'TYR30_slow', 'TYR30_fast', 'USMC30_slow', 'USMC30_fast', 'CORE30_slow', 'CORE30_fast',...
	'SORD30_slow', 'SORD30_fast'};

% Rearrange the names so they are in correct order
for t = 1:length(subjects)
	
	% Determine conditions for that subject
	subjectFields = fieldnames(ID_metrics.(subjects{t}));
	subjectName = subjects{t};
	Key   = '_';
	Index = strfind(subjectName, Key);
	Value = sscanf(subjectName(Index(1) + length(Key):end), '%g', 1);
	subjectNumber(t, 1) = Value;
	
	if length(subjectFields) ~= 26
		% Find which conditions are missing
		diffFields = setdiff(newOrder, subjectFields);
		
		% Remove those from new order variable
		for k = 1:length(diffFields)
			if k == 1
			subject_new_order = newOrder(cellfun('isempty', strfind(newOrder, diffFields(k))));
			else
				subject_new_order = subject_new_order(cellfun('isempty', strfind(subject_new_order, diffFields(k))));
			end
		end
		
		% Reorder conditions
		ID_metrics_ordered.(subjects{t}) = orderfields(ID_metrics.(subjects{t}), subject_new_order);
		IK_metrics_ordered.(subjects{t}) = orderfields(IK_metrics.(subjects{t}), subject_new_order);
	else
		ID_metrics_ordered.(subjects{t}) = orderfields(ID_metrics.(subjects{t}), newOrder);
		IK_metrics_ordered.(subjects{t}) = orderfields(IK_metrics.(subjects{t}), newOrder);
	end
end

% Determine how many condition fields exist
conditions = fieldnames(ID_metrics_ordered.(subjects{1}))';

% Create structure to store outputs
metricsID = fieldnames(ID_metrics_ordered.(subjects{1}).(conditions{1}));
metricsIK = fieldnames(IK_metrics_ordered.(subjects{1}).(conditions{1}));
metricsID(ismember(metricsID,{'POWER_POS_TOTAL', 'POWER_NEG_TOTAL',...
	'lumbar_bending_moment', 'lumbar_rotation_moment'}))=[];

%% Extract the variables of interest from the data
%[ID_metrics_for_analysis, IK_metrics_for_analysis] = extractVariables(subjects, metricsID, metricsIK, conditions, ID_metrics_ordered, IK_metrics_ordered);

[ID_metrics_for_analysis, IK_metrics_for_analysis] = extractVariables_02(subjects, metricsID, metricsIK, conditions, ID_metrics_ordered, IK_metrics_ordered);

% Quick loop through variables to convert zeros to NaNs
variable_names_ID = fieldnames(ID_metrics_for_analysis);
variable_names_IK = fieldnames(IK_metrics_for_analysis);

% ID data
for vName = 1:length(variable_names_ID)
	conditionsToNan = fieldnames(ID_metrics_for_analysis.(variable_names_ID{vName}));
	for c = 1:length(conditionsToNan)
		ID_metrics_for_analysis.(variable_names_ID{vName}).(conditionsToNan{c})...
			(~ID_metrics_for_analysis.(variable_names_ID{vName}).(conditionsToNan{c})) = nan;
	end
end

% IK data
for vName = 1:length(variable_names_IK)
	conditionsToNan = fieldnames(IK_metrics_for_analysis.(variable_names_IK{vName}));
	for c = 1:length(conditionsToNan)
		IK_metrics_for_analysis.(variable_names_IK{vName}).(conditionsToNan{c})...
			(~IK_metrics_for_analysis.(variable_names_IK{vName}).(conditionsToNan{c})) = nan;
	end
end

% Name condition labels
conditionLabels = {'TBAS'; 'cARM1'; 'cARM2'; 'pARM1'; 'pARM2'; 'pARM3'};

%% Moment plotting

% Moment trajectories analysis
% Obtain the mean and SD for the hip, knee, and ankle moments of
% interest
% Also cleans up the trajectories

combinedMetrics_ID = combineEnsembles(ID_metrics_for_analysis, subjects, conditions);

combinedMetrics_IK = combineEnsembles(IK_metrics_for_analysis, subjects, conditions);

% Redo joint powers for plots
ID_metrics_plots = workAndPowerCalcs(combinedMetrics_IK, combinedMetrics_ID, conditions, ID_metrics_for_analysis.step_time);

%% Plots - have to re-write this so it works for the code
plot_trajectories_diffTBAS(metrics_fromTBAS);

% plot trajectories
plot_trajectories(combinedMetrics_ID);

%% Save peaks in .mat files and .csv files

% Save moments first
% Make dir if it's not there
if ~isdir([BasePath, filesep, 'results_moments'])
	mkdir(BasePath, 'results_moments');
end

cd(fullfile(BasePath, 'results_moments'));
save('Moments.mat', '-struct', 'ID_metrics_for_analysis');

% Then kinematics
% Make dir if it's not there
if ~isdir([BasePath, filesep, 'results_kinematics'])
	mkdir(BasePath, 'results_kinematics');
end

cd(fullfile(BasePath, 'results_kinematics'));
save('kinematics.mat', '-struct', 'IK_metrics_for_analysis');

%% Create a table of data for export to R

% First put all conditions into one array
ID_metrics_final = combineConditionsPeaks(ID_metrics_for_analysis);
IK_metrics_final = combineConditionsPeaks(IK_metrics_for_analysis);

% Create paths to save data
savePathID = fullfile(BasePath, 'results_moments');
savePathIK = fullfile(BasePath, 'results_kinematics');

% Create tables and export csv for R
createTableCSV(ID_metrics_final, subjectNumber, conditions, savePathID);
createTableCSV(IK_metrics_final, subjectNumber, conditions, savePathIK);

%% Peak values of moments
% Obtain the mean and SD for the hip, knee, and ankle peak moments of
% interest. Data output is each column represents a different joint and
% each row represents a different armour condition

[peak_moments, peak_angles] = jointMomentsPeaksAnalysis(ID_metrics_final, IK_metrics_final);

%% Power percentages

% Obtain the mean and SD for the hip, knee, and ankle percent power
% contributions. Data output is each column represent a different joint and
% each row represents a different armour condition

power_percent = jointPowerPercentageAnalysis(ID_metrics_final);

% subplots for the joint power percentage analysis
plot_power_percentages(power_percent, conditionLabels)

%% Subplots for peak moments
plot_peak_moments(peak_moments, conditionLabels);
