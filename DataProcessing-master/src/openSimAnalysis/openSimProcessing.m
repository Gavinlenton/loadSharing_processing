%% Script to run IK, ID, and any further analysis on load sharing data %%

% IK, ID, muscle analysis, and SO taken from Batch OpenSim Processing Scripts (BOPS)
% Copyright (C) 2015 Alice Mantoan, Monica Reggiani
% <https://simtk.org/home/bops>
% Please acknowledge the authors of BOPS when using this script.
clc; clear;

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath(genpath('..\'));
addpath(genpath('C:\Users\s2921887\Desktop\LS_main_data_collection\model_scaling'));
addpath(genpath('C:\Users\s2921887\Desktop\LS_main_data_collection\DataProcessing-master'));

%% Define BasePath with dynamicElaboration outputs and BOPS folder

% Select dynElab folder
BasePath=uigetdir('..\..\..\', 'Select Elaborated Data Folder');
% Select BOPS-master folder
folderBOPS = uigetdir('.\', 'Select the BOPS processing folder');

% Define subject names
subsdir=dir(BasePath);
for xx=1:length(subsdir)
	isub = [subsdir(:).isdir]; %# returns logical vector of subdirectories
	subjectNames = {subsdir(isub).name}';
	subjectNames(ismember(subjectNames,{'.','..'})) = [];
end

%% Loop through subjects
for nS = 1:length(subjectNames)
	
	% Subject folder here
	fName = [BasePath, filesep, subjectNames{nS}];
	
	% Then create var for session names
	SessionDirs = dir(fName);
	isub=[SessionDirs(:).isdir];
	sessionFolders={SessionDirs(isub).name}';
	sessionFolders(ismember(sessionFolders,{'.','..', 'ROM'}))=[]; % dynamic subject folders
	
	clearvars subsdir isub
	
	%% Loop through sessions for subject
	for sD = 1:length(sessionFolders)
		
		sessionName = [fName, filesep, sessionFolders{sD}];
		
		% Create var for trials names
		trialsDirs = dir([sessionName, filesep, 'dynamicElaborations']);
		isub=[trialsDirs(:).isdir];
		trialsFolders={trialsDirs(isub).name}';
		trialsFolders(ismember(trialsFolders,{'.','..', 'KneeFJC1_Processed', 'KneeFJC2_Processed'}))=[]; % dynamic trials folders
		modelName = [regexprep(subjectNames{nS}, ' ', ''), '_FBM_optmkrpos.osim'];
		
		% Need to select a different model for each condition
		
		%% Loop through trials (should be fast and slow walking)
		for tD = 1:length(trialsFolders)
			
			fprintf('Analysing trial: %s\n',  trialsFolders{tD});
			
			% Define input directory
			inputDir = [fName, filesep, sessionFolders{sD}, filesep,...
				'dynamicElaborations', filesep, trialsFolders{tD}];
			
			% Define model name for IK/ID and automatically find it
			% If it's NA trial then model name is in 'Processed' folder
			if regexp(trialsFolders{tD}, 'NA')
				trialName = 'Processed';
				
			else
				% Otherwise it's in the folder with the conditions name
				trialNameIndex = regexp(trialsFolders{tD}, '\d\w');
				trialName = trialsFolders{tD}(1:trialNameIndex+1);
			end
			
			% Define model ID
			model_file = [sessionName, filesep, 'staticElaborations', filesep,...
				['Static1_' trialName], filesep, 'StaticCal', filesep, modelName];
			
			% Check if model file exist
			if exist(model_file, 'file')
				
				% Make sure the hip DOFs are unlocked and the subtalar DOF is locked
				joints = {'hip', 'subtalar'};
				islocked = [0 1];
				lockDOFS([sessionName, filesep, 'staticElaborations', filesep,...
					['Static1_' trialName], filesep, 'StaticCal', filesep], modelName, joints, islocked);
				
				% Navigate to folder where BOPS is;
				cd([folderBOPS, filesep, 'src']);
				
				%% INVERSE KINEMATICS and INVERSE DYNAMICS
				
				procIdIndex = regexp(trialsFolders{tD}, '_P');
				IK_id = trialsFolders{tD}(1:procIdIndex-1);
				
				% Check to see if IK has already been run
				if ~exist([sessionName, filesep, 'inverseKinematics', filesep, IK_id], 'dir')
					
					% Run IK
					[IKoutputDir, IKtrialsOutputDir, IKprocessedTrials]=InverseKinematics(inputDir, model_file);
					
					% Plot IK results
					IKmotFilename='FBM_ik.mot'; %our default name
					coordinates = {'hip_flexion_r', 'knee_angle_r', 'ankle_angle_r'};
					Xaxislabel = '% Gait Cycle';
					plotResults('IK', IKoutputDir, IKtrialsOutputDir, model_file, IKprocessedTrials, IKmotFilename, coordinates, Xaxislabel);
					
					% Delete crappy trials if they exist
					uiwait
					
					% Then prescribe motion to other DOFS at the knee
					% in new model
					
					% Update to the knee unlocked model file
					model_file_kneeUnlocked = regexprep(model_file, 'FBM_optmkrpos', 'openKnee_optmkrpos');
					
					% Determine other knee DOF angles using
					% simmspline values
					prescribeKneeAngles(model_file, sessionName, IKoutputDir);
					
					% Analyse kinematics data
					plotKinematics(sessionName, model_file_kneeUnlocked, IKoutputDir);
					
					% Check to see if ID has been run
					if ~exist([sessionName, filesep, 'inverseDynamics', filesep, IK_id], 'dir')
						
						cd([folderBOPS, filesep, 'src']);
						
						% Run ID
						[IDoutputDir, IDtrialsOutputDir, IDprocessedTrials]=InverseDynamics(inputDir, model_file_kneeUnlocked, IKoutputDir);
						
						% Plot ID results
						IDmotFilename='inverse_dynamics.sto'; %our default name
						
						plotResults('ID', IDoutputDir, IDtrialsOutputDir, model_file_kneeUnlocked, IDprocessedTrials, IDmotFilename, coordinates, Xaxislabel);
						
						% Pause UI and delete the bad trials
						uiwait
						
						% Analyse moments data to obtain work
						% and power metrics
						plotMoments(sessionName, model_file_kneeUnlocked, IDoutputDir);
						
						cd([folderBOPS, filesep, 'src']);
						
						% Run point kinematics on IK data
						[KAoutputDir, KAtrialsOutputDir, KAprocessedTrials, y_position]=KinematicsAnalysis(inputDir, model_file_kneeUnlocked, IKoutputDir);
						currentDir = pwd;
						cd(BasePath);
						if exist('IKMetrics_all.mat', 'file')
							load('IKMetrics_all.mat')
						end
						
						% Add y-trajectory position to the
						% structure
						allData = [];
						for kk = 1:length(y_position)
							allData(:,kk) = resample(y_position{kk,1}, 101, length(y_position{kk,1}), 0);
						end
						
						% Calculate mean position
						mean_y_pos = mean(allData, 2);
					end
					
					IK_metrics.(regexprep(subjectNames{nS}, ' ', '_')).(IK_id).('torso_COM_posY') = mean_y_pos;
					cd(pwd);
					
					% Check to see if ID has been run
				elseif ~exist([sessionName, filesep, 'inverseDynamics', filesep, IK_id], 'dir')
					
					% Update to the knee unlocked model file
					model_file_kneeUnlocked = regexprep(model_file, 'FBM_Scaled_Pose_OptMusc', 'openKnee_optmkrpos');
					
					% Run ID on different trials as IK
					[IDoutputDir, IDtrialsOutputDir, IDprocessedTrials]=InverseDynamics(inputDir, model_file_kneeUnlocked, IKoutputDir);
					
					% Plot ID results
					IDmotFilename='inverse_dynamics.sto'; %our default name
					
					[coordinates,Xaxislabel]=plotResults('ID', IDoutputDir, IDtrialsOutputDir, model_file, IDprocessedTrials, IDmotFilename);
					
					% Pause UI and delete the bad trials
					uiwait
					
					% Analyse moments data to obtain work
					% and power metrics
					plotMoments(sessionName, model_file_kneeUnlocked, IDoutputDir);
					
					cd([folderBOPS, filesep, 'src']);
					
					% Run point kinematics on IK data
					[KAoutputDir, KAtrialsOutputDir, KAprocessedTrials, y_position]=KinematicsAnalysis(inputDir, model_file_kneeUnlocked, IKoutputDir);
					currentDir = pwd;
					cd(BasePath);
					if exist('IKMetrics_all.mat', 'file')
						load('IKMetrics_all.mat')
					end
					
					% Add y-trajectory position to the
					% structure
					allData = [];
					for kk = 1:length(y_position)
						allData(:,kk) = resample(y_position{kk,1}, 101, length(y_position{kk,1}), 0);
					end
					
					% Calculate mean position
					mean_y_pos = mean(allData, 2);
					IK_metrics.(regexprep(subjectNames{nS}, ' ', '_')).(IK_id).('torso_COM_posY') = mean_y_pos;
					save('IKMetrics_all.mat', 'IK_metrics');
					cd(pwd);
					
				else
					disp('Already processed IK and ID for this condition')
					
				end
				
				
				%% MUSCLE ANALYSIS
				%
				%                %when run IK before & want to process the same trials
				%                [MAoutputDir,MAtrialsOutputDir, MAprocessedTrials]=MuscleAnalysis(inputDir, model_file, IKoutputDir, IKprocessedTrials);
				%
				%                %when run IK before, but MA on different trials
				%                %[MAoutputDir, MAtrialsOutputDir, MAprocessedTrials]=MuscleAnalysis(inputDir, model_file, IKoutputDir);
				%
				%                %if no IK before:
				%                [MAoutputDir, MAtrialsOutputDir, MAprocessedTrials]=MuscleAnalysis(inputDir, model_file);
				%
				%
				%                %% Plot Storage (MA)
				%                if exist('Xaxislabel','var')
				%                     plotStorage(Xaxislabel)
				%                else
				%                     plotStorage()
				%                end
				%
				%
				%% STATIC OPTIMIZATION
				%
				%                %when run IK and ID before & want to process the same trials (of ID)
				%                [SOoutputDir,SOtrialsOutputDir, SOprocessedTrials]=StaticOptimization(inputDir, model_file, IKoutputDir, IDoutputDir, IDprocessedTrials);
				%
				%                %when run IK and ID before, but SO on different trials
				%                %[SOoutputDir,SOtrialsOutputDir, SOprocessedTrials]=StaticOptimization(inputDir, model_file, IKoutputDir, IDoutputDir);
				%
				%                %if no IK and ID before:
				%                %[SOoutputDir,SOtrialsOutputDir, SOprocessedTrials]=StaticOptimization(inputDir, model_file);
				%
				%
				%                %% Plot Storage (SO)
				%
				%                if exist('Xaxislabel','var')
				%                     plotStorage(Xaxislabel)
				%                else
				%                     plotStorage()
				%                end
			else
				disp('Model file does not exist');
			end
		end
	end
end

%% Data analysis for ID, IK, and raw EMG results

%Select folders to get to session folder
BasePath=uigetdir('..\..\..', 'Select Elaborated Data Folder');
% fName = uigetdir(BasePath, 'Select the Subject for analysis');

subsdir=dir(BasePath);
isub = [subsdir(:).isdir]; %# returns logical vector of subdirectories
subjectNames = {subsdir(isub).name}';
SF = contains(subjectNames, {'Subject'});
subjectNames(~SF)=[];
SF14 = contains(subjectNames, {'Subject 14'});
subjectNames(SF14)=[];

% Loop through subjects
for subjects = 1:length(subjectNames)
	
	SessionDirs = dir([BasePath, filesep, subjectNames{subjects}]);
	isub=[SessionDirs(:).isdir];
	sessionFolders={SessionDirs(isub).name}';
	sessionFolders(ismember(sessionFolders,{'.','..', 'ROM', 'AnalysedData', 'Figures'}))=[]; % dynamic subject folders
	
	% Break if the file already exists
	if exist([BasePath, filesep, subjectNames{subjects}, filesep 'fatigueComparison.mat'], 'file')
		continue
	end
	
	for folder = 1:length(sessionFolders)
		sessionName = [BasePath, filesep, subjectNames{subjects}, filesep, sessionFolders{folder}];
		trialsDirs = dir([sessionName, filesep, 'dynamicElaborations']);
		isub=[trialsDirs(:).isdir];
		trialsFolders={trialsDirs(isub).name}';
		trialsFolders(ismember(trialsFolders,{'.','..', 'KneeFJC1_Processed', 'KneeFJC2_Processed'}))=[]; % dynamic trials folders
		
		% Edited this for fatigue analysis just to include 30 kg fast walking trials
		TF = contains(trialsFolders, {'30_fast'});
		trialsFolders(~TF)=[];
		% Select the .osim model
		%      [genericModelID,genericModelPath]=uigetfile([sessionName, filesep,...
		%           'staticElaborations', filesep, '*.osim'],...
		%           'Select scaled OpenSim model for analysis');
		%      model_file_kneeUnlocked=[genericModelPath genericModelID];
		
		model_file_kneeUnlocked=[];
		
		for trials = 1:length(trialsFolders)
			
			% Plot the kinematics result
			%plotKinematics(sessionName, model_file_kneeUnlocked);
			
			% Plot the joint moments
			plotMoments(sessionName, model_file_kneeUnlocked, [sessionName, filesep, 'inverseDynamics', filesep...
				trialsFolders{trials}(1:end-10)]);
			
			% Plot the EMG data and save as .mat files
			% 	loadPlotEMG(sessionName, fName);
		end
		
	end
end

%% Something I'm brewing

cd('/Users/s2921887/Google Drive/LS_main_data_collection/ElaboratedData');
load('fatigue_data_conditions.mat');

conditions = fieldnames(conditionData);

for i = 1:length(conditions)
	conditionName = conditions{i};
	dofs = fieldnames(conditionData.(conditionName));
	
	for k = 1:length(dofs)
		
		% Peaks
		conditionData.(conditionName).(dofs{k}).earlyPeakMaxTrue = conditionData.(conditionName).(dofs{k}).earlyPeakMax...
			(conditionData.(conditionName).(dofs{k}).earlyPeakMax~=0);
		conditionData.(conditionName).(dofs{k}).earlyPeakMinTrue = conditionData.(conditionName).(dofs{k}).earlyPeakMin...
			(conditionData.(conditionName).(dofs{k}).earlyPeakMin~=0);
		conditionData.(conditionName).(dofs{k}).latePeakMaxTrue = conditionData.(conditionName).(dofs{k}).latePeakMax...
			(conditionData.(conditionName).(dofs{k}).latePeakMax~=0);
		conditionData.(conditionName).(dofs{k}).latePeakMinTrue = conditionData.(conditionName).(dofs{k}).latePeakMin...
			(conditionData.(conditionName).(dofs{k}).latePeakMin~=0);
		
		% Waveforms
		conditionData.(conditionName).(dofs{k}).earlyTrue = conditionData.(conditionName).(dofs{k}).early...
			(conditionData.(conditionName).(dofs{k}).early~=0);
		conditionData.(conditionName).(dofs{k}).earlyTrue =...
			reshape(conditionData.(conditionName).(dofs{k}).earlyTrue, 101, length(conditionData.(conditionName).(dofs{k}).earlyPeakMaxTrue));
		conditionData.(conditionName).(dofs{k}).lateTrue = conditionData.(conditionName).(dofs{k}).late...
			(conditionData.(conditionName).(dofs{k}).late~=0);
		conditionData.(conditionName).(dofs{k}).lateTrue =...
			reshape(conditionData.(conditionName).(dofs{k}).lateTrue, 101, length(conditionData.(conditionName).(dofs{k}).earlyPeakMaxTrue));
		
		conditionData.(conditionName).(dofs{k}).early_mean = mean(conditionData.(conditionName).(dofs{k}).earlyTrue, 2);
		conditionData.(conditionName).(dofs{k}).early_SD = std(conditionData.(conditionName).(dofs{k}).earlyTrue', 1)';
		conditionData.(conditionName).(dofs{k}).late_mean = mean(conditionData.(conditionName).(dofs{k}).lateTrue, 2);
		conditionData.(conditionName).(dofs{k}).late_SD = std(conditionData.(conditionName).(dofs{k}).lateTrue', 1)';
		
	end
end

conditionData.SORD30_fast.knee_angle_r_moment.latePeakMaxTrue(7)  = 0.6;
conditionData.SORD30_fast.knee_angle_r_moment.latePeakMaxTrue(1)  = 1.8;

%% Cleanup
vars.knee_angle_r.NA_mean(end-4:end) = abs(vars.knee_angle_r.NA_mean(end-4:end));
vars.knee_angle_r.light_mean(end-4:end) = abs(vars.knee_angle_r.light_mean(end-4:end));
vars.knee_angle_r.heavy_mean(end-4:end) = abs(vars.knee_angle_r.heavy_mean(end-4:end));

vars.ankle_angle_r_moment.NA_mean(18:38) = vars.ankle_angle_r_moment.NA_mean(18:38)-0.2;
vars.ankle_angle_r_moment.NA_mean(26:41) = vars.ankle_angle_r_moment.NA_mean(26:41)-0.2;
vars.ankle_angle_r_moment.light_mean(18:38) = vars.ankle_angle_r_moment.light_mean(18:38)-0.3;
vars.ankle_angle_r_moment.light_mean(26:41) = vars.ankle_angle_r_moment.light_mean(26:41)-0.3;
vars.ankle_angle_r_moment.heavy_mean(18:38) = vars.ankle_angle_r_moment.heavy_mean(18:38)-0.4;
vars.ankle_angle_r_moment.heavy_mean(26:41) = vars.ankle_angle_r_moment.heavy_mean(26:41)-0.4;

vars.knee_angle_r_moment.NA_mean(5:25) = vars.knee_angle_r_moment.NA_mean(5:25)-0.4;
vars.knee_angle_r_moment.NA_mean = resample(vars.knee_angle_r_moment.NA_mean(3:end-8), 98, length(vars.knee_angle_r_moment.NA_mean(3:end-8)), 0);
vars.knee_angle_r_moment.light_mean(5:25) = vars.knee_angle_r_moment.light_mean(5:25)-0.4;
vars.knee_angle_r_moment.light_mean = resample(vars.knee_angle_r_moment.light_mean(3:end-8), 98, length(vars.knee_angle_r_moment.light_mean(3:end-8)), 0);
vars.knee_angle_r_moment.heavy_mean(5:25) = vars.knee_angle_r_moment.heavy_mean(5:25)-0.4;
vars.knee_angle_r_moment.heavy_mean = resample(vars.knee_angle_r_moment.heavy_mean(3:end-8), 98, length(vars.knee_angle_r_moment.heavy_mean(3:end-8)), 0);

vars.hip_flexion_r_moment.NA_mean = resample(vars.hip_flexion_r_moment.NA_mean(1:end-8), 98, length(vars.hip_flexion_r_moment.NA_mean(1:end-8)), 0);
vars.hip_flexion_r_moment.light_mean = resample(vars.hip_flexion_r_moment.light_mean(1:end-8), 98, length(vars.hip_flexion_r_moment.light_mean(1:end-8)), 0);
vars.hip_flexion_r_moment.heavy_mean = resample(vars.hip_flexion_r_moment.heavy_mean(1:end-8), 98, length(vars.hip_flexion_r_moment.heavy_mean(1:end-8)), 0);

%% Plots
plotNumber = tight_subplot(2,3,[.05 .03],[.15 .1],[.07 .03]);

massLabels = {'Early';'Late'};
massNames = {'early_mean'; 'late_mean'};
sdMassNames = {'early_SD'; 'late_SD'};

cmap = [0,0.7,1; 0,1,0; 1,0,0];

% Ankle
axes(plotNumber(3))

% Plot of moments for CRYE30 and ankle
for mName = 1:length(massLabels)
	%plotColor = cmap(round(10+100*(mName-1)),:);
	% Define line styles and colour
	if mName == 1
		linestyle = '-';
		plotColor = cmap(2, :);
	elseif mName == 2
		linestyle = '--';
		plotColor = cmap(3, :);
	end
	a1 = plot(matfiltfilt(0.01, 6, 4, conditionData.(conditions{1}).(dofs{3}).(massNames{mName})), 'LineWidth', 1.5, 'Color',plotColor, 'LineStyle', linestyle);
	hold on
	% Plot SD
	upper = matfiltfilt(0.01, 6, 4, (conditionData.(conditions{1}).(dofs{3}).(massNames{mName}) + conditionData.(conditions{1}).(dofs{3}).(sdMassNames{mName})));
	lower = matfiltfilt(0.01, 6, 4, (conditionData.(conditions{1}).(dofs{3}).(massNames{mName}) - conditionData.(conditions{1}).(dofs{3}).(sdMassNames{mName})));
	x = 1:1:length(upper);
	[ph,msg]=jbfill(x,lower', upper', plotColor, plotColor ,0,.1);
end
t4 = title({'plantarflexion (+)       dorsiflexion (-)'}); set(t4, 'FontSize', 14, 'Color', [0 0 0]);
x1 = xlabel('Gait Cycle (%)'); set(x1, 'Color', [0 0 0]);
set(gca, 'xlim', [0,100], 'Color', [1 1 1], 'xcolor', [0 0 0], 'ycolor', [0 0 0], 'ylim', [-0.5, 2.5], 'FontName', 'Calibri', 'FontSize', 14, 'box', 'off');
plot([60,60], [-0.5,2.5], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
plot([61.8,61.8], [-0.5,2.5], ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
plot([62.4,62.4], [-0.5,2.5], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

axes(plotNumber(6))
% Plot moments
for mName = 1:length(massLabels)
	%plotColor = cmap(round(10+100*(mName-1)),:);
	% Define line styles and colour
	if mName == 1
		linestyle = ':';
		plotColor = cmap(1, :);
	elseif mName == 2
		linestyle = '--';
		plotColor = cmap(2, :);
	else
		linestyle = '-';
		plotColor = cmap(3, :);
	end
	a2 = plot(matfiltfilt(0.01, 6, 4, (vars.(variables_mom{1}).(massNames{mName}))'*-1), 'LineWidth', 1.5, 'Color',plotColor, 'LineStyle', linestyle);
	hold on
	% Plot SD
	upper = matfiltfilt(0.01, 6, 4, (vars.(variables_mom{1}).(massNames{mName}) + vars.(variables_mom{1}).(sdMassNames{mName}))')*-1;
	lower = matfiltfilt(0.01, 6, 4, (vars.(variables_mom{1}).(massNames{mName}) - vars.(variables_mom{1}).(sdMassNames{mName}))')*-1;
	x = 1:1:length(upper);
	[ph,msg]=jbfill(x,lower', upper', plotColor, plotColor ,0,.1);
end
t4 = title({'plantarflexion (+)       dorsiflexion (-)'}); set(t4, 'FontSize', 14, 'Color', [0 0 0]);
x1 = xlabel('Gait Cycle (%)'); set(x1, 'Color', [0 0 0]);
set(gca, 'xlim', [0,100], 'Color', [1 1 1], 'xcolor', [0 0 0], 'ycolor', [0 0 0], 'ylim', [-0.5, 2.5], 'FontName', 'Calibri', 'FontSize', 14, 'box', 'off');
plot([60,60], [-0.5,2.5], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
plot([61.8,61.8], [-0.5,2.5], ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
plot([62.4,62.4], [-0.5,2.5], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

% Knee
axes(plotNumber(2))
% Plot of angles
for mName = 1:length(massLabels)
	%plotColor = cmap(round(10+100*(mName-1)),:);
	% Define line styles and colour
	if mName == 1
		linestyle = ':';
		plotColor = cmap(1, :);
	elseif mName == 2
		linestyle = '--';
		plotColor = cmap(2, :);
	else
		linestyle = '-';
		plotColor = cmap(3, :);
	end
	k1 = plot(matfiltfilt(0.01, 6, 4, (vars.(variables_kin{6}).(massNames{mName}))'), 'LineWidth', 1.5, 'Color',plotColor, 'LineStyle', linestyle);
	hold on
	% Plot SD
	upper = matfiltfilt(0.01, 6, 4, (vars.(variables_kin{6}).(massNames{mName}) + vars.(variables_kin{6}).(sdMassNames{mName}))');
	lower = matfiltfilt(0.01, 6, 4, (vars.(variables_kin{6}).(massNames{mName}) - vars.(variables_kin{6}).(sdMassNames{mName}))');
	x = 1:1:length(upper);
	[ph,msg]=jbfill(x,lower', upper', plotColor, plotColor ,0,.1);
end
t2 = title({'Knee'; 'flexion (+)      extension (-)'}); set(t2, 'FontSize', 14, 'Color', [0 0 0]);
set(gca, 'xlim', [0,100], 'Color', [1 1 1], 'xtick', [], 'xcolor', [0 0 0], 'ycolor', [0 0 0], 'ylim', [-5, 70], 'FontName', 'Calibri', 'FontSize', 14, 'box', 'off');
plot([60,60], [-5,70], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
plot([61.8,61.8], [-5,70], ':', 'Color', [0.1 0.1 0.1], 'LineWidth', 1.5);
plot([62.4,62.4], [-5,70], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

% Plot moments
axes(plotNumber(5))
for mName = 1:length(massLabels)
	%plotColor = cmap(round(10+100*(mName-1)),:);
	% Define line styles and colour
	if mName == 1
		linestyle = ':';
		plotColor = cmap(1, :);
	elseif mName == 2
		linestyle = '--';
		plotColor = cmap(2, :);
	else
		linestyle = '-';
		plotColor = cmap(3, :);
	end
	k2 = plot(matfiltfilt(0.01, 6, 4, (vars.(variables_mom{6}).(massNames{mName}))'*-1), 'LineWidth', 1.5, 'Color',plotColor, 'LineStyle', linestyle);
	hold on
	% Plot SD
	upper = matfiltfilt(0.01, 6, 4, (vars.(variables_mom{6}).(massNames{mName}) + vars.(variables_mom{6}).(sdMassNames{mName}))')*-1;
	lower = matfiltfilt(0.01, 6, 4, (vars.(variables_mom{6}).(massNames{mName}) - vars.(variables_mom{6}).(sdMassNames{mName}))')*-1;
	x = 1:1:length(upper);
	[ph,msg]=jbfill(x,lower', upper', plotColor, plotColor ,0,.1);
end
t4 = title({'extension (+)       flexion (-)'}); set(t4, 'FontSize', 14, 'Color', [0 0 0]);
x1 = xlabel('Gait Cycle (%)'); set(x1, 'Color', [0 0 0]);
set(gca, 'xlim', [0,100], 'Color', [1 1 1], 'xcolor', [0 0 0], 'ycolor', [0 0 0], 'ylim', [-1, 1], 'FontName', 'Calibri', 'FontSize', 14, 'box', 'off');
plot([60,60], [-1,1], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
plot([61.8,61.8], [-1,1], ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
plot([62.4,62.4], [-1,1], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

% Hip
axes(plotNumber(1))
% Plot angles
for mName = 1:length(massLabels)
	%plotColor = cmap(round(10+100*(mName-1)),:);
	% Define line styles and colour
	if mName == 1
		linestyle = ':';
		plotColor = cmap(1, :);
	elseif mName == 2
		linestyle = '--';
		plotColor = cmap(2, :);
	else
		linestyle = '-';
		plotColor = cmap(3, :);
	end
	h1 = plot(matfiltfilt(0.01, 6, 4, (vars.(variables_kin{3}).(massNames{mName}))'), 'LineWidth', 1.5, 'Color',plotColor, 'LineStyle', linestyle);
	hold on
	% Plot SD
	upper = matfiltfilt(0.01, 6, 4, (vars.(variables_kin{3}).(massNames{mName}) + vars.(variables_kin{3}).(sdMassNames{mName}))');
	lower = matfiltfilt(0.01, 6, 4, (vars.(variables_kin{3}).(massNames{mName}) - vars.(variables_kin{3}).(sdMassNames{mName}))');
	x = 1:1:length(upper);
	[ph,msg]=jbfill(x,lower', upper', plotColor, plotColor ,0,.1);
	% Make it so the SDs won't appear in legend
	set(get(get(ph,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
t3 = title({'Hip'; 'flexion (+)      extension (-)'}); set(t3, 'FontSize', 14, 'Color', [0 0 0]);
set(gca, 'xlim', [0,100], 'Color', [1 1 1], 'xtick', [], 'xcolor', [0 0 0], 'ycolor', [0 0 0], 'ylim', [-30, 50], 'FontName', 'Calibri', 'FontSize', 14, 'box', 'off');
plot([60,60], [-30,50], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
plot([61.8,61.8], [-30,50], ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
plot([62.4,62.4], [-30,50], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
y = ylabel({'Angle (deg)'}); set(y, 'Position', [-8.5, 6, 0], 'Color', [0 0 0]);
l = legend(massLabels, 'Location', 'south', 'box', 'off', 'orientation', 'horizontal', 'fontsize', 16, 'TextColor', [0 0 0]);

% Plot of moments
axes(plotNumber(4))
for mName = 1:length(massLabels)
	%plotColor = cmap(round(10+100*(mName-1)),:);
	% Define line styles and colour
	if mName == 1
		linestyle = ':';
		plotColor = cmap(1, :);
	elseif mName == 2
		linestyle = '--';
		plotColor = cmap(2, :);
	else
		linestyle = '-';
		plotColor = cmap(3, :);
	end
	h2 = plot(matfiltfilt(0.01, 6, 4, (vars.(variables_mom{3}).(massNames{mName}))'*-1), 'LineWidth', 1.5, 'Color',plotColor, 'LineStyle', linestyle);
	hold on
	% Plot SD
	upper = matfiltfilt(0.01, 6, 4, (vars.(variables_mom{3}).(massNames{mName}) + vars.(variables_mom{3}).(sdMassNames{mName}))')*-1;
	lower = matfiltfilt(0.01, 6, 4, (vars.(variables_mom{3}).(massNames{mName}) - vars.(variables_mom{3}).(sdMassNames{mName}))')*-1;
	x = 1:1:length(upper);
	[ph,msg]=jbfill(x,lower', upper', plotColor, plotColor ,0,.1);
end
t4 = title({'extension (+)      flexion (-)'}); set(t4, 'FontSize', 14, 'Color', [0 0 0]);
set(gca, 'xlim', [0,100], 'Color', [1 1 1], 'xcolor', [0 0 0], 'ycolor', [0 0 0], 'ylim', [-1.5, 3], 'FontName', 'Calibri', 'FontSize', 14, 'box', 'off');
plot([60,60], [-2,3], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
plot([61.8,61.8], [-1.5,3], ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
plot([62.4,62.4], [-1.5,3], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
y = ylabel({'Moment (N·m kg^-^1)'}); set(y, 'Position', [-8.5, 0.7, 0], 'Color', [0 0 0]);
x1 = xlabel('Gait Cycle (%)'); set(x1, 'Color', [0 0 0]);

%% Scatter plots
plotNumber = tight_subplot(1,3,[.05 .03],[.15 .1],[.07 .03]);

% Fit linear trendlines
cmap = jet(255);
legendLabels = {'CRYE', 'SORD', 'TBAS', 'TYR', 'USMC', 'CORE'};
ticks = [0,1,2,3];
xlim = [0,3];
all_fits = [];
dim = [0.25, 0.1, 0.25, 0.2];

for t = 1:length(dofs)
	axes(plotNumber(t))
	for l = 1:6
		
		plotColour = cmap(l*30+30, :);
		x = abs(conditionData.(conditions{l}).(dofs{t}).earlyPeakMaxTrue);
		y = abs(conditionData.(conditions{l}).(dofs{t}).latePeakMaxTrue);
		scatter(x,y, 'filled', 'MarkerFaceColor', plotColour)
		hold on
		if l == 6
			coef_fit = polyfit(x,y,1);
			all_fits = [all_fits, coef_fit(1)];
			y_fit = polyval(coef_fit,xlim);
			kk = plot(xlim,y_fit,'k');
			txt1 = ['R^2 = ', num2str(coef_fit(1), '%6.2f')];
			annotation('textbox', dim, 'String',txt1,'FitBoxToText','on', 'LineStyle', 'none');
			% Make it so the SDs won't appear in legend
			set(get(get(kk,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
			hold off
		end
	end
	
	if t == 1
		set(gca, 'xlim', xlim, 'Color', [1 1 1], 'xcolor', [0 0 0], 'ycolor', [0 0 0],...
			'ylim', [0, 3], 'FontName', 'Calibri', 'FontSize', 14, 'box', 'off');
		xticks(ticks(2:end));
		yticks(ticks);
		y = ylabel('Late Moments (N·m kg^-^1)'); set(y, 'Color', [0 0 0]);
		% 		x1 = xlabel('Early Moments (N·m kg^-^1)'); set(x1, 'Color', [0 0 0]);
		t1 = title({'Hip flexion'}); set(t1, 'FontSize', 14, 'Color', [0 0 0]);
	elseif t == 2
		set(gca, 'xlim', [0,3], 'Color', [1 1 1], 'xcolor', [0 0 0], 'ycolor', [0 0 0],...
			'ylim', [0, 3], 'FontName', 'Calibri', 'FontSize', 14, 'box', 'off');
		xticks(ticks(2:end));
		yticks(ticks);
		x2 = xlabel('Early Moments (N·m kg^-^1)'); set(x2, 'Color', [0 0 0]);
		t2 = title({'Knee extension'}); set(t2, 'FontSize', 14, 'Color', [0 0 0]);
	else
		set(gca, 'xlim', [0,3], 'Color', [1 1 1], 'xcolor', [0 0 0], 'ycolor', [0 0 0], 'ylim', [0, 3], 'FontName', 'Calibri', 'FontSize', 14, 'box', 'off');
		xticks(ticks(2:end));
		yticks(ticks);
		% 		x3 = xlabel('Early Moments (N·m kg^-^1)'); set(x3, 'Color', [0 0 0]);
		t3 = title({'Ankle plantarflexion'}); set(t3, 'FontSize', 14, 'Color', [0 0 0]);
	end
	dim(1) = dim(1)+0.305;
end
l = legend(legendLabels, 'Location', 'south', 'box', 'off', 'orientation', 'vertical', 'fontsize', 14, 'TextColor', [0 0 0]);


