function plotPointKinematics(sessionName, PKoutputDir)
%Plot multiple moments from the .sto files generated from OpenSim's ID
%analysis
%   Input model file directory and session name directory to generate individual
%   joint moment figures of the DOFs of interest (e.g., hip flexion).
%   Multiple trials can be plotted on the same figure for each DOF

if nargin < 2
     % Folder where the Point kinematics Results files are stored
     PKoutputDir = uigetdir([sessionName, filesep, trialName], 'Select folder with POINT KINEMATICS results to use');
end

% Generate list of trials
trials=dir(PKoutputDir);
j=1;
for k = 3:length(trials)
     trialsList{j}=trials(k).name;
     j = j + 1;
end
trialsList(ismember(trialsList,{'Figures','IDMetrics.mat', '.DS_Store'}))=[];

 % Be selective if you want to
% [trialsIndex,~] = listdlg('PromptString','Select ID trials to plot:',...
%      'SelectionMode','multiple',...
%      'ListString',trialsList);

inputTrials=trialsList;

% inputTrials=trialsList;

% Define subject weight for normalisation
% Folder containing acquisition xml
acquisitionFolder = [regexprep(sessionName, 'ElaboratedData', 'InputData'), filesep];
acquisitionInfo=xml_read(fullfile(acquisitionFolder, 'acquisition.xml'));
subject_weight = acquisitionInfo.Subject.Weight;
subject_name = acquisitionInfo.Subject.Code;
condition_nameIndex = regexp(inputTrials{1}, '_P');
condition_name = inputTrials{1}(1:condition_nameIndex-1);

% File names for point kinematics analysis
PKAnalysis={'Torso_PointKinematics_COM_acc.sto'; 'Torso_PointKinematics_COM_vel.sto';  ...
	'Torso_PointKinematics_COM_pos.sto'};

% Assign x-axis label
xaxislabel = '% Gait Cycle';

% Directory to eventually save all the data
elabDataFolder = sessionName(1:end-10);

% Plot multiple results on a figure per DOF
% The moments plotted from OpenSim are the inverse of what is typically
% seen, you can choose to invert the results if desired.
plotResultsPK_LS(PKoutputDir, elabDataFolder, inputTrials, xaxislabel, PKAnalysis, subject_weight, subject_name, condition_name)

end

