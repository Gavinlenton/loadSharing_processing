function plotMoments(sessionName, model_file, IDoutputDir)
%Plot multiple moments from the .sto files generated from OpenSim's ID
%analysis
%   Input model file directory and session name directory to generate individual
%   joint moment figures of the DOFs of interest (e.g., hip flexion).
%   Multiple trials can be plotted on the same figure for each DOF

if nargin < 3
     % Folder where the ID Results files are stored
     IDoutputDir = uigetdir([sessionName, filesep, trialName], 'Select folder with INVERSE DYNAMICS results to use');
end

% Generate list of trials
trials=dir(IDoutputDir);
j=1;
for k = 3:length(trials)
     trialsList{j}=trials(k).name;
     j = j + 1;
end
trialsList(ismember(trialsList,{'Figures','IDMetrics.mat', '.DS_Store'}))=[];

 % Be selective if you want to
[trialsIndex,~] = listdlg('PromptString','Select ID trials to plot:',...
     'SelectionMode','multiple',...
     'ListString',trialsList);

inputTrials=trialsList(trialsIndex);

% inputTrials=trialsList;

% Define subject weight for normalisation
% Folder containing acquisition xml
acquisitionFolder = [regexprep(sessionName, 'ElaboratedData', 'InputData'), filesep];
acquisitionInfo=xml_read(fullfile(acquisitionFolder, 'acquisition.xml'));
subject_weight = acquisitionInfo.Subject.Weight;
subject_name = acquisitionInfo.Subject.Code;
condition_nameIndex = regexp(inputTrials{1}, '_P');
condition_name = inputTrials{1}(1:condition_nameIndex-1);

% Plot multiple trials
IDfilename='inverse_dynamics.sto';

% Automatic selection of DOFs 
dofsToPlot = {'hip_flexion_r'; 'hip_adduction_r'; 'hip_rotation_r'; 'knee_angle_r';...
     'knee_adduction_r'; 'knee_rotation_r'; 'ankle_angle_r'; 'lumbar_extension';...
     'lumbar_bending'; 'lumbar_rotation'};
 

% UNCOMMENT THIS TO MANUALLY SELECT WHICH DOFS TO PLOT
% [selectedDofsIndex,v] = listdlg('PromptString','Select dofs for plots:',...
%      'SelectionMode','multiple',...
%      'ListString',dofs);
% dofs=getDofsFromModel(model_file);
% Assign dofs to plot
% dofsToPlot=dofs(selectedDofsIndex)';

% Pre-allocate
momentsToPlot = cell(length(dofsToPlot),1);

% Convert dofs to names of moments output in .sto file
for i=1:length(dofsToPlot)
     momentsToPlot{i}=strcat(dofsToPlot{i}, '_moment');
end

% Assign x-axis label
xaxislabel = '% Gait Cycle';

% Directory to eventually save all the data
elabDataFolder = sessionName(1:end-10);

% Plot multiple results on a figure per DOF
% The moments plotted from OpenSim are the inverse of what is typically
% seen, you can choose to invert the results if desired.
plotResultsMultipleTrials_LS(IDoutputDir, elabDataFolder, inputTrials,  IDfilename, xaxislabel, momentsToPlot, subject_weight, subject_name, condition_name)


