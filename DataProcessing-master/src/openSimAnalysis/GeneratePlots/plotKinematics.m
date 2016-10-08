function plotKinematics(sessionName, model_file)
%Plot multiple moments from the .sto files generated from OpenSim's ID
%analysis
%   Input model file directory and session name directory to generate individual
%   joint moment figures of the DOFs of interest (e.g., hip flexion).
%   Multiple trials can be plotted on the same figure for each DOF

% Folder where the IK Results files are stored
IKresultsDir = uigetdir(sessionName, 'Select folder with INVERSE KINEMATICS results to use');

% Generate list of trials
trials=dir(IKresultsDir);
j=1;

for k = 3:length(trials)
     trialsList{j}=trials(k).name;
     j = j + 1;
end
trialsList(ismember(trialsList,{'Figures','IDMetrics.mat', 'out.log', 'error.log', 'AnalysedData'}))=[];

% Be selective if you want to
[trialsIndex,~] = listdlg('PromptString','Select trials to plot:',...
     'SelectionMode','multiple',...
     'ListString',trialsList);

inputTrials=trialsList(trialsIndex);

% Generate full path to files
for n = 1:length(inputTrials)
     trialsListDir{n} = [IKresultsDir filesep inputTrials{n}];
end

% Define subject weight for normalisation
% Folder containing acquisition xml
acquisitionFolder = [regexprep(sessionName, 'ElaboratedData', 'InputData'), filesep];
acquisitionInfo=xml_read(fullfile(acquisitionFolder, 'acquisition.xml'));
subject_weight = acquisitionInfo.Subject.Weight;
subject_name = acquisitionInfo.Subject.Code;
condition_nameIndex = regexp(inputTrials{1}, '_P');
condition_name = inputTrials{1}(1:condition_nameIndex-1);

% Plot multiple trials
IKfilename='openKnee_ik.mot';
dofs=getDofsFromModel(model_file);
[selectedDofsIndex,v] = listdlg('PromptString','Select dofs for plots:',...
     'SelectionMode','multiple',...
     'ListString',dofs);

% Assign dofs to plot
dofsToPlot=dofs(selectedDofsIndex)';

% Assign x-axis label
xaxislabel = '% Gait Cycle';

% Directory to eventually save all the data
elabDataFolder = sessionName(1:end-10);

% Plot multiple results on a figure per DOF
% The angles plotted from OpenSim are the inverse of what is typically
% seen, you can choose to invert the results if desired.
plotResultsMultipleTrials_LS(IKresultsDir, elabDataFolder, inputTrials, IKfilename, xaxislabel, dofsToPlot, subject_weight, subject_name, condition_name)


