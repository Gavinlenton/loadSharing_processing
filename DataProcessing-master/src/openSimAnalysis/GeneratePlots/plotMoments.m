function plotMoments(sessionName, model_file)
%Plot multiple moments from the .sto files generated from OpenSim's ID
%analysis
%   Input model file directory and session name directory to generate individual
%   joint moment figures of the DOFs of interest (e.g., hip flexion).
%   Multiple trials can be plotted on the same figure for each DOF

% Folder where the IK Results files are stored
IDresultsDir = uigetdir('..\', 'Select folder with INVERSE DYNAMICS results to use');

% Generate list of trials
trials=dir(IDresultsDir);
j=1;

for k = 3:length(trials)
     trialsList{j}=trials(k).name;
     j = j + 1;
end

% Be selective if you want to
[trialsIndex,~] = listdlg('PromptString','Select trials to plot:',...
     'SelectionMode','multiple',...
     'ListString',trialsList);

inputTrials=trialsList(trialsIndex);

% Generate full path to files
for n = 1:length(inputTrials)
     trialsListDir{n} = [IDresultsDir filesep inputTrials{n}];
end

% Define subject weight for normalisation
% Folder containing acquisition xml
acquisitionFolder = [regexprep(sessionName, 'ElaboratedData', 'InputData'), filesep];
acquisitionInfo=xml_read(fullfile(acquisitionFolder, 'acquisition.xml'));
subject_weight = acquisitionInfo.Subject.Weight;

% Plot multiple trials
IDfilename='inverse_dynamics.sto';
dofs=getDofsFromModel(model_file);
[selectedDofsIndex,v] = listdlg('PromptString','Select dofs for plots:',...
     'SelectionMode','multiple',...
     'ListString',dofs);

% Assign dofs to plot
dofsToPlot=dofs(selectedDofsIndex)';

% Convert dofs to names of moments output in .sto file
for i=1:length(dofsToPlot)
     momentsToPlot{i}=strcat(dofsToPlot{i}, '_moment');
end

% Assign x-axis label
xaxislabel = '% Gait Cycle';

% Plot multiple results on a figure per DOF
% The moments plotted from OpenSim are the inverse of what is typically
% seen, you can choose to invert the results if desired.
plotResultsMultipleTrials_LS(IDresultsDir, inputTrials, IDfilename, xaxislabel, momentsToPlot, subject_weight)


