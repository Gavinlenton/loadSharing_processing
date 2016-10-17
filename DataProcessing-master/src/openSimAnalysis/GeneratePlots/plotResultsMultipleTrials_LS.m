function []=plotResultsMultipleTrials_LS(resultsPath, elabDataFolder, trialsList, filename,x, Yquantities, varargin)
% Function to plot results from multiple trials. Also determines joint
% angular velocity and joint powers.

% This file is part of Batch OpenSim Processing Scripts (BOPS).
% Copyright (C) 2015 Alice Mantoan, Monica Reggiani
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
% Alice Mantoan, Monica Reggiani
% <ali.mantoan@gmail.com>, <monica.reggiani@gmail.com>

% Edited by Gavin Lenton October 2016 for load sharing project

%%

close all

subject_weight = varargin{1};
subject_name = regexprep(varargin{2}, ' ', '_');
condition_name = varargin{3};

% Define path names to save data - load if it alrady exists.
if strcmp(filename,'inverse_dynamics.sto')
     % Paths for ID results
     figurePath=[elabDataFolder 'Figures' filesep 'Torques' filesep];
     metricsPath = [elabDataFolder 'AnalysedData' filesep 'Torques' filesep];
     if exist(metricsPath,'dir') ~= 7
          mkdir(metricsPath);
     end
     % Load ID file if it exists
     cd(metricsPath);
     if exist('IDMetrics.mat', 'var')
          load('IDMetrics.mat');
     end
     AnalysisName = [condition_name, '_ID_raw_data'];
else
     % IK data
     figurePath=[elabDataFolder 'Figures' filesep 'Angles' filesep];
     metricsPath=[elabDataFolder 'AnalysedData' filesep 'Angles' filesep];
     % Load IK file if it exists
     if exist(metricsPath,'dir') ~= 7
          mkdir(metricsPath);
     end
     cd(metricsPath);
     if exist('IKMetrics.mat', 'var')
          load('IKMetrics.mat');
     end
     AnalysisName = [condition_name, '_IK_raw_data'];
end

% Make directory to save figures
if exist(figurePath,'dir') ~= 7
	mkdir(figurePath);
end

%Load data
for k=1:length(trialsList)
	
	file=importdata([resultsPath filesep trialsList{k} filesep filename]);
	
	if nargin>4
		
		coord_idx=findIndexes(file.colheaders,Yquantities);
	else
		Yquantities=file.colheaders(2:end); %take all columns except time
		coord_idx=[2:size(file.colheaders,2)];
	end
	
	%% CHECK THE FILENAMES FOR BOTH IK AND ID
	
     results = file.data;
     
	% Normalising moments to body weight
	if strcmp(filename,'inverse_dynamics.sto')
		results(:, 2:end) = results(:, 2:end)./subject_weight;
	end
	
	for j =1: length(coord_idx)
		
		coordCol=coord_idx(j);
		
		y{k,j} = results(:,coordCol);
		
	end
	
	timeVector{k}=getXaxis(x, results);
	
end

%Save data in mat format
save([metricsPath, AnalysisName], 'y');

% save([figurePath, 'plottedData'], 'y')
% 
% % Settings for plot
% plotLabels=regexprep(Yquantities, '_', ' ');
% legendLabels=regexprep(trialsList, '_', ' ');
% cmap = colormap(parula(180));
% plotTitle = filename;
% 
% for k=1:size(y,1)
% 	
% 	plotColor = cmap(round(1+5.5*(k-1)),:);
% 	
% 	for j=1:size(y,2)
% 		
% 		h(j)=figure(j);
% 		
% 		plot(timeVector{k}, y{k,j},'Color',plotColor)
% 		hold on
% 		
% 		xlabel(x)
% 		ylabel([plotLabels(j)])
% 		warning off
% 		legend(legendLabels)
% 		title(filename)
% 		
% 		saveas(h(j),[figurePath  Yquantities{j} '.fig'])
% 	end
% end

% After deleting bad data calculate the means and save to a structure
headers = Yquantities;

% Time in frames for angular velocity calc and filtering (change this is you didn't
% collect data at 100Hz)
dt = 1/100;

% Put the metrics in a structure
for dof = 1:length(headers)
	
    allData.(headers{dof}) = [];
	% Loop through good trials
     for trial = 1:length(trialsList)
          
          % FILTER PROCESSED DATA
          y{trial,dof} = lpfilt(y{trial,dof}, 8, dt, 'butter');
          
          % Combine data into array - resample if necessary
          % If it doesn't equal 101 then resample
          if length(y{trial, dof}) ~= 101
               allData.(headers{dof})(:, trial) = resample(y{trial, dof}, 101, length(y{trial, dof}));
          else
               allData.(headers{dof})(:, trial) = y{trial, dof};
          end
     end
     
     % Summary statistics
     meanOfTrials = mean(allData.(headers{dof}), 2);
     standardDev = std(allData.(headers{dof}), 0, 2);
     variance = var(allData.(headers{dof}), 0, 2);
     maxResult = max(meanOfTrials);
     minResult = min(meanOfTrials);
     rangeResult = abs(maxResult-minResult);
     
	% Name for ID metrics
	if strcmp(filename,'inverse_dynamics.sto')
		ID_metrics.(subject_name).(condition_name).(headers{dof}).('summary') = [maxResult, minResult, rangeResult];
          ID_metrics.(subject_name).(condition_name).(headers{dof}).('sd') = standardDev;
          ID_metrics.(subject_name).(condition_name).(headers{dof}).('var') = variance;
          ID_metrics.(subject_name).(condition_name).(headers{dof}).('mean') = meanOfTrials;
          
          % Load IK metrics
          metricsPathIK=[elabDataFolder 'AnalysedData' filesep 'Angles' filesep];
          cd(metricsPathIK);
          if exist('IKMetrics.mat', 'var')
               load('IKMetrics.mat');
               
               % Calculate joint power and add it to structure
               ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER') = meanOfTrials...
                    .* IK_metrics.(subject_name).(condition_name).(headers{dof}).('ANGULAR_VEL');
               
               % Determine when joint power polarities are changing -  can
               % also do for joint work
               
               abovePositiveThreshold = (ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER') > 0);
               spanLocs = bwlabel(abovePositiveThreshold);   %identify contiguous ones
               spanLength = regionprops(abovePositiveThreshold, 'area');  %length of each span
               spanLength = [ spanLength.Area];
               goodSpans = find(spanLength>=5);   %get only spans of 5+ points
               positiveInSpans = ismember(spanLocs, goodSpans);  %indices of these spans
               %                positiveInSpans = find(ismember(spanLocs, goodSpans));  %indices of these spans
               
               % Check this code to make sure it works.
               positivePower = ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER')(positiveInSpans);
               negativePower = ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER')(~positiveInSpans);
               
               % Integrate power to obtain joint work
               ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_WORK') = trapz(dt, ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER'));
               ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER_POS') = positivePower;
               ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER_NEG') = negativePower;
               
               % Determine when joint work polarities are changing
               abovePositiveThreshold = (ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_WORK') > 0);
               spanLocs = bwlabel(abovePositiveThreshold);   %identify contiguous ones
               spanLength = regionprops(abovePositiveThreshold, 'area');  %length of each span
               spanLength = [ spanLength.Area];
               goodSpans = find(spanLength>=5);   %get only spans of 5+ points
               positiveInSpans = ismember(spanLocs, goodSpans);  %indices of these spans
               %                positiveInSpans = find(ismember(spanLocs, goodSpans));  %indices of these spans
               
               % Check this code to make sure it works.
               positiveWork = ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_WORK')(positiveInSpans);
               negativeWork = ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_WORK')(~positiveInSpans);
               
               ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_WORK_POS') = positiveWork;
               ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_WORK_NEG') = negativeWork;
               
          else 
               disp('IK data file missing, check to see if you need to process it first');
          end
          
	elseif strcmp(filename,'ik.mot')
		% Name for IK metrics
		IK_metrics.(subject_name).(condition_name).(headers{dof}).('max_min_range') = [maxResult, minResult, rangeResult];
          IK_metrics.(subject_name).(condition_name).(headers{dof}).('sd') = standardDev;
          IK_metrics.(subject_name).(condition_name).(headers{dof}).('var') = variance;
          IK_metrics.(subject_name).(condition_name).(headers{dof}).('mean') = meanOfTrials;	
          IK_metrics.(subject_name).(condition_name).(headers{dof}).('ANGULAR_VEL') = gradient(meanOfTrials,dt)*(pi/180);
	end
end

% Save files
if strcmp(filename,'inverse_dynamics.sto')
     save([metricsPath, condition_name, '_IDMetrics.mat'], 'ID_metrics');
elseif strcmp(filename,'ik.mot')
     save([metricsPath, condition_name, '_IKMetrics.mat'], 'IK_metrics'); 
end

