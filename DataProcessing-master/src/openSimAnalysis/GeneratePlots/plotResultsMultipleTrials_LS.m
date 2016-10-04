function []=plotResultsMultipleTrials_LS(resultsPath, trialsList, filename,x, Yquantities, varargin)
% Function to plot results from multiple trials

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


%%

close all

subject_weight = varargin{1};
subject_name = regexprep(varargin{2}, ' ', '_');
condition_name = varargin{3};

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
	
	% Normalising moments to body weight
	if strcmp(filename,'inverse_dynamics.sto')
		results = file.data;
		results(:, 2:end) = results(:, 2:end)./subject_weight;
	end
	
	for j =1: length(coord_idx)
		
		coordCol=coord_idx(j);
		
		y{k,j} = results(:,coordCol);
		
	end
	
	timeVector{k}=getXaxis(x, results);
	
end

% Path names to save data
metricsPath = [resultsPath filesep];

if strcmp(filename,'inverse_dynamics.sto')
	figurePath=[resultsPath filesep 'Figures' filesep 'Torques' filesep];
	% Load ID file if it exists
	cd(metricsPath);
	if exist('IDMetrics.mat', 'var')
		load('IDMetrics.mat');
	end
else
	figurePath=[resultsPath filesep 'Figures' filesep];
	% Load IK file if it exists
	cd(metricsPath);
	if exist('IKMetrics.mat', 'var')
		load('IKMetrics.mat');
	end
end

% Make directories
if exist(figurePath,'dir') ~= 7
	mkdir(figurePath);
end

if exist(metricsPath,'dir') ~= 7
	mkdir(metricsPath);
end

%Save data in mat format
save([figurePath, 'plottedData'], 'y')


% Settings for plot
plotLabels=regexprep(Yquantities, '_', ' ');
legendLabels=regexprep(trialsList, '_', ' ');
cmap = colormap(parula(180));
%plotTitle = filename;

for k=1:size(y,1)
	
	plotColor = cmap(round(1+5.5*(k-1)),:);
	
	for j=1:size(y,2)
		
		h(j)=figure(j);
		
		plot(timeVector{k}, y{k,j},'Color',plotColor)
		hold on
		
		xlabel(x)
		ylabel([plotLabels(j)])
		warning off
		legend(legendLabels)
		%title(filename)
		
		saveas(h(j),[figurePath  Yquantities{j} '.fig'])
	end
end

% After deleting bad data calculate the means and save to a structure
headers = file.colheaders(:, 2:end);

% Put the metrics in a structure
for dof = 1:length(headers)
	
	% Loop through good trials
	for trial = 1:length(trialsList)
		
		% Extract max, min, and range from each data set
		allData.ensemble(trial, :) = results(:, 2:end);
		allData.maxResult(trial,:) = max(results(:, 2:end));
		allData.minResult(trial,:) = min(results(:, 2:end));
		allData.rangeResult(trial,:) = allData.maxResult(trial,:)-allData.minResult(trial,:);
	end
	
	% Name for ID metrics
	if strcmp(filename,'inverse_dynamics.sto')
		ID_metrics.(subject_name).(condition_name).(headers{dof}).max = mean(maxResult(:, dof));
		ID_metrics.(subject_name).(condition_name).(headers{dof}).min = mean(minResult(:, dof));
		ID_metrics.(subject_name).(condition_name).(headers{dof}).range = mean(rangeResult(:, dof));
		save([metricsPath, 'IDMetrics.mat'], 'ID_metrics');
		
	elseif strcmp(filename,'ik.mot')
		
		% Name for IK metrics
		IK_metrics.(subject_name).(condition_name).(headers{dof}).max = mean(maxResult(:, dof));
		IK_metrics.(subject_name).(condition_name).(headers{dof}).min = mean(minResult(:, dof));
		IK_metrics.(subject_name).(condition_name).(headers{dof}).range = mean(rangeResult(:, dof));
		save([metricsPath, 'IKMetrics.mat'], 'IK_metrics');
	end
end



