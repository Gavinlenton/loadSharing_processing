function []=plotResultsGRF_LS(resultsPath, elabDataFolder, trialsList, x, Yquantities, varargin)
% Function to plot results from point kinematics analysis

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

dirRemoved = 0;

%Load data
for k=1:length(trialsList)
	
	folder_trial = [resultsPath filesep trialsList{k} filesep];
	
	% Define filename
	filename = [Yquantities{1}, num2str(k), '_grf.mot'];
	
	if exist([folder_trial, filename], 'file')
		
		% Import data
		file = load_sto_file([folder_trial, filename]);
		forceName = fieldnames(file);
		
		% Only extract 3 forces, 3 COPs and free (z) moment
		
		TF = contains(forceName, {'time','ground_force1', 'ground_torque1_z'});
		forceName(~TF)=[];
		
		%% CHECK THE FILENAMES FOR BOTH IK AND ID
		
		for j =1: length(forceName)-1
			
			coordCol=forceName{j+1};
			
			grfData{k,j} = file.(coordCol);
			
		end
		
		% 			file = load_sto_file([folder_trial, filename]);
		%
		% 			stateNames = fieldnames(file);
		%
		% 			% Extract data into cell. k corresponds to trial number and j
		% 			% to coordinate (column 1 is X, 2 is Y, and 3
		% 			% is Z).
		%
		% 			for j =1: length(stateNames)-1
		%
		% 				coordCol=stateNames{j+1};
		%
		% 				% Acceleration COM
		% 				if pk_result == 1
		%
		% 					pk_acc{k,j} = file.(coordCol);
		%
		% 					% Velocity COM
		% 				elseif pk_result == 2
		%
		% 					pk_vel{k,j} = file.(coordCol);
		%
		% 					% Displacement COM
		% 				else
		% 					pk_pos{k,j} = file.(coordCol);
		% 				end
		%
		% 			end
		
		timeVector{k}= file.('time');
	else
		rmdir([resultsPath, filesep, trialsList{k}], 's');
		dirRemoved = dirRemoved + 1;
	end
end

% % Remove the deleted trials from the cell matrix
% % Acc
% [sizeR, sizeC] = size(pk_acc);
% sizeR_deletedRows = sizeR - dirRemoved;
%
% % Remove any nonzero
% pk_accNew = pk_acc(~cellfun('isempty',pk_acc));
% pk_accNewLength = length(pk_accNew);
% difference = pk_accNewLength/length(Yquantities) - sizeR_deletedRows;
%
% if difference ~= 0
%      sizeR_deletedRows = sizeR_deletedRows + difference;
% end
%
% y = reshape(yNew, [sizeR_deletedRows, sizeC]);
%
% % Vel
% [sizeR, sizeC] = size(pk_vel);
% sizeR_deletedRows = sizeR - dirRemoved;
%
%
% % Disp
% [sizeR, sizeC] = size(pk_pos);
% sizeR_deletedRows = sizeR - dirRemoved;
%
% % Remove any nonzero
% yNew = y(~cellfun('isempty',y));
% yLength = length(yNew);
% difference = yLength/length(Yquantities) - sizeR_deletedRows;
%
% if difference ~= 0
%      sizeR_deletedRows = sizeR_deletedRows + difference;
% end
%
% y = reshape(yNew, [sizeR_deletedRows, sizeC]);

%Save data in mat format
save([resultsPath, '_grf_results'], 'grfData');

% Define some parameters for combining data
[trialNum, grfNum] = size(grfData);
headers = forceName(2:8);

% Frequency of GRF data (1000Hz)
dt = 1/1000;

% Mass of person + armour - changes based on condition type
if contains(condition_name, '15')
	
	massSystem = subject_weight + 15;
	
elseif contains(condition_name, '30')
	massSystem= subject_weight + 30;
	
else
	massSystem = subject_weight;
end

elabFolder = elabDataFolder(1:end-10);

% Load existing file
cd(elabDataFolder);
if exist([elabFolder, 'GRF_metrics_all.mat'], 'file')
	load([elabFolder, 'GRF_metrics_all.mat']);
else
	allData = [];
end

% Put the metrics in a structure
for dof = 1:grfNum
	
	% 	stepTime = [];
	% Loop through good trials
	for trial = 1:trialNum
		
		% Combine data into array - resample if necessary
		% If it doesn't equal 1010 then resample
		if length(grfData{trial, dof}) ~= 101
			
			% Get 101 points evenly spaced
			indexArray = floor(linspace(1,1010, 101));
			% Up sample data to 1010 points first
			test2 = lpfilter(resample(grfData...
				{trial, dof}, 1010, length(grfData{trial, dof}), 0), 8, dt, 'butter');
			% Then extract 101 points
			allData.(headers{dof})(:, trial) = test2(indexArray);
			
		else
			allData.(headers{dof})(:, trial) = grfData.(headers{dof}){trial, dof};
		end
		
		% If it's the final trial then compute the mean
		if trial == trialNum
			
			% Get summary stats of all trials
			meanOfTrials = mean(allData.(headers{dof}), 2);
			[maxResult, maxIndex] = max(meanOfTrials);
			[minResult, minIndex] = min(meanOfTrials);
			rangeResult = abs(maxResult-minResult);
			
			GRF_metrics.(subject_name).(condition_name).(headers{dof}).('max_min_range')...
				= [maxResult, minResult, rangeResult];
			GRF_metrics.(subject_name).(condition_name).(headers{dof}).('sd')...
				= std(allData.(headers{dof}), 0, 2);
			GRF_metrics.(subject_name).(condition_name).(headers{dof}).('var')...
				= var(allData.(headers{dof}), 0, 2);
			GRF_metrics.(subject_name).(condition_name).(headers{dof}).('mean')...
				= meanOfTrials;
		end
	end
end

% Define variables first
data4stiffy = GRF_metrics.(subject_name).(condition_name).(headers{3}).(axes{2}).mean;
vPosEquil = GRF_metrics.(subject_name).(condition_name).(headers{3}).(axes{2}).('vPosMean');

% Get indexes to start and end of period of interest
maxIndex = GRF_metrics.(subject_name).(condition_name).(headers{3}).(axes{2}).maxIndex;
minIndex = GRF_metrics.(subject_name).(condition_name).(headers{3}).(axes{2}).minIndex;

% Filter and differentiate to get accelerations
vVel = gradient(GRF_metrics.(subject_name).(condition_name).(headers{3}).(axes{2}).mean...
	,dt)*(pi/180);
vAcc = lpfilter((gradient(vVel...
	,dt)*(pi/180) *1000), 20, dt, 'damped'); % *1000 to account for discrepancy

% Initiate stiffness array
stiffnessVert = zeros(abs(minIndex-maxIndex), 1);
indexVert = 1;

% Loop through period of interest
for t = maxIndex:minIndex
	
	vPosCurrent = data4stiffy(t);
	
	stiffnessVert(indexVert, 1) = massSystem * vAcc(t)/(vPosCurrent - vPosEquil);
	indexVert = indexVert + 1;
end

% Vert stiffness save to person and condition
meanVertStiff = mean(stiffnessVert);

GRF_metrics.(subject_name).(condition_name).('jointStiffnessMean') = abs(meanVertStiff);
GRF_metrics.(subject_name).(condition_name).('jointStiffness') = stiffnessVert;
save([elabFolder, 'PK_metrics_all.mat'], 'PK_metrics');

