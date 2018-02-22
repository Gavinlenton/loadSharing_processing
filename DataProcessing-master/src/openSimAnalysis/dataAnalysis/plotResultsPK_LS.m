function []=plotResultsPK_LS(resultsPath, elabDataFolder, trialsList, x, Yquantities, varargin)
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
	
	folder_trial = [resultsPath, filesep, trialsList{k} filesep];
	
	if ~exist(folder_trial, 'dir')
		
		disp(['no folder: ' folder_trial]);
		dirRemoved = dirRemoved + 1;
		continue
	end
		
	for pk_result = 1:length(Yquantities)
		
		% Define filename
		filename = Yquantities{pk_result};
		
		if exist([folder_trial, filename], 'file')
			
			% Import data
			
			file = load_sto_file([folder_trial, filename]);
			
			stateNames = fieldnames(file);
			
			% Extract data into cell. k corresponds to trial number and j
			% to coordinate (column 1 is X, 2 is Y, and 3
			% is Z).
			
			for j =1: length(stateNames)-1
				
				coordCol=stateNames{j+1};
				
				% Acceleration COM
				if pk_result == 1
					
					pk_acc{k,j} = file.(coordCol);
					
					% Velocity COM
				elseif pk_result == 2
					
					pk_vel{k,j} = file.(coordCol);
					
					% Displacement COM
				else
					pk_pos{k,j} = file.(coordCol);
				end
				
			end
			
			timeVector{k}= file.('time');
		else
			error(['No folder specified for: ', subject_name, ', ', condition_name]);
		end
	end
end

% Remove the deleted trials from the cell matrix
% Acc
[sizeR, sizeC] = size(pk_acc);
sizeR_deletedRows = sizeR - dirRemoved;

% Remove any nonzero
pk_accNew = pk_acc(~cellfun('isempty',pk_acc));
pk_accNewLength = length(pk_accNew);
difference = pk_accNewLength/sizeC - sizeR_deletedRows;

if difference ~= 0
     sizeR_deletedRows = sizeR_deletedRows + difference;
end

pk_acc = reshape(pk_accNew, [sizeR_deletedRows, sizeC]);

% Vel
[sizeR, sizeC] = size(pk_vel);
sizeR_deletedRows = sizeR - dirRemoved;
% Remove any nonzero
pk_velNew = pk_vel(~cellfun('isempty',pk_vel));
pk_velNewLength = length(pk_velNew);
difference = pk_velNewLength/sizeC - sizeR_deletedRows;

if difference ~= 0
     sizeR_deletedRows = sizeR_deletedRows + difference;
end

pk_vel = reshape(pk_velNew, [sizeR_deletedRows, sizeC]);

% Disp
[sizeR, sizeC] = size(pk_pos);
sizeR_deletedRows = sizeR - dirRemoved;

% Remove any nonzero
pk_posNew = pk_pos(~cellfun('isempty',pk_pos));
pk_posLength = length(pk_posNew);
difference = pk_posLength/sizeC - sizeR_deletedRows;

if difference ~= 0
     sizeR_deletedRows = sizeR_deletedRows + difference;
end

pk_pos = reshape(pk_posNew, [sizeR_deletedRows, sizeC]);

%Save data in mat format
y.acc = pk_acc; y.vel = pk_vel; y.pos = pk_pos;
save([resultsPath, '_pk_results'], 'y');

headers = fieldnames(y);
% Time in frames for angular velocity calc and filtering (change this is you didn't
% collect data at 100Hz)
dt = 1/100;

% Mass of person + armour - changes based on condition type
if contains(condition_name, '15')
	
	massSystem = subject_weight + 15;
	
elseif contains(condition_name, '30')
	massSystem= subject_weight + 30;
	
else
	massSystem = subject_weight;
end

% Changed based on whether subject name has 2 numbers or one
if length(subject_name) == 10
	elabFolder = elabDataFolder(1:end-11);
else
elabFolder = elabDataFolder(1:end-10);
end

% Load existing file
cd(elabDataFolder);
if exist([elabFolder, 'PK_metrics_all.mat'], 'file')
	load([elabFolder, 'PK_metrics_all.mat']);
else
	allData = [];
end
	
% Vertical stiffness calculation
% stiffnessVert = massSystem * vAcc/(vPosCurrent - vPosEquil);

% Where massSystem = person + armour mass, vAcc is vertical acceleration of
% COM, vPosCurrent is current position of COM, and vPosEquil is mean
% vertical position of the COM between max and min VPos in stance.


axes = {'X', 'Y', 'Z'};
% Put the metrics in a structure
for dof = 1:length(headers)
	
	% 	stepTime = [];
	% Loop through good trials
	for trial = 1:sizeR
		
		for axis = 1:3
			% Combine data into array - resample if necessary
			% If it doesn't equal 101 then resample
			if length(y.(headers{dof}){trial, axis}) ~= 101
				allData.(headers{dof}).(axes{axis})(:, trial) = lpfilter(resample(y.(headers{dof})...
					{trial, axis}, 101, length(y.(headers{dof}){trial, axis}), 0), 8, dt, 'butter');
			else
				allData.(headers{dof}).(axes{axis})(:, trial) = y.(headers{dof}){trial, dof};
			end
			
			% If it's the final trial then compute the mean
			if trial == sizeR
				
				% Get summary stats of all trials
				meanOfTrials = mean(allData.(headers{dof}).(axes{axis}), 2);
				[maxResult, maxIndex] = max(meanOfTrials);
				[minResult, minIndex] = min(meanOfTrials);
				rangeResult = abs(maxResult-minResult);
				
				PK_metrics.(subject_name).(condition_name).(headers{dof}).(axes{axis}).('max_min_range')...
					= [maxResult, minResult, rangeResult];
				PK_metrics.(subject_name).(condition_name).(headers{dof}).(axes{axis}).('sd')...
					= std(allData.(headers{dof}).(axes{axis}), 0, 2);
				PK_metrics.(subject_name).(condition_name).(headers{dof}).(axes{axis}).('var')...
					= var(allData.(headers{dof}).(axes{axis}), 0, 2);
				PK_metrics.(subject_name).(condition_name).(headers{dof}).(axes{axis}).('mean')...
					= meanOfTrials;
				
				% Determine vPosEquil
				if dof == 3 && axis == 2 % Means that data is positional
					
					PK_metrics.(subject_name).(condition_name).(headers{dof}).(axes{axis}).('vPosMean')...
						= mean(meanOfTrials(maxIndex: minIndex));
					PK_metrics.(subject_name).(condition_name).('range_COM_disp_vert') = rangeResult;
					
					% Correct max index if it occurs after min index
					if maxIndex>minIndex
						disp('Max index occurs after min index')
						[peaks, loc] = findpeaks(meanOfTrials, 'minPeakProminence', 0.03);
						maxIndex = loc(1);
						PK_metrics.(subject_name).(condition_name).('maxIndex')...
							= maxIndex;
						[peaks2, loc2] = findpeaks(meanOfTrials*-1, 'minPeakProminence', 0.03);
						
						% For some reason if displacement is small then
						% reduce peak prominence
						if isempty(loc2)
							[peaks2, loc2] = findpeaks(meanOfTrials*-1, 'minPeakProminence', 0.02);
						end
						
						minIndex = loc2(1);
						
						if maxIndex<minIndex
							disp('Correction worked, good job team!')
							PK_metrics.(subject_name).(condition_name).('minIndex')...
								= minIndex;
							PK_metrics.(subject_name).(condition_name).(headers{dof}).(axes{axis}).('vPosMean')...
								= mean(meanOfTrials(maxIndex: minIndex));
						else
							disp('correction failed')
						end
						
					else
						PK_metrics.(subject_name).(condition_name).('maxIndex')...
							= maxIndex;
						PK_metrics.(subject_name).(condition_name).('minIndex')...
							= minIndex;
					end
				end
			end
		end
	end
end

% Define variables first
data4stiffy = PK_metrics.(subject_name).(condition_name).(headers{3}).(axes{2}).mean;
vPosEquil = PK_metrics.(subject_name).(condition_name).(headers{3}).(axes{2}).('vPosMean');

% Get indexes to start and end of period of interest
maxIndex = PK_metrics.(subject_name).(condition_name).maxIndex;
minIndex = PK_metrics.(subject_name).(condition_name).minIndex;

% Filter and differentiate to get accelerations
vVel = gradient(PK_metrics.(subject_name).(condition_name).(headers{3}).(axes{2}).mean...
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

if isnan(meanVertStiff)
	error(['Vertical stiffness is NaN for ', subject_name, ', ', condition_name], 'Error in stiffness calculation');
end

PK_metrics.(subject_name).(condition_name).('jointStiffnessMean') = abs(meanVertStiff);
PK_metrics.(subject_name).(condition_name).('jointStiffness') = stiffnessVert;
save([elabFolder, 'PK_metrics_all.mat'], 'PK_metrics');

