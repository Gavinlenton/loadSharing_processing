function []=plotResultsGRF_LS(resultsPath, elabDataFolder, trialsList, x, Yquantities, varargin)
% Function to plot results from GRF analysis

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
	motFile=dir([folder_trial, '*.mot']);
	
	% If there is not mot file then skip to next iteration
	if isempty(motFile) || contains(motFile.name, 'NFU')
		disp(['no .mot file in directory: ', folder_trial]);
		dirRemoved = dirRemoved + 1;
		continue
	else
		filename = motFile.name;
		
	end
	
	if exist([folder_trial, filename], 'file')
		
		% Import data
		file = load_sto_file([folder_trial, filename]);
		forceName = fieldnames(file);
		
		% Only extract 3 forces, 3 COPs and free (z) moment
		
		TF = contains(forceName, {'time','ground_force1', 'ground_torque1_y'});
		forceName(~TF)=[];
		
		for j =1: length(forceName)-1
			
			coordCol=forceName{j+1};
			
			grfData{k,j} = file.(coordCol);
			
		end
		
		timeVector{k}= file.('time');
	else
		disp(['no files in: ', filename ' in directory: ', folder_trial]);
	end
end

% Remove the deleted trials from the cell matrix
[sizeR, sizeC] = size(grfData);
sizeR_deletedRows = sizeR - dirRemoved;

% Remove any nonzero
grfNew = grfData(~cellfun('isempty',grfData));
grfNewLength = length(grfNew);
difference = grfNewLength/sizeC - sizeR_deletedRows;

if difference ~= 0
	sizeR_deletedRows = sizeR_deletedRows + difference;
end

grfData = reshape(grfNew, [sizeR_deletedRows, sizeC]);

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

% Changed based on whether subject name has 2 numbers or one
if length(subject_name) == 10
	elabFolder = elabDataFolder(1:end-11);
else
	elabFolder = elabDataFolder(1:end-10);
end

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
			upSample = lpfilter(resample(grfData...
				{trial, dof}, 1010, length(grfData{trial, dof}), 0), 8, dt, 'butter');
			% Then extract 101 points
			allData.(headers{dof})(:, trial) = upSample(indexArray);
			
		else
			allData.(headers{dof})(:, trial) = grfData.(headers{dof}){trial, dof};
		end
		
		% If it's the final trial then compute the mean
		if trial == trialNum
			
			% Get summary stats of all trials
			meanOfTrials = mean(allData.(headers{dof}), 2);
			[maxResult, maxIndex] = max(meanOfTrials);
			[minResult, minIndex] = min(meanOfTrials);
			rangeResult = abs(maxResult) + abs(minResult);
			
			GRF_metrics.(subject_name).(condition_name).(headers{dof}).('max_min_range')...
				= [maxResult, minResult, rangeResult];
			GRF_metrics.(subject_name).(condition_name).(headers{dof}).('sd')...
				= std(allData.(headers{dof}), 0, 2);
			GRF_metrics.(subject_name).(condition_name).(headers{dof}).('var')...
				= var(allData.(headers{dof}), 0, 2);
			GRF_metrics.(subject_name).(condition_name).(headers{dof}).('mean')...
				= meanOfTrials;
			
			% For vertical GRF do a check of the data
			if dof == 2 && maxResult < 9.81 * subject_weight
				disp(['Vertical GRF max is less than gravity X subject mass, something is wrong for: ', subject_name, ', ', condition_name]);
			end
			
			% Need to obtain values normalised to body mass
			GRF_metrics.(subject_name).(condition_name).(headers{dof}).('meanNormalised')...
				= meanOfTrials./subject_weight;
			GRF_metrics.(subject_name).(condition_name).(headers{dof}).('max_min_rangeNormalised')...
				= [max(meanOfTrials./subject_weight), min(meanOfTrials./subject_weight),...
				abs(max(meanOfTrials./subject_weight))+ abs(min(meanOfTrials./subject_weight))];
			
			
		end
	end
end

% Compute COM power to add to structure
% Individual limb method (Donelan et al. 2002) is the sum of the dot product of
% each force (e.g., GRF under each foot) with COM velocity.

GRFCombinedVector = [GRF_metrics.(subject_name).(condition_name).(headers{1}).mean, GRF_metrics.(subject_name).(condition_name).(headers{2}).mean,...
	GRF_metrics.(subject_name).(condition_name).(headers{3}).mean];

% Load existing point kinematics and GRF data data
cd(elabDataFolder);
load([elabFolder, 'PK_metrics_all.mat']);
load([elabFolder, 'IDmetrics_all.mat']);

% Load ID-derived powers and check against COM power
hip_power = (ID_metrics.(subject_name).(condition_name).hip_flexion_r_moment.JOINT_POWER +...
	ID_metrics.(subject_name).(condition_name).hip_adduction_r_moment.JOINT_POWER+...
	ID_metrics.(subject_name).(condition_name).hip_rotation_r_moment.JOINT_POWER) * subject_weight;
knee_power = (ID_metrics.(subject_name).(condition_name).knee_angle_r_moment.JOINT_POWER + ...
	ID_metrics.(subject_name).(condition_name).knee_adduction_r_moment.JOINT_POWER +...
	ID_metrics.(subject_name).(condition_name).knee_rotation_r_moment.JOINT_POWER) * subject_weight;
ankle_power = ID_metrics.(subject_name).(condition_name).ankle_angle_r_moment.JOINT_POWER * subject_weight;
lumbar_power = (ID_metrics.(subject_name).(condition_name).lumbar_extension_moment.JOINT_POWER +...
	ID_metrics.(subject_name).(condition_name).lumbar_bending_moment.JOINT_POWER +...
	ID_metrics.(subject_name).(condition_name).lumbar_rotation_moment.JOINT_POWER) * subject_weight * -1;

% Calculating only for one limb so probably have to double the output and
% assume symmetry.
g = 9.81;

velocityCOMVector = [PK_metrics.(subject_name).(condition_name).vel.X.mean, PK_metrics.(subject_name).(condition_name).vel.Y.mean,...
	PK_metrics.(subject_name).(condition_name).vel.Z.mean];

% Correct horizontal velocity because it has negative values - make it so
% average horizontal velocity = walking speed (e.g., 1.51 m/s)
offsetX = abs(min(PK_metrics.(subject_name).(condition_name).vel.X.mean));
correctedVelX = PK_metrics.(subject_name).(condition_name).vel.X.mean + offsetX;
meanVelX = mean(correctedVelX);

% Change movement speed value based on trial name
% If it ends with w then it's the sloW trial
if condition_name(end) == 'w'
	% Slow trial
	move_speed = 1.53; % 5.5 km/h
else
	% Fast trial
	move_speed = 1.81; % 6.5 km/h
end

% Final corrected velocity for movement speed
velXFinal = (correctedVelX + move_speed) - meanVelX;

% Obtain power values
powerCOMx = GRFCombinedVector(:,1) .* velXFinal;
powerCOMy = GRFCombinedVector(:,2) .* velocityCOMVector(:,2);
powerCOMz = GRFCombinedVector(:,3) .* velocityCOMVector(:,3);
powerCOM = powerCOMx + powerCOMy + powerCOMz;

% ID power
ID_power = hip_power + knee_power + ankle_power + lumbar_power;

GRF_metrics.(subject_name).(condition_name).('powerCOM') = powerCOM;

timeVec = (1:101)';

% Plot soft tissue power to check for consistency
% p = plot(timeVec, softTissuePower,  'r-', timeVec, ID_power, 'g--', timeVec, powerCOM,  'b:'); hold on; plot(timeVec, zeros(101), 'k-');
% set(p, 'LineWidth', 2); legend('Soft tissue', 'ID-based', 'COM-based', 'Location', 'SouthOutside',...
% 	'Orientation', 'horizontal'); set(gca, 'fontsize', 16, 'FontName', 'Gravity', 'Box', 'off');
% xlabel('% gait cycle'); ylabel('Mechanical Power (W)');

%% Get positive and negative regions of curve for work comparison over the
% stance phase. DO THIS FOR ID POWER AND COM POWER USING SAME TIME STAMPS

% Get time stamps of different regions from COM work

% Determine when joint power polarities are changing to first get positive
% powers only

abovePositiveThreshold = powerCOM(1:70) >= 0;
spanLocs = bwlabel(abovePositiveThreshold);   %identify contiguous ones
spanLength = regionprops(abovePositiveThreshold, 'area');  %length of each span
spanLength = [spanLength.Area];
goodSpans = find(spanLength>=5);   %get only spans of 5+ frames

% Rebound
positiveSpansRebound = ismember(spanLocs, goodSpans(1));  %indices of these spans
power_index_rebound = positiveSpansRebound;

% Toe off
% If there are not two power phases then skip this
if length(goodSpans) ~= 2
	powerCOM = powerCOM - 50;
	abovePositiveThreshold = powerCOM(1:70) >= 0;
	spanLocs = bwlabel(abovePositiveThreshold);   %identify contiguous ones
	spanLength = regionprops(abovePositiveThreshold, 'area');  %length of each span
	spanLength = [spanLength.Area];
	goodSpans = find(spanLength>=5);   %get only spans of 5+ frames
	
	% Rebound
	positiveSpansRebound = ismember(spanLocs, goodSpans(1));  %indices of these spans
	power_index_rebound = positiveSpansRebound;
	% Toe off
	positiveSpansToeOff = ismember(spanLocs, goodSpans(2));  %indices of these spans
	power_index_toeOff = positiveSpansToeOff;
	
else
	positiveSpansToeOff = ismember(spanLocs, goodSpans(2));  %indices of these spans
	power_index_toeOff = positiveSpansToeOff;
end

ID_toeOff_index = ID_metrics.(subject_name).(condition_name).ankle_angle_r_moment.positive_power_index;

% Negative periods of interest
belowPositiveThreshold = powerCOM(1:70) < 0;
spanLocsBelow = bwlabel(belowPositiveThreshold);   %identify contiguous ones
spanLengthBelow = regionprops(belowPositiveThreshold, 'area');  %length of each span
spanLengthBelow = [ spanLengthBelow.Area];
goodSpansBelow = find(spanLengthBelow>=5);   %get only spans of 5+ frames

% Collision
negativeSpansCollision = ismember(spanLocsBelow, goodSpansBelow(1));  %indices of these spans
power_index_collision = negativeSpansCollision;

% Make correction to ID_power data cause it's offset at the beginning
[mag, loc] = findpeaks(ID_power);

if loc(1) > 8
	disp('May have taken the wrong peak');
	
end

%plot(timeVec, ID_power); hold on; plot(loc, mag, 'o');

ID_power_corrected = ID_power;

if mag(1) == max(ID_power(1:10))
	ID_power_corrected(power_index_collision) = ID_power_corrected(power_index_collision) - (mag(1) +10) ;
	
elseif mag(2) == max(ID_power(1:10))
	ID_power_corrected(power_index_collision) = ID_power_corrected(power_index_collision) - (mag(2) +10) ;
	
end


% Compute soft tissue power as difference between COM power and joint power
softTissuePower = powerCOM - ID_power_corrected;

GRF_metrics.(subject_name).(condition_name).('powerID') = ID_power_corrected;
GRF_metrics.(subject_name).(condition_name).('powerSoftTissue') = softTissuePower;

plot(timeVec, ID_power, timeVec, ID_power_corrected);

% Pre load
negativeSpansPreLoad = ismember(spanLocsBelow, goodSpansBelow(2));  %indices of these spans
power_index_preLoad = negativeSpansPreLoad;

% Normalise curves to body weight and stride time, x2 to estimate output
% from both limbs
powerCOM_norm = (powerCOM/subject_weight)*ID_metrics.(subject_name).(condition_name).step_time*0.01*2;
powerID_norm = (ID_power_corrected/subject_weight)*ID_metrics.(subject_name).(condition_name).step_time*0.01*2;
softTissuePower = (softTissuePower/subject_weight)*ID_metrics.(subject_name).(condition_name).step_time*0.01*2;
%% Determine work in each region

% Region 1 = Collision
work_COM_collision =  trapz(powerCOM_norm(power_index_collision,1));
work_ID_collision =  trapz(powerID_norm(power_index_collision,1));
work_ST_collision =  work_COM_collision - work_ID_collision;

% Save to structure
GRF_metrics.(subject_name).(condition_name).('workCollisionID') = work_ID_collision;
GRF_metrics.(subject_name).(condition_name).('workCollisionCOM') = work_COM_collision;
GRF_metrics.(subject_name).(condition_name).('workCollisionST') = work_ST_collision;

% Region 2 = Rebound
work_COM_rebound =  trapz(powerCOM_norm(power_index_rebound,1));
work_ID_rebound =  trapz(powerID_norm(power_index_rebound,1));
work_ST_rebound =  work_COM_rebound - work_ID_rebound;

% Save to structure
GRF_metrics.(subject_name).(condition_name).('workReboundID') = work_ID_rebound;
GRF_metrics.(subject_name).(condition_name).('workReboundCOM') = work_COM_rebound;
GRF_metrics.(subject_name).(condition_name).('workReboundST') = work_ST_rebound;

% Region 3 = Preload
work_COM_preload =  trapz(powerCOM_norm(power_index_preLoad,1));
work_ID_preload =  trapz(powerID_norm(power_index_preLoad,1));
work_ST_preload =  work_COM_preload - work_ID_preload;

% Save to structure
GRF_metrics.(subject_name).(condition_name).('workPreloadID') = work_ID_preload;
GRF_metrics.(subject_name).(condition_name).('workPreloadCOM') = work_COM_preload;
GRF_metrics.(subject_name).(condition_name).('workPreloadST') = work_ST_preload;

% Region 4 = Push off
work_COM_toeOff =  trapz(powerCOM_norm(power_index_toeOff,1));
work_ID_toeOff =  trapz(powerID_norm(power_index_toeOff,1));
work_ID_toeOff_corrected = trapz(powerID_norm(ID_toeOff_index,1))+0.2;
work_ST_toeOff =  work_COM_toeOff - work_ID_toeOff_corrected;

% Save to structure
GRF_metrics.(subject_name).(condition_name).('workToeoffID') = work_ID_toeOff_corrected;
GRF_metrics.(subject_name).(condition_name).('workToeoffCOM') = work_COM_toeOff;
GRF_metrics.(subject_name).(condition_name).('workToeoffST') = work_ST_toeOff;

% Positive COM work and negative COM work
posCOMwork = work_COM_rebound + work_COM_toeOff; GRF_metrics.(subject_name).(condition_name).('posCOMwork') = posCOMwork;
negCOMwork = work_COM_collision + work_COM_preload; GRF_metrics.(subject_name).(condition_name).('negCOMwork') = negCOMwork;
netCOMwork = posCOMwork + negCOMwork; GRF_metrics.(subject_name).(condition_name).('netCOMwork') = netCOMwork;

% Can compare this value to positive power derived from

% Compute COM power from difference in kinetic and potential energy
% powerCOMeY = gradient(1/2*massSystem*((velocityCOMVector(:,2)*-1).^2) + subject_weight*g*PK_metrics.(subject_name).(condition_name).pos.Y.mean, 0.01 *(pi/180));

% Save structure
save([elabFolder, 'GRF_metrics_all.mat'], 'GRF_metrics');


