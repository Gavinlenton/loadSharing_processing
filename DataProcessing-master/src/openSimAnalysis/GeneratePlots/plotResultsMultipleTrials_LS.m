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

% Elab folder name
elab_folder_index = regexp(metricsPath, 'Subject');
elab_folder = metricsPath(1:elab_folder_index-1);

% Put the metrics in a structure
for dof = 1:length(headers)
     
     allData.(headers{dof}) = [];
     stepTime = [];
     % Loop through good trials
     for trial = 1:length(trialsList)
          
          % FILTER PROCESSED DATA
          % Determine average stride time
          if dof == 1
               stepTime(trial) = length(y{trial,dof});
          end
          
          % Combine data into array - resample if necessary
          % If it doesn't equal 101 then resample
          if length(y{trial, dof}) ~= 101
               allData.(headers{dof})(:, trial) = lpfilter(resample(y{trial, dof}, 101, length(y{trial, dof})), 8, dt, 'butter');
          else
               allData.(headers{dof})(:, trial) = lpfilter(y{trial,dof}, 8, dt, 'butter');
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
          ID_metrics.(subject_name).(condition_name).(headers{dof}).('max_min_range') = [maxResult, minResult, rangeResult];
          ID_metrics.(subject_name).(condition_name).(headers{dof}).('sd') = standardDev;
          ID_metrics.(subject_name).(condition_name).(headers{dof}).('var') = variance;
          ID_metrics.(subject_name).(condition_name).(headers{dof}).('mean') = meanOfTrials;
          
          % Only get step time once
          if dof == 1
               step_time = round(mean(stepTime), 0) * dt;
               ID_metrics.(subject_name).(condition_name).('step_time') = step_time;
          end
          
          % Load IK metrics
          cd(elab_folder);
          if exist('IKMetrics_all.mat', 'file')
               load('IKMetrics_all.mat');
               
               % Calculate joint power and add it to structure
               ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER') = (meanOfTrials...
                    .* IK_metrics.(subject_name).(condition_name).(headers{dof}(1:end-7)).('ANGULAR_VEL'));
               
               % Determine when joint power polarities are changing to
               % first get positive powers only
               
               abovePositiveThreshold = ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER') >= 0;
               spanLocs = bwlabel(abovePositiveThreshold);   %identify contiguous ones
               spanLength = regionprops(abovePositiveThreshold, 'area');  %length of each span
               spanLength = [ spanLength.Area];
               goodSpans = find(spanLength>=5);   %get only spans of 5+ frames
               positiveInSpans = ismember(spanLocs, goodSpans);  %indices of these spans
               ID_metrics.(subject_name).(condition_name).(headers{dof}).('positive_power_index') = positiveInSpans;
               
               % Negative periods of interest
               belowPositiveThreshold = ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER') < 0;
               spanLocsBelow = bwlabel(belowPositiveThreshold);   %identify contiguous ones
               spanLengthBelow = regionprops(belowPositiveThreshold, 'area');  %length of each span
               spanLengthBelow = [ spanLengthBelow.Area];
               goodSpansBelow = find(spanLengthBelow>=5);   %get only spans of 10+ frames
               negativeInSpans = ismember(spanLocsBelow, goodSpansBelow);  %indices of these spans
               ID_metrics.(subject_name).(condition_name).(headers{dof}).('negative_power_index') = negativeInSpans;
               
               % POWER BURSTS - define burst periods as described in Winter
               % 1987 text.
               
               % HIP
               if strcmp(headers{dof}, 'hip_flexion_r_moment')
                    
                    % Positive work
                    H1 = ismember(spanLocs, goodSpans(1));
                    work_H1 = trapz(ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER')(H1,1));
                    H3 = ismember(spanLocs, goodSpans(2));
                    work_H3 = trapz(ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER')(H3,1));
                    
                    % Sum positive work to get total positive work
                    work_pos_hip =  work_H1 + work_H3;
                    % Divide by sampling freq to get normalised work
                    work_pos_hip_norm = work_pos_hip*dt;
                    
                    % Multiply by 2 to get both limbs (assumed symmetry)
                    work_pos_hip_bothLimbs = work_pos_hip_norm *2;
                    power_pos_hip = work_pos_hip_bothLimbs / ID_metrics.(subject_name).(condition_name).('step_time');
                    
                    % Negative work
                    H2 = ismember(spanLocsBelow, goodSpansBelow(1));
                    work_neg_hip = trapz(ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER')(H2,1));
                    work_neg_hip_norm = work_neg_hip * dt;
                    
                    % Multiply by 2 to get both limbs (assumed symmetry)
                    work_neg_hip_bothLimbs = work_neg_hip_norm *2;
                    power_neg_hip = work_neg_hip_bothLimbs / ID_metrics.(subject_name).(condition_name).('step_time');
                    
                    % Net work
                    net_work_hip = work_pos_hip_norm + work_neg_hip_norm;
                    
                    % Add to structure
                    ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_WORK_NET') = net_work_hip;
                    ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_WORK_POS') = work_pos_hip_norm;
                    ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_WORK_NEG') = work_neg_hip_norm;
                    ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER_POS') = power_pos_hip;
                    ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER_NEG') = power_neg_hip;
                    
                    % KNEE
               elseif strcmp(headers{dof}, 'knee_angle_r_moment')
                    
                    % Negative work
                    K1 = ismember(spanLocsBelow, goodSpansBelow(1));
                    work_K1 = trapz(ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER')(K1,1));
                    
                    % Check to see how many spans were detected and assign
                    % appropriately.
                    if length(goodSpansBelow) == 4
                         K3 = ismember(spanLocsBelow, goodSpansBelow(3));
                         work_K3 = trapz(ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER')(K3,1));
                         K4 = ismember(spanLocsBelow, goodSpansBelow(4));
                         work_K4 = trapz(ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER')(K4,1));
                    
                    elseif length(goodSpansBelow) == 3
                         K3 = ismember(spanLocsBelow, goodSpansBelow(2));
                         work_K3 = trapz(ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER')(K3,1));
                         K4 = ismember(spanLocsBelow, goodSpansBelow(3));
                         work_K4 = trapz(ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER')(K4,1));            
                    end
                    
                    % Sum positive work to get total positive work
                    work_neg_knee =  work_K1 + work_K3 + work_K4;
                    work_neg_knee_norm = work_neg_knee*dt;
                    
                    % Multiply by 2 to get both limbs (assumed symmetry)
                    work_neg_knee_bothLimbs = work_neg_knee_norm *2;
                    power_neg_knee = work_neg_knee_bothLimbs / ID_metrics.(subject_name).(condition_name).('step_time');
                    
                    % Positive work 
                    A = find(goodSpans == 2, 1);
                    % If positive span was not detected (i.e., was less
                    % than 5 frames)
                    if isempty(A)
                         goodSpans = find(spanLength>=1);   %get spans of 1 frame or greater
                         positiveInSpans = ismember(spanLocs, goodSpans);  %indices of these spans
                         ID_metrics.(subject_name).(condition_name).(headers{dof}).('positive_power_index') = positiveInSpans;
                    end
                    
                     K2 = ismember(spanLocs, goodSpans(2));
                     b = find(K2 > 0);
                     
                     % If the assignment is in the range I'm looking for
                     % (e.g., 15-40% of stride) then perform calcs
                     if b(1) < 40
                          work_pos_knee = trapz(ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER')(K2,1));
                          work_pos_knee_norm = work_pos_knee * dt;
                     else
                          work_pos_knee_norm = 0;
                     end
                     
                    % Multiply by 2 to get both limbs (assumed symmetry)
                    work_pos_knee_bothLimbs = work_pos_knee_norm  *2;
                    power_pos_knee = work_pos_knee_bothLimbs / ID_metrics.(subject_name).(condition_name).('step_time');
                    
                    % Net work
                    net_work_knee = work_pos_knee_norm + work_neg_knee_norm;
                    
                    % Add to structure
                    ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_WORK_NET') = net_work_knee;
                    ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_WORK_POS') = work_pos_knee_norm;
                    ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_WORK_NEG') = work_neg_knee_norm;
                    ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER_POS') = power_pos_knee;
                    ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER_NEG') = power_neg_knee;
                    
                    % ANKLE
               elseif strcmp(headers{dof}, 'ankle_angle_r_moment')
                    
                    % Positive work
                    A2 = ismember(spanLocs, goodSpans(end-1));
                    c = find(A2 > 0);
                    
                    % Make sure it's the correct time stamp
                    if c(1) < 40
                         A2 = ismember(spanLocs, goodSpans(end));
                    elseif c(end) > 70
                         A2 = ismember(spanLocs, goodSpans(end-2));
                    end
                    
                    work_A2 = trapz(ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER')(A2,1));
                    
                    % Sum positive work to get total positive work
                    work_pos_ankle =  work_A2;
                    work_pos_ankle_norm = work_pos_ankle*dt;
                    
                    % Multiply by 2 to get both limbs (assumed symmetry)
                    work_pos_ankle_bothLimbs = work_pos_ankle_norm  *2;
                    power_pos_ankle = work_pos_ankle_bothLimbs / ID_metrics.(subject_name).(condition_name).('step_time');
                    
                    % Negative work
                    A1_half = ismember(spanLocsBelow, goodSpansBelow(1));
                    A1_half2 = ismember(spanLocsBelow, goodSpansBelow(2));
                    
                    % Again checking for correct timestamps
                    c = find(A1_half2 > 0);
                    if c(1) > 40
                         work_neg_ankle = trapz(ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER')(A1_half,1));
                    else
                    work_neg_ankle = trapz(ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER')(A1_half,1)) +...
                         trapz(ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER')(A1_half2,1));
                    end
                    work_neg_ankle_norm = work_neg_ankle * dt;
                    
                    % Multiply by 2 to get both limbs (assumed symmetry)
                    work_neg_ankle_bothLimbs = work_neg_ankle_norm  *2;
                    power_neg_ankle = work_neg_ankle_bothLimbs / ID_metrics.(subject_name).(condition_name).('step_time');
                    
                    % Net work
                    net_work_ankle = work_pos_ankle_norm + work_neg_ankle_norm;
                    
                    % Add to structure
                    ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_WORK_NET') = net_work_ankle;
                    ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_WORK_POS') = work_pos_ankle_norm;
                    ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_WORK_NEG') = work_neg_ankle_norm;
                    ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER_POS') = power_pos_ankle;
                    ID_metrics.(subject_name).(condition_name).(headers{dof}).('JOINT_POWER_NEG') = power_neg_ankle;
                    
               end
          else
               disp('IK data file missing, check to see if you need to process it first');
          end
          
     elseif strcmp(filename,'openKnee_ik.mot')
          
          % Load IK metrics if it exists
          cd(elab_folder);
          if exist('IKMetrics_all.mat', 'var')
               load('IKMetrics_all.mat');
          end
          
          % Name for IK metrics
          IK_metrics.(subject_name).(condition_name).(headers{dof}).('max_min_range') = [maxResult, minResult, rangeResult];
          IK_metrics.(subject_name).(condition_name).(headers{dof}).('sd') = standardDev;
          IK_metrics.(subject_name).(condition_name).(headers{dof}).('var') = variance;
          IK_metrics.(subject_name).(condition_name).(headers{dof}).('mean') = meanOfTrials;
          IK_metrics.(subject_name).(condition_name).(headers{dof}).('ANGULAR_VEL') = gradient(meanOfTrials,dt)*(pi/180);
     end
end

% Calculate total work and power metrics here
cd(elab_folder);
if strcmp(filename,'inverse_dynamics.sto')
     
     % Positive work at all joints
     work_pos_total = ID_metrics.(subject_name).(condition_name).('hip_flexion_r_moment').('JOINT_WORK_POS') + ...
          ID_metrics.(subject_name).(condition_name).('knee_angle_r_moment').('JOINT_WORK_POS') + ...
          ID_metrics.(subject_name).(condition_name).('ankle_angle_r_moment').('JOINT_WORK_POS');
     % Negative work at all joints
     work_neg_total = ID_metrics.(subject_name).(condition_name).('hip_flexion_r_moment').('JOINT_WORK_NEG') + ...
          ID_metrics.(subject_name).(condition_name).('knee_angle_r_moment').('JOINT_WORK_NEG') + ...
          ID_metrics.(subject_name).(condition_name).('ankle_angle_r_moment').('JOINT_WORK_NEG');
     % Positive power at all joints
     power_pos_total = ID_metrics.(subject_name).(condition_name).('hip_flexion_r_moment').('JOINT_POWER_POS') + ...
          ID_metrics.(subject_name).(condition_name).('knee_angle_r_moment').('JOINT_POWER_POS') + ...
          ID_metrics.(subject_name).(condition_name).('ankle_angle_r_moment').('JOINT_POWER_POS');
     % Negative power at all joints
     power_neg_total = ID_metrics.(subject_name).(condition_name).('hip_flexion_r_moment').('JOINT_POWER_NEG') + ...
          ID_metrics.(subject_name).(condition_name).('knee_angle_r_moment').('JOINT_POWER_NEG') + ...
          ID_metrics.(subject_name).(condition_name).('ankle_angle_r_moment').('JOINT_POWER_NEG');
     
     % Determine percent contribution of each joint to total power
     powerP_perc_hip = (ID_metrics.(subject_name).(condition_name).('hip_flexion_r_moment').('JOINT_POWER_POS') / power_pos_total) * 100;
     powerP_perc_knee = (ID_metrics.(subject_name).(condition_name).('knee_angle_r_moment').('JOINT_POWER_POS') / power_pos_total) * 100;
     powerP_perc_ankle = (ID_metrics.(subject_name).(condition_name).('ankle_angle_r_moment').('JOINT_POWER_POS') / power_pos_total) * 100;
     
     powerN_perc_hip = (ID_metrics.(subject_name).(condition_name).('hip_flexion_r_moment').('JOINT_POWER_NEG') / power_neg_total) * 100;
     powerN_perc_knee = (ID_metrics.(subject_name).(condition_name).('knee_angle_r_moment').('JOINT_POWER_NEG') / power_neg_total) * 100;
     powerN_perc_ankle = (ID_metrics.(subject_name).(condition_name).('ankle_angle_r_moment').('JOINT_POWER_NEG') / power_neg_total) * 100;
     
     % Add to structure
     ID_metrics.(subject_name).(condition_name).('POWER_POS_TOTAL') = work_pos_total;
     ID_metrics.(subject_name).(condition_name).('POWER_NEG_TOTAL') = work_neg_total;
     ID_metrics.(subject_name).(condition_name).('WORK_POS_TOTAL') = power_pos_total;
     ID_metrics.(subject_name).(condition_name).('WORK_NEG_TOTAL') = power_neg_total;
     ID_metrics.(subject_name).(condition_name).('PERC_POS_HIP') = powerP_perc_hip;
     ID_metrics.(subject_name).(condition_name).('PERC_POS_KNEE') = powerP_perc_knee;
     ID_metrics.(subject_name).(condition_name).('PERC_POS_ANKLE') = powerP_perc_ankle;
     ID_metrics.(subject_name).(condition_name).('PERC_NEG_HIP') = powerN_perc_hip;
     ID_metrics.(subject_name).(condition_name).('PERC_NEG_KNEE') = powerN_perc_knee;
     ID_metrics.(subject_name).(condition_name).('PERC_NEG_ANKLE') = powerN_perc_ankle;
     
end

% Save files
if strcmp(filename,'inverse_dynamics.sto')
     save([elab_folder, 'IDMetrics_all.mat'], 'ID_metrics');
elseif strcmp(filename,'openKnee_ik.mot')
     save([elab_folder, 'IKMetrics_all.mat'], 'IK_metrics');
end

