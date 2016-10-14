function processROMData(romDataFolder)
%Process the ROM data for load carriage trials to obtain means, sds, and
%CIs for all the conditions
%   Input path to the folder containing the ROM data

% Load the file containing all data
cd(romDataFolder)
load('romData.mat')

% Determine how many subject fields there are
subjects = fieldnames(anglesJointMeans);
conditions = fieldnames(anglesJointMeans.(subjects{1}));
allAnglesCalculated = fieldnames(anglesJointMeans.(subjects{1}).CRYE30);

% Create empty struct for results
anglesStatistics = struct();

% Loop through conditions
for j = 1:length(conditions)
	
	% Specify condition name
	conditionName = conditions{j};
	
	% Loop through subjects
	for i = 1:length(subjects)
		
		% Check if subject has that condition
		if any(strcmp(fieldnames(anglesJointMeans.(subjects{i})), conditionName))
			
			% Determine the names of the conditions which have ROM data for the
			% subject
			anglesCalculated = fieldnames(anglesJointMeans.(subjects{i}).(conditionName));
			
			% Check how many ROM trials were processed
			if length(anglesCalculated) == 4
				anglesData(i,1) = anglesJointMeans.(subjects{i}).(conditionName).(anglesCalculated{1});
				anglesData(i,2) = anglesJointMeans.(subjects{i}).(conditionName).(anglesCalculated{2});
				anglesData(i,3) = anglesJointMeans.(subjects{i}).(conditionName).(anglesCalculated{3});
				anglesData(i,4) = anglesJointMeans.(subjects{i}).(conditionName).(anglesCalculated{4});
				
				% If there's less than 4 then loop through them and
				% determine what trials were processed and assign to
				% appropriate variable
			elseif length(anglesCalculated) < 4
				
				% Loop through ROMs
				for k = 1:length(anglesCalculated)
					if strcmp(anglesCalculated{k}, 'HF')
						anglesData(i,1) = anglesJointMeans.(subjects{i}).(conditionName).(anglesCalculated{k});
					elseif strcmp(anglesCalculated{k}, 'ShoulderFF')
						anglesData(i,2) = anglesJointMeans.(subjects{i}).(conditionName).(anglesCalculated{k});
					elseif strcmp(anglesCalculated{k}, 'UUA')
						anglesData(i,3) = anglesJointMeans.(subjects{i}).(conditionName).(anglesCalculated{k});
					elseif  strcmp(anglesCalculated{k}, 'TF')
						anglesData(i,4) = anglesJointMeans.(subjects{i}).(conditionName).(anglesCalculated{k});
					end
				end
			end
		end
	end
	
	% Find when there are no values and convert them to NaNs
	[A, B] = find(anglesData(:,:) == 0);
	anglesData(A,B) = NaN;
	% Then assign all subjects data for that condition to the structure
	anglesStatistics.(conditionName) = anglesData;
	
	clearvars anglesData
	
end
	
% Save data
	save('romDataGrouped.mat', 'anglesStatistics');
	
	% Use grouped data to run repeated measures stats
	conditionNamez = fieldnames(anglesStatistics);
	
	for ii = 1:length(conditionNamez)
		HF(:, ii) = anglesStatistics.(conditionNamez{ii})(:,1);
		ShoulderFF(:, ii) = anglesStatistics.(conditionNamez{ii})(:,2);
		UUA(:, ii) = anglesStatistics.(conditionNamez{ii})(:,3);
		TF(:, ii) = anglesStatistics.(conditionNamez{ii})(:,4);
	end
	
	Participants = [1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20];
	
	% Create tables of the data
	tableHF = table(Participants, HF(:,1), HF(:,2), HF(:,3), HF(:,4), HF(:,5), HF(:,6)...
		, HF(:,7), HF(:,8), HF(:,9), HF(:,10), HF(:,11), HF(:,12), HF(:,13), 'VariableNames', {'Participants', conditionNamez{1}...
		conditionNamez{2},conditionNamez{3},conditionNamez{4},conditionNamez{5},conditionNamez{6},conditionNamez{7},...
		conditionNamez{8},conditionNamez{9},conditionNamez{10}, conditionNamez{11}, conditionNamez{12}, conditionNamez{13}});
	
	% Means and SD
	
	save('HF_table.mat', 'tableHF');
	save('ShoulderFF_table.mat', 'tableShoulderFF');
	save('UUA_table.mat', 'tableUUA');
	save('TF_table.mat', 'tableTF');
	
	
	tableShoulderFF = table(Participants, ShoulderFF(:,1), ShoulderFF(:,2), ShoulderFF(:,3), ShoulderFF(:,4), ShoulderFF(:,5), ShoulderFF(:,6)...
		, ShoulderFF(:,7), ShoulderFF(:,8), ShoulderFF(:,9), ShoulderFF(:,10), ShoulderFF(:,11), ShoulderFF(:,12), ShoulderFF(:,13), 'VariableNames', {'Participants', conditionNamez{1}...
		conditionNamez{2},conditionNamez{3},conditionNamez{4},conditionNamez{5},conditionNamez{6},conditionNamez{7},...
		conditionNamez{8},conditionNamez{9},conditionNamez{10}, conditionNamez{11}, conditionNamez{12}, conditionNamez{13}});
	
	tableUUA = table(Participants, UUA(:,1), UUA(:,2), UUA(:,3), UUA(:,4), UUA(:,5), UUA(:,6)...
		, UUA(:,7), UUA(:,8), UUA(:,9), UUA(:,10), UUA(:,11), UUA(:,12), UUA(:,13), 'VariableNames', {'Participants', conditionNamez{1}...
		conditionNamez{2},conditionNamez{3},conditionNamez{4},conditionNamez{5},conditionNamez{6},conditionNamez{7},...
		conditionNamez{8},conditionNamez{9},conditionNamez{10}, conditionNamez{11}, conditionNamez{12}, conditionNamez{13}});
	
	tableTF = table(Participants, TF(:,1), TF(:,2), TF(:,3), TF(:,4), TF(:,5), TF(:,6)...
		, TF(:,7), TF(:,8), TF(:,9), TF(:,10), TF(:,11), TF(:,12), TF(:,13), 'VariableNames', {'Participants', conditionNamez{1}...
		conditionNamez{2},conditionNamez{3},conditionNamez{4},conditionNamez{5},conditionNamez{6},conditionNamez{7},...
		conditionNamez{8},conditionNamez{9},conditionNamez{10}, conditionNamez{11}, conditionNamez{12}, conditionNamez{13}});
	
	
end


