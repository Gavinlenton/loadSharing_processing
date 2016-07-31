function loadPlotEMG(sessionName, fName)
%Load and plot EMG linear envelopes for 9 muscles

% Folder where the EMG Results files are stored
dynamicElabFolder = [sessionName, filesep, 'dynamicElaborations'];
conditionDirs = dir(dynamicElabFolder);
isub=[conditionDirs(:).isdir];
conditionFolders={conditionDirs(isub).name}';
conditionFolders(ismember(conditionFolders,{'.','..'}))=[]; % dynamic subject folders

% Loop through conditions
for conditions = 1:length(conditionFolders)
	
	EMGresultsDir = [dynamicElabFolder, filesep, conditionFolders{conditions}, filesep, 'EMGs'];
	
	% Generate list of trials
	trials=dir([EMGresultsDir, filesep, '*.mot']);
	j=1;
	
	for k = 1:length(trials)
		trialsList{j}=trials(k).name;
		j = j + 1;
	end
	
	% Be selective if you want to
	% 	[trialsIndex,~] = listdlg('PromptString','Select trials to plot:',...
	% 		'SelectionMode','multiple',...
	% 		'ListString',trialsList);
	
	% 	inputTrials=trialsList(trialsIndex);
	
	inputTrials=trialsList;
	
	
	xaxislabel = '% Gait Cycle';
	
	% Plot here - these plots are for each muscle over many trials
	
	for i=1:length(inputTrials)
		
		try
			file=importdata([EMGresultsDir filesep inputTrials{i}]);
		catch me
			continue
		end
		
		Yquantities=file.colheaders(2:end); %take all columns except time
		coord_idx=[2:size(file.colheaders,2)];
		
		results = file.data;
		
		for j =1: length(coord_idx)
			
			coordCol=coord_idx(j);
			
			y{i,j} = results(:,coordCol);
			
		end
		
		timeVector{i}=getXaxis(xaxislabel, results);
		
	end
	
	% If file does not exist we know that it didn't load correctly
	if ~isempty(file) || exist(file, 'var')
		% Save the EMG cell as a mat file.
		armourIndex = regexp(EMGresultsDir, '\w*\d\d\w');
		
		% If armourIndex is less than 100 than I know it indexed wrong, which is
		% what will happen for the NA conditions.
		if armourIndex(end) < 100
			armourIndex = regexp(EMGresultsDir, 'NA');
		end
		
		% Save cell array of all trials
		saveName = [EMGresultsDir(armourIndex(end):end-14), 'EMGs.mat'];
		cd(EMGresultsDir);
		save(saveName, 'y');
		
		% Normalise to 101 data points
		NewSample = 101;
		
		%Resample to 101 points
		for i = 1:size(y,1)
			Orig_sample = size(y{i,1},1);
			TA(:,i) = resample(y{i,1}, NewSample, Orig_sample);
			MG(:,i) = resample(y{i,3}, NewSample, Orig_sample);
			LG(:,i) = resample(y{i,4}, NewSample, Orig_sample);
			Sol(:,i) = resample(y{i,10}, NewSample, Orig_sample);
			MH(:,i) = resample(y{i,11}, NewSample, Orig_sample);
			BF(:,i) = resample(y{i,6}, NewSample, Orig_sample);
			VM(:,i) = resample(y{i,7}, NewSample, Orig_sample);
			VL(:,i) = resample(y{i,8}, NewSample, Orig_sample);
			RF(:,i) = resample(y{i,9}, NewSample, Orig_sample);
		end
		
		muscleNames = {'TA', 'MG', 'LG', 'Sol', 'MH', 'BF', 'VM', 'VL', 'RF'};
		
		conditionName = EMGresultsDir(armourIndex(end):end-15);
		
		% DEFINE EMPTY STRUCTURE NAME (if it doesn't exist)
		cd(fName);
		if ~exist('EMGMetricsAll.mat', 'file')
			muscles.(conditionName) = struct();
			
			% If it does exist then load file
		else
			load('EMGMetricsAll.mat');
			
		end
		
		% Assign mean and SD to structure
		muscles.(conditionName).(muscleNames{1}) = struct('mean', mean(TA, 2), 'sd',  std(TA, 0, 2));
		muscles.(conditionName).(muscleNames{2}) = struct('mean', mean(MG, 2), 'sd',  std(MG, 0, 2));
		muscles.(conditionName).(muscleNames{3}) = struct('mean', mean(LG, 2), 'sd',  std(LG, 0, 2));
		muscles.(conditionName).(muscleNames{4}) = struct('mean', mean(Sol, 2), 'sd',  std(Sol, 0, 2));
		muscles.(conditionName).(muscleNames{5}) = struct('mean', mean(MH, 2), 'sd',  std(MH, 0, 2));
		muscles.(conditionName).(muscleNames{6}) = struct('mean', mean(BF, 2), 'sd',  std(BF, 0, 2));
		muscles.(conditionName).(muscleNames{7}) = struct('mean', mean(VM, 2), 'sd',  std(VM, 0, 2));
		muscles.(conditionName).(muscleNames{8}) = struct('mean', mean(VL, 2), 'sd',  std(VL, 0, 2));
		muscles.(conditionName).(muscleNames{9}) = struct('mean', mean(RF, 2), 'sd',  std(RF, 0, 2));
		
		% Save metrics of all trials
		save('EMGMetricsAll.mat', 'muscles');
		
		% 	% UNCOMMENT TO PLOT RESULTS
		% 	plotLabels=regexprep(Yquantities, '_', ' ');
		% 	legendLabels=regexprep(trialsList(1:end-4), '_', ' ');
		% 	cmap = colormap(parula(300));
		% 	%plotTitle = filename;
		%
		% 	for k=1:size(y,1)
		%
		% 		plotColor = cmap(round(1+5.5*(k-1)),:);
		%
		% 		for j=1:size(y,2)
		%
		% 			h(j)=figure(j);
		%
		% 			plot(timeVector{k}, y{k,j},'Color',plotColor)
		% 			hold on
		%
		% 			xlabel(xaxislabel)
		% 			ylabel([plotLabels(j)])
		% 			warning off
		% 			legend(legendLabels)
		% 		end
		% 	end
	end
	clearvars -except conditionFolders fName sessionName conditions dynamicElabFolder
end