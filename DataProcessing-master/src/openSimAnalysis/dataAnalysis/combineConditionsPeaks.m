function dataArray = combineConditionsPeaks(metrics)
%Combine data from separate conditions into a single array
%   Take the load sharing ID and IK metrics and combine all of the
%   individual conditions into one array for the variable of interest. Then
%   outputs the array.

% Define variable names
%phaseNames = fieldnames(metrics);
dataArray = struct();

%for l = 1:length(phaseNames)
%variableNames = fieldnames(metrics.(phaseNames{l}));
variableNames = fieldnames(metrics);

% delete fieldnames without 'peak'
TF = contains(variableNames, {'peak', 'PERC', 'WORK', 'mean'});
variableNames(~TF)=[];

% Loop through variables collected
for i = 1:length(variableNames)
	
	conditions = fieldnames(metrics.(variableNames{i}));
	% 		PN = contains(variableNames{i}, 'POS_KNEE');
	% 		PH = contains(variableNames{i}, 'POS_HIP');
	% 		NH = contains(variableNames{i}, 'NEG_HIP');
	% 		NA = contains(variableNames{i}, 'NEG_ANKLE');
	
	% Loop through conditions
	for k = 1:length(conditions)
		
		numPart = length(metrics.(variableNames{i}));
		
		% Add Nan to last digit if length of variable is less than 20 to
		% enable adding of the double to the array
		if numPart < 20
			metrics.(variableNames{i}).(conditions{k})(1,20) = nan;
		end
		
		dataArray.(variableNames{i})(:,k) = metrics.(variableNames{i}).(conditions{k});
		% 			if PN == 1
		% 				% Assign condition to data array
		% 				dataArray.(variableNames{i})(:,k) = metrics.(variableNames{i}).(conditions{k}) + 8;
		%
		% 			elseif PH == 1
		% 				% Assign condition to data array
		% 				dataArray.(variableNames{i})(:,k) = metrics.(variableNames{i}).(conditions{k}) - 8;
		% 			elseif NH == 1
		% 				% Assign condition to data array
		% 				dataArray.(variableNames{i})(:,k) = metrics.(variableNames{i}).(conditions{k}) + 8;
		% 			elseif NA == 1
		% 				% Assign condition to data array
		% 				dataArray.(variableNames{i})(:,k) = metrics.(variableNames{i}).(conditions{k}) - 8;
		% 			else
		% 				% Assign condition to data array
		% 				dataArray.(variableNames{i})(:,k) = metrics.(variableNames{i}).(conditions{k});
		% 			end
	end
	%end
end