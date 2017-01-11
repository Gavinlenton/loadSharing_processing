function dataArray = combineConditionsPeaks(metrics)
%Combine data from separate conditions into a single array
%   Take the load sharing ID and IK metrics and combine all of the
%   individual conditions into one array for the variable of interest. Then
%   outputs the array.

% Define variable names
variableNames = fieldnames(metrics);

% delete fieldnames without 'peak'
TF = contains(variableNames, 'peak');
variableNames(~TF)=[];

dataArray = struct();

% Loop through variables collected
for i = 1:length(variableNames)
	
	conditions = fieldnames(metrics.(variableNames{i}));
	
	% Loop through conditions
	for k = 1:length(conditions)
		% Assign condition to data array
		dataArray.(variableNames{i})(:,i) = metrics.(variableNames{i}).(conditions{k});
	end
end


end

