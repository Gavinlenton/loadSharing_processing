function metrics_fromTBAS = trajectoriesDiffFromTBAS(metrics, conditions)
% Subtract the values from the armour conditions from TBAS values to
% determine how much they differ from the control condition

%   INPUT - metrics: structure containing the metrics from an analysis
%         - conditions: cell contain all the condition names

%   OUTPUT - structure containing the metrics subtracted from TBAS

% Initialize output var
metrics_fromTBAS = struct();

% Remove NA conditions
conditions(ismember(conditions,{'NA_fast', 'NA_slow'}))=[];

% Get the names of the metrics
metricNames = fieldnames(metrics);

% Loop through metrics
for mName = 1:length(metricNames)
	
	metricName = metricNames{mName};
	
	% Loop through conditions
	for cName = 1:length(conditions)
				
		% Define the condition name
		conditionName = conditions{cName};
		
		% Check to see if condition is not TBAS
		if ~strcmp(conditionName(1:4), 'TBAS')
			
			% Define TBAS condition name
			tbasName = regexprep(conditionName, conditionName(1:end-7), 'TBAS');
			
			% Subtract values from TBAS
			metrics_fromTBAS.(metricName).(conditionName).mean = (metrics.(metricName).(conditionName).mean)...
				- (metrics.(metricName).(tbasName).mean);
			
		else
			% If it is TBAS then skip
			disp('condition is TBAS so skipping');
		end
	end
	
end
