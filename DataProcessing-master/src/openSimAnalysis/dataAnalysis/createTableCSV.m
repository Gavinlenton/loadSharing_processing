function createTableCSV(metrics_for_stats, subjectNumber, conditions, savePath)
%createTableCSV takes the input joint moments structure and save certain
%parameters to a .csv file
%   Input a structure containing joint moment data, an array of the subject numbers,
%   and the condition names, and output the relevant peak moments. First
%   convert the data into tables and then output the table in a .csv file.

% Save peaks .csv files
cd(savePath);

if contains(savePath, 'powers')
	
	phases = fieldnames(metrics_for_stats);
	
	for tt = 1:length(phases)
		
		variables = fieldnames(metrics_for_stats.(phases{tt}));
		
		%TF = contains(variables, {'peak', 'diffTBAS'});
		%variables(TF)=[];
		
		for k = 1:length(variables)
			variableName = variables{k};
			
			% If it's diff from TBAS data
			VN = contains(variableName, 'diffTBAS');
			
			if VN ~= 1
				% Create table with data
				table = array2table([subjectNumber, metrics_for_stats.(phases{tt}).(variableName)], 'VariableNames', ['Participants', conditions]);
				
			else
				% Create table with modified variable names - removed NA labels
				conditions4TBAS = conditions;
				V4TBAS = contains(conditions4TBAS, {'TBAS', 'NA'});
				conditions4TBAS(V4TBAS) = [];
				table = array2table([subjectNumber, metrics_for_stats.(variableName)], 'VariableNames', ['Participants', conditions4TBAS]);
				
			end
			% Sort first column so participants are in ascending order
			table_final = sortrows(table);
			
			table_final_export = standardizeMissing(table_final, 0);
			
			if ~isdir([savePath, filesep, phases{tt}])
				mkdir([savePath, filesep, phases{tt}]);
			end
			
			cd([savePath, filesep, phases{tt}])
			
			% Write the table to a .csv file
			writetable(table_final_export, [variableName '.csv']);
			
		end
		
		cd ../
		
	end
	
else
	variables = fieldnames(metrics_for_stats);
	
	%TF = contains(variables, {'peak', 'diffTBAS'});
	%variables(TF)=[];
	
	for k = 1:length(variables)
		variableName = variables{k};
		
		% If it's not a scalar (i.e., waveform)
		sizeOfData = size(metrics_for_stats.(variables{k}));
		
		if length(sizeOfData) < 3
			
			% If it's diff from TBAS data
			VN = contains(variableName, 'diffTBAS');
			
			if VN ~= 1
				% Create table with data
				table = array2table([subjectNumber, metrics_for_stats.(variableName)], 'VariableNames', ['Participants', conditions]);
				
			else
				% Create table with modified variable names - removed NA labels
				conditions4TBAS = conditions;
				V4TBAS = contains(conditions4TBAS, {'TBAS', 'NA'});
				conditions4TBAS(V4TBAS) = [];
				table = array2table([subjectNumber, metrics_for_stats.(variableName)], 'VariableNames', ['Participants', conditions4TBAS]);
				
			end
			
			
			% Sort first column so participants are in ascending order
			table_final = sortrows(table);
			
			% Write the table to a .csv file
			writetable(table_final, [variableName '.csv']);
		end
	end
end

cd ../