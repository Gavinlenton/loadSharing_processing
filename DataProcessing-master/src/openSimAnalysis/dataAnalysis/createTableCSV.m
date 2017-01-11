function createTableCSV(metrics_for_stats, subjectNumber, conditions, savePath)
%createTableCSV takes the input joint moments structure and save certain
%parameters to a .csv file
%   Input a structure containing joint moment data, an array of the subject numbers,
%   and the condition names, and output the relevant peak moments. First 
%   convert the data into tables and then output the table in a .csv file.

% Save peaks .csv files
cd(savePath);

variables = fieldnames(metrics_for_stats);

TF = contains(variables, {'peak', 'diffTBAS'});
variables(TF)=[];

for k = 1:length(variables)
	variableName = variables{k};
	
	% Create table with data
	table = array2table([subjectNumber, metrics_for_stats.(variableName)], 'VariableNames', ['Participants', conditions]);
	
	% Sort first column so participants are in ascending order
	table_final = sortrows(table);
	
	% Write the table to a .csv file
	writetable(table_final, [variableName '.csv']);
	
end

cd ../