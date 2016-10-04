function  croppedTrialsProcessing(pname, fName, motoDir)
%Evaluates the cropped load sharing trials
%   Input the name of directory containing the c3d files and evaluate
%   the files to generate .mot and .trc files for analysis in OpenSim.

% Re-set folder as that chosen above to include new files
croppedSessionDirs = dir([pname, filesep, '*.c3d']);
isub2=[croppedSessionDirs(:).bytes]';
% Only include files above 500000 bytes as these are walking trials
a = isub2 < 500000;

% Delete files I don't want to analyse
c3dFilesCropped = {croppedSessionDirs(a).name}';
c3dFilesCropped = selectWalkingTrials(c3dFilesCropped, 0);

% Run c3d2mat again on cropped trials
% Navigate to directory where function is
% cd([motoDir, filesep, 'src' filesep, 'C3D2MAT_btk']);

% Run modified c3d2mat
% C3D2MAT_cropped(fName, c3dFilesCropped, pname);

% 	 % Figure properties
% cmap = colormap(parula(350));
% legendLabels=regexprep(c3dFilesCropped, '_', ' ');

%% Loop through gait cycle trials
for croppedTrialNum = 1:length(c3dFilesCropped)
	
	fileName = c3dFilesCropped{croppedTrialNum,1};
	
	%Load the cropped acquisition
	data1 = btk_loadc3d([pname, filesep, fileName], 20);
	
	% Assign force to feet, stitch forces together, and output .trc
	% and .mot files for further analysis.
	
	[dataFinal, force_data2] = assignForceOutputTrcMot(data1, fileName);
	
	
	if ~isempty(force_data2)
		
		% Uncomment to check to see if forces assigned correctly
% 		plotColor = cmap(round(1+1.5*(croppedTrialNum-1)),:);
% 		
% 		plot(dataFinal.fp_data.Time(:), force_data2(:,2),...
% 			dataFinal.fp_data.Time(:), force_data2(:,8), 'Color', plotColor)
% 		hold on
% 		
% 		xlabel('Time (s)')
% 		ylabel('Force (N)')
% 		title('Vertical GRF')
% 		legend(legendLabels, 'Location', 'eastoutside')
% 		legend boxoff
		
		% Save output for future use
		outputDir = [pname, filesep, 'matData'];
		
		if ~isdir(outputDir)
			mkdir(outputDir);
		end
		
		save([outputDir, filesep, fileName(1:end-4), '.mat'], 'dataFinal');
		
	end
	% Close vars to save memory
	clearvars dataFinal force_data2 data1 fileName outputDir
end

end
