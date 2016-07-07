function  croppedTrialsProcessing(pname, fName, motoDir)
%Evaluates the cropped load sharing trials
%   Input the name of directory containing the c3d files and evaluate
%   the files to generate .mot and .trc files for analysis in OpenSim.

  % Re-set folder as that chosen above to include new files
     croppedSessionDirs = dir([pname, filesep, '*.c3d']);
     isub2=[croppedSessionDirs(:).bytes]';
     % Only include files above 2000000 bytes as these are walking trials
     a = isub2 < 2000000;
     
     % Delete files I don't want to analyse
     c3dFilesCropped = {croppedSessionDirs(a).name}';
     c3dFilesCropped = selectWalkingTrials(c3dFilesCropped, 0);
     
     % Run c3d2mat again on cropped trials.
     % Navigate to directory where function is
     cd([motoDir, filesep, 'src' filesep, 'C3D2MAT_btk']);
     
     % Run modified c3d2mat
     C3D2MAT_cropped(fName, c3dFilesCropped, pname);
     
     %% Loop through gait cycle trials
     for croppedTrials = 1:length(c3dFilesCropped)
          
          fileName = c3dFilesCropped{croppedTrials,1};
          
          %Load the cropped acquisition
          data1 = btk_loadc3d([pname, filesep, fileName], 50);
          
          % Assign force to feet, stitch forces together, and output .trc
          % and .mot files for further analysis.
          [dataFinal, force_data2] = assignForceOutputTrcMot(data1);
          
          % Check to see if forces assigned correctly
%           f = figure('Name', fileName);
%           plot(dataFinal.fp_data.Time(:), force_data2(:,2),...
%                dataFinal.fp_data.Time(:), force_data2(:,8))
%           xlabel('Time (s)'); ylabel('Force (N)'); title('Vertical GRF');
%           legend('Right foot', 'Left foot');
%           legend boxoff;
%           
%           uicontrol('Position',[0 0 200 40],'String','Continue',...
%                'Callback','uiresume(gcbf)');
%           fprintf('\nMot file printed, click continue if you''re happy with the output\n');
%           uiwait(gcf, 5);
%           close(f);
 
          % Save output for future use
          outputDir = [pname, filesep, 'matData'];
          
          if ~isdir(outputDir)
              mkdir(outputDir);
          end
          
          save([outputDir, filesep, fileName(1:end-4), '.mat'], 'dataFinal');

          % Close vars and figure to save memory
          close(gcf);
          clearvars dataFinal force_data2 data1
     end

end
