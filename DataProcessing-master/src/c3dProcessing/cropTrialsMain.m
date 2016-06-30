function [times] = cropTrialsMain(c3dFile_name, physFolderName, acqLS, dynamicCropFolders)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


          % Navigate to InputData folder on C: because I cannot write
          % new c3dfiles to Google Drive folder.
          cd(physFolderName);
          
          % Insert Events and crop trials into conscutive gait cycles
          % that start on right heel-strike
          % Only run this for walking trials, not kneeFJC or static trials.
          walkingTrial = strcmp(dynamicCropFolders, c3dFile_name(1:end-4));
          
          % Function to crop
          if any(walkingTrial) == 1
               [rightHS, rightTO] = cropTrials(acqLS, c3dFile_name,data);
          
          % Now have to copy the cropped c3d files from C: to Google Drive.
          % Loop through c3d files
          for newFiles = 1:length(rightHS)-1
               fileSource = [physFolderName, filesep, [c3dFile_name(1:end-4), ...
                    num2str(newFiles), '.c3d']];
               try
                    % Copy file to Google Drive folder
                    copyfile(fileSource, c3dFile_folder)
               catch
                    
                    % Switch to the next file
                    fileSource = [physFolderName, filesep, [c3dFile_name(1:end-4), ...
                         num2str(newFiles+1), '.c3d']];
                    
                    % Try again, with modified file
                    try
                         copyfile(fileSource, c3dFile_folder)
                         
                    catch
                         % Re-throw original error
                         uiwait
                    end
               end
          end
          
          % Create a times variable with right heel-strike and right toe-off
          % to use in emg analysis.
          % Make sure HS corresponds with TO in length otherwise error will
          % be thrown.
          if length(rightHS) > length (rightTO)
               rightHS(end) = [];
          elseif length(rightHS) < length (rightTO)
               rightTO(end) = [];
          end
          
          times = [rightHS, rightTO];
          
          end

end

