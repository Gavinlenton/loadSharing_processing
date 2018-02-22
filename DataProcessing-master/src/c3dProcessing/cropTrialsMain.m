function times = cropTrialsMain(pname, c3dFile_name, physFolderName, acqLS, dynamicCropFolders, data)
%Main function to crop trials based on right heel-strike times
%   Crop the walking trials into consecutive gait cycles from right
%   heel-strike to next right heel-strike.


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
%     for newFiles = 1:length(rightHS)-1
%         fileSource = [physFolderName, filesep, [c3dFile_name(1:end-4), ...
%             num2str(newFiles), '.c3d']];
%         try
%             % Copy file to Google Drive folder
%             copyfile(fileSource, pname)
%         catch
%             
%             % Switch to the next file
%             fileSource = [physFolderName, filesep, [c3dFile_name(1:end-4), ...
%                 num2str(newFiles+1), '.c3d']];
%             
%             % Try again, with modified file
%             try
%                 copyfile(fileSource, pname)
%                 
%             catch
%                 % Re-throw original error
%                 disp('Cannot copy file, stop code and check it');
%                 uiwait
%             end
%         end
%         
%         % Delete file from physical drive
%         delete(fileSource);
%         
%     end
    
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
else
    times = [];
end

clearvars -except times

end

