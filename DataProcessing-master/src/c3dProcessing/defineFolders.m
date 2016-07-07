function [fName, physFolder, motoDir, subjectFolders, subjectName] = defineFolders(systemUsed)
%DEtermines the folder structure for analysis of load sharing data based
%upon the system used (i.e., mac or pc or linux)
%  Input the system used and the code will define directories for future
%  analysis.

if systemUsed == 1
     
     % Choose subject folder here
     fName = uigetdir('Z:\s2921887\Google Drive\Load Sharing Main Data Collection\InputData', 'Select the Subject for analysis');
     
     % Select physical folder directory
     physFolder = uigetdir('C:\Users\s2921887\Documents\', 'Select the input data folder on your physical drive');
     
     % Auto defines motonms directory (may need to change if your folder
     % structure is different - but it shouldn't be)
     motoDir = [fName(1:(regexp(fName, '\WInput'))), 'DataProcessing-master', filesep...
          'src', filesep, 'c3dProcessing', filesep, 'MOtoNMS-master'];
     
else
     % Choose subject folder here
     fName = uigetdir('/Users/s2921887/Google Drive/Load Sharing Main Data Collection/InputData', 'Select the Subject for analysis');
     
     % Select physical folder directory
     physFolder = uigetdir('/Users/s2921887/Documents/PhD_Griffith/', 'Select the input data folder on your physical drive');
     
     % Auto defines motonms directory (may need to change if your folder
     % structure is different - but it shouldn't be)
     motoDir = [fName(1:(regexp(fName, '\WInput'))), 'DataProcessing-master', filesep...
          'src', filesep, 'c3dProcessing', filesep, 'MOtoNMS-master'];
     
end

% Create directory cell array of session dates for chosen subject
subjectDirs = dir(fName);
isub=[subjectDirs(:).isdir];
subjectFolders={subjectDirs(isub).name}';
subjectFolders(ismember(subjectFolders,{'.','..'}))=[]; % dynamic subject folders
% Subject name
subjectName = regexp(fName, 'Subject\s\d*', 'match');

end

