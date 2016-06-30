function [newPathName, dynamicFolders, dynamicCropFolders, maxc3dFile_name,...
     sessionData, maxc3dFileOther, maxName] = initialiseForAnalysis(pname)
%Defines paths and files for processing of EMG and dynamic walking trials
%in LS data
%   Input the pathname to determine max EMG trial name, path names, and
%   dynamic folder names

% Prompt to choose the files for EMG normalisation
% SQUAT TRIAL = KFJC1/2

% UNCOMMENT TO ANALYSE MAX TRIALS FOR NORMALISATION

prompt = {'Enter squat trial for EMG max normalisation:'};
dlg_title = 'Trial for EMG normalisation'; num_lines = 1;
def = {'KneeFJC'};
maxName = inputdlg(prompt, dlg_title, num_lines, def);
maxName = [maxName{1}, '_Processed'];

% Specify name of max file
maxc3dFile_name = [maxName, '.c3d'];

% Specify name of other max file if it's in there.
if strcmp(maxc3dFile_name, 'KneeFJC1_Processed.c3d') == 1
     maxc3dFileOther = 'KneeFJC2_Processed.c3d';
else
     maxc3dFileOther = 'KneeFJC1_Processed.c3d';
end

% Path to eventually store elaboratedData
newPathName = [strrep(pname, 'InputData', 'ElaboratedData'), filesep, 'dynamicElaborations'];

% Specify sessionData and dynamic folders
sessionData = [newPathName(1:end-19), 'sessionData'];
sessionDirs = dir(sessionData);
isub=[sessionDirs(:).isdir];
dynamicFolders={sessionDirs(isub).name}';
dynamicCropFolders={sessionDirs(isub).name}';
dynamicFolders = selectWalkingTrials(dynamicFolders, 1); % dynamic subject folders
dynamicCropFolders = selectWalkingTrials(dynamicCropFolders, 0); % dynamic subject folders without KneeFJC

end

