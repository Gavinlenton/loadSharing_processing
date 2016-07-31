function [newPathName, dynamicFolders, dynamicCropFolders, maxc3dFile_name,...
     sessionData, maxc3dFileOther, maxName] = initialiseForAnalysis(pname)
%Defines paths and files for processing of EMG and dynamic walking trials
%in LS data
%   Input the pathname to determine max EMG trial name, path names, and
%   dynamic folder names

% Choose the files for EMG normalisation
% SQUAT TRIAL = KFJC1/2

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
maxTrials = selectMaxTrials(dynamicFolders); % Only KneeFJC trials

% Specify name of max file
maxc3dFile_name = [maxTrials{1},'.c3d'];

% Specify name of other max file if it's in there.
if ~exist([pname, filesep, maxc3dFile_name], 'file')
     maxName = maxTrials{2};
     maxc3dFileOther = [];
     if ~exist([pname, filesep, maxName, '.c3d'], 'file')
         disp('No max trials exist in this acquisition');
         maxc3dFileOther = [];
     end
else
     maxName = maxTrials{1};
     maxc3dFileOther = maxTrials{2};
end

end

