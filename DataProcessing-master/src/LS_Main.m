%% -------------------------------------------------------------%%

% This script will:
% - Generate the appropriate files to run an OpenSim simulation sequence for the Load Sharing Data
% - Process EMG data for eventual CEINMS analysis
% - Process ROM trials and output joint angles

% Please acknowledge Glen Lichtwark from the University of Queensland
% for the use of his OpenSim pipeline tools

% Written by Gavin Lenton June 2016
% gavin.lenton@griffithuni.edu.au

% ---------------------------------------------------------------%

clear;
clc;

% Specify folder directories - may need to modify default directories as I
% set this up to quickly navigate to my folders.

if ispc
     
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
     physFolder = uigetdir('/Users/s2921887/Documents/', 'Select the input data folder on your physical drive');
     
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

%% Loop through all sessions for chosen subject.

for i = 1:length(subjectFolders)
     
     % Initialise and define directories
     pname = [fName, filesep, subjectFolders{i}];
     c3dFile_folder = pname;
     c3dFiles=dir([c3dFile_folder,'\*.c3d']);
     txtFiles=dir([c3dFile_folder,'\*.txt']);
     physFolderName = [physFolder, filesep, subjectName, filesep, subjectFolders{i}];
     [sessionConditions] = conditionNames(c3dFiles);
     
     %% RUN ACQUISITION INTERFACE
     
     % If the acqusition xml does not exist then generate one
     if ~exist(fullfile(pname, 'acquisition.xml'), 'file')
          disp('acquisition.xml file does not exist, running AcquisitionInterface...');
          
          % Nav to file directory and run modified interface function
          cd([motoDir, filesep, 'src', filesep, 'AcquisitionInterface', filesep]);
          AcquisitionInterface_LS(subjectName, pname, subjectFolders{i,1});
          
     else
          % It exists so continue happily
          message = [sprintf('acquisition.xml already exist in folder: %s', pname), ',\n continuing with analysis...'];
          disp(message)
     end
     
     %% MERGE EMG IF COLLECTED INTO A TXT FILE
     
     % Prompt asking if EMG data needs merging.
     mergeEMG = questdlg('Do you need to merge the EMG data with the c3d?',...
          'EMG Merging', 'Yes', 'No', 'No');
     
     if strcmp(mergeEMG, 'Yes') == 1
          % Merge EMG data
          mergeEmgMain(c3dFiles, txtFiles, physFolderName)
          
     else
          % If not merging then continue with analysis
          disp('EMG not merged with c3d file, continuing with analysis...');
     end
     
     %% RUN MOtoNMS C3D2BTK AND RENAME EMG CHANNELS
     
     % Navigate to directory where function is
     cd([motoDir, filesep, 'src' filesep, 'C3D2MAT_btk']);
     % Run c3d2mat
     C3D2MAT(fName);
     % Replace emg analog labels
     replaceAnalogLabels(pname);
     
     %% LOAD AND PROCESS C3D FILES IN THE ACQUISITION SESSION
     
     % Define all required path names
     [newPathName, dynamicFolders, dynamicCropFolders, maxc3dFile_name,...
          sessionData, maxc3dFileOther, maxName] = initialiseForAnalysis(pname);
     
     %Create directory to store elaboratedData
     mkdir(newPathName, c3dFile_name(1:end-4));
     
     % Check if EMG was captured in the session
     % Pick first condition in the session
     dynamicTrialsName = dynamicCropFolders{1,1};
     emgCaptured = checkSessionEMG(subjectName, dynamicTrialsName(1:end-15));
     
     w = waitbar(0,'Processing your data, be patient!');
     
     % Loop through c3d files that aren't ROM trials
     for t_trial = 1:length(dynamicFolders)
          
          % Load the c3d file and data using btk
          c3dFile_name = [dynamicFolders{t_trial,1}, '.c3d'];
          acqLS = btkReadAcquisition([c3dFile_folder, filesep, c3dFile_name]);
          data = btk_loadc3d([c3dFile_folder, filesep, c3dFile_name]);
          
          %% --- CROP TRIALS --- %%
          [times] = cropTrialsMain(c3dFile_name, physFolderName, acqLS, dynamicCropFolders);
          
          %% --- EMG PROCESSING --- %%
          % Only run EMG processing if EMG was collected
          if  emgCaptured == 1
               
               emgProcessingMain(c3dFile_name, maxc3dFile_name, maxc3dFileOther, sessionData, maxName, motoDir)
               
               % If no EMG captured in the session
          else
               disp('No EMG captured in this armour condition, continuing with analysis...');
          end
          waitbar(t_trial/length(dynamicFolders));
     end
     close(w);
     
     %% --- PROCESSING CROPPED TRIALS FOR USE IN OPENSIM --- %%
     croppedTrialsProcessing(pname)
     
     %% --- ROM TRIALS PROCESSING --- %%
     ROMTrialsProcessing(pname, i, sessionConditions)
end
