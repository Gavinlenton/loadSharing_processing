%% -------------------------------------------------------------%%

% Main processing script for load sharing data
% - Generate the appropriate files to run an OpenSim simulation sequence for the Load Sharing Data
% - Process EMG data and output excitations in MOT format
% - Process ROM trials and output joint angles

% Please acknowledge Glen Lichtwark from the University of Queensland and
% Alice Mantoan from Universita di Padua for the use of their OpenSim pipeline tools

% Written by Gavin Lenton June 2016
% gavin.lenton@griffithuni.edu.au

%% --------------------------------------------------------------- %%

clear; clc; close all;

% Specify folder directories based on system used - may need to modify default directories as I
% set this up to quickly navigate to my folders.

if ispc
     % For PC use
    [fName, physFolder, motoDir, subjectFolders, subjectName] = defineFolders(1);
	addpath('Z:\s2921887\Google Drive\LS_main_data_collection\DataProcessing-master\src');
	addpath('Z:\s2921887\Google Drive\LS_main_data_collection\DataProcessing-master\src\c3dProcessing');
	addpath('Z:\s2921887\Google Drive\LS_main_data_collection\DataProcessing-master\src\emgProcessing');
	addpath('Z:\s2921887\Google Drive\LS_main_data_collection\DataProcessing-master\src\matlabOpenSimPipeline_LS');
else
     % Mac or Linux
     [fName, physFolder, motoDir, subjectFolders, subjectName] = defineFolders(2);
end

%% --- MAIN DATA ANALYSIS --- %%

for i = 1:length(subjectFolders)
     
     % Initialise and define directories
     pname = [fName, filesep, subjectFolders{i}];
     c3dFile_folder = pname;
     c3dFiles=dir([c3dFile_folder, filesep, '*.c3d']);
     txtFiles=dir([c3dFile_folder, filesep, '*.txt']);
     physFolderName = [physFolder, filesep, subjectName{1}, filesep, subjectFolders{i}];
     rootDir = fName(1:regexp(fName, '\WI'));
     [sessionConditions] = conditionNames(c3dFiles);
     
     %% --- RUN ACQUISITION INTERFACE --- %%
     
     % If the acqusition xml does not exist then generate one
     if ~exist(fullfile(pname, 'acquisition.xml'), 'file') && i == 1
          disp('acquisition.xml file does not exist, running AcquisitionInterface...');
          
          % Nav to file directory and run modified interface function
          cd([motoDir, filesep, 'src', filesep, 'AcquisitionInterface', filesep]);
          AcquisitionInterface_LS(subjectName, pname, subjectFolders{i,1});
          
          % If the acquisition doesn't exist then just copy it from other
          % file and updated parameter    
     elseif ~exist(fullfile(pname, 'acquisition.xml'), 'file') && i ~= 1
          fileSource = [fName, filesep, subjectFolders{1}, filesep, 'acquisition.xml'];
          copyfile(fileSource, c3dFile_folder)
          cd(pname);
          acquisitionInfo=xml_read(fullfile(pname, 'acquisition.xml'));
          % Update acquisition date and save
          acquisitionInfo.AcquisitionDate = subjectFolders{i};
          xml_write('acquisition.xml' , acquisitionInfo);
          
     else
          % It exists so continue happily
         fprintf('acquisition.xml already exist in folder: %s, \n Continuing with analysis... \n', pname);
     end
     
     %% --- MERGE EMG IF COLLECTED INTO A TXT FILE --- %%
     
     % Search inputData folder for txt files
     if ~isempty(fieldnames(txtFiles))
          % Merge EMG data if txt files exist
          mergeEmgMain(pname, c3dFiles, txtFiles, physFolderName)
     else
          % If not merging then continue with analysis
          disp('No text files exist in this session, continuing with analysis...');
     end
     
     %% --- RUN MOtoNMS C3D2BTK AND RENAME EMG CHANNELS --- %%
     
     % Navigate to directory where function is
     cd([motoDir, filesep, 'src' filesep, 'C3D2MAT_btk']);
     % Run c3d2mat
     C3D2MAT(fName, c3dFiles, pname);
     % Replace emg analog labels
     replaceAnalogLabels(pname);
     
     %% --- LOAD AND PROCESS C3D FILES IN THE ACQUISITION SESSION --- %%
     
     % Define all required path names
     [newPathName, dynamicFolders, dynamicCropFolders, maxc3dFile_name,...
          sessionData, maxc3dFileOther, maxName] = initialiseForAnalysis(pname);
     
     % Check if EMG was captured in the session with first condition
     dynamicTrialsName = dynamicCropFolders{1,1};
     emgCaptured = checkSessionEMG(subjectName, dynamicTrialsName(1:end-15), rootDir);
     
     w = waitbar(0,'Processing your data, be patient!');
     
     %% --- LOOP THROUGH DYNAMIC TRIALS --- %%
     
     for t_trial = 1:length(dynamicFolders)
          
          % Load the c3d file and data using btk
          c3dFile_name = [dynamicFolders{t_trial,1}, '.c3d'];
          acqLS = btkReadAcquisition([c3dFile_folder, filesep, c3dFile_name]);
          data = btk_loadc3d([c3dFile_folder, filesep, c3dFile_name]);
          
          %Create directory to store elaboratedData
          if ~exist([newPathName, c3dFile_name(1:end-4)], 'dir')
          mkdir(newPathName, c3dFile_name(1:end-4));
          end
          
          %% --- CROP TRIALS --- %%
          
          times = cropTrialsMain(pname, c3dFile_name, physFolderName, acqLS, dynamicCropFolders, data);
          
          %% --- EMG PROCESSING --- %%
          
          % Only run EMG processing if EMG was collected
          if  emgCaptured == 1
               
               emgProcessingMain(pname, c3dFile_name, maxc3dFile_name,...
                    maxc3dFileOther, sessionData, maxName,...
                    motoDir, dynamicFolders, times);
               
               % If no EMG captured in the session
          else
               disp('No EMG captured in this armour condition, continuing with analysis...');
          end
          
          % Clear for next trial
          clearvars acqLS data
          waitbar(t_trial/length(dynamicFolders));
          
     end
     close(w);
     
     %% --- PROCESSING CROPPED TRIALS FOR USE IN OPENSIM --- %%
     
     croppedTrialsProcessing(pname, fName, motoDir);
     
     %% --- ROM TRIALS PROCESSING --- %%
     
     ROMTrialsProcessing(pname, sessionConditions, fName);
     
     % Clear for next session
     clearvars -except fName motoDir physFolder subjectFolders subjectName
end