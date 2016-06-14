% This script will create the appropriate files to run an OpenSim
% simulation sequence for the Load Sharing Data

% Written by Gavin Lenton June 2016
% gavin.lenton@griffithuni.edu.au

clear;
clc;

% Choose subject folder here
fName = uigetdir('Z:\s2921887\Google Drive\Load Sharing Main Data Collection\InputData', 'Select the Subject for analysis');
subjectDirs = dir(fName);
isub=[subjectDirs(:).isdir];
subjectFolders={subjectDirs(isub).name}';
subjectFolders(ismember(subjectFolders,{'.','..'}))=[]; % dynamic subject folders
subjectName = fName(70:end);

%% Loop through all sessions for chosen subject.

for i = 1:length(subjectFolders)
     
     % Assign session name
     pname = [fName, filesep, subjectFolders{i}];
     
     % Initialise directories
     c3dFile_folder = pname;
     c3dFiles=dir([c3dFile_folder,'\*.c3d']);
     txtFiles=dir([c3dFile_folder,'\*.txt']);
     
     %% RUN ACQUISITION INTERFACE IF .XML DOES NOT EXIST
     % generate one.
     
     if ~exist(fullfile(pname, 'acquisition.xml'), 'file')
          disp('acquisition.xml file does not exist, running AcquisitionInterface...');
          
          % Run Acquisition Interface for each session to determine subject name,
          % weight, and height for use in LinScale
          % Select MOtoNMS directory
          motoDir = uigetdir('Z:\s2921887\Google Drive\Load Sharing Main Data Collection\',...
               'Select the folder corresponding to your MOtoNMS directory');
          
          % Nav to file directory
          cd([motoDir, '\src\AcquisitionInterface\']);
          AcquisitionInterface_LS(subjectName, pname, subjectFolders{i,1})
          
     else
          message = [sprintf('acquisition.xml already exist in folder: %s', pname), ',\n continuing with analysis...'];
          disp(message)
          
     end
     
     %% MERGE EMG IF COLLECTED INTO A TXT FILE
     
     % Prompt asking if EMG data needs merging.
     mergeEMG = questdlg('Do you need to merge the EMG data with the c3d?',...
          'EMG Merging', 'Yes', 'No', 'No');
     
     if strcmp(mergeEMG, 'Yes') == 1
          
          % Loop through all trials
          for trial = 1:length(c3dFiles)
               
               % Define input file/s
               inputc3d = c3dFiles(trial,1).name(1:end-4);
               txtFile = txtFiles(trial,1).name(1:end-4);
               
               if txtFile == [];
                    
                    % If no text file exists
                    disp('txt file does not exist for this trial');
                    
               else
                    
                    % Ensure processing only occurs when .txt file
                    % corresponds with .c3d file.
                    if inputc3d(1:end-10) ~= txtFile
                         
                         % Find the correct c3d file to match txt file
                         t = struct2cell(c3dFiles)';
                         index = strfind(t(:,1), txtFile);
                         inputc3d = c3dFiles((find(not(cellfun('isempty', index)))),1).name(1:end-14);
                         
                         % Run EMG analysis to insert .txt file EMG into c3d
                         try
                              emgAsciiAnalysis(inputc3d, txtFile);
                         catch me
                              error(['Text file (', txtFile, ') still does not match c3d File (', inputc3d(1:end-10), ')']);
                              disp('Please ensure ALL txt files are in the same folder as c3d files');
                              break
                         end
                         
                    else
                         % If they are the same just run processing
                         emgAsciiAnalysis(inputc3d, txtFile);
                    end
               end
          end
     else
          % If not merging then continue with analysis
          disp('EMG not merged with c3d file, continuing with analysis');
     end
     
     %% RUN MOtoNMS C3D2BTK AND RENAME EMG CHANNELS
     
     % Navigate to directory where function is
     cd([motoDir, 'src\C3D2MAT_btk']);
     % Run c3d2mat
     C3D2MAT(fName);
     
     % Navigate to sessionData folder and select AnalogDataLabels file
     cd([strrep(pname, 'InputData', 'ElaboratedData'), filesep, 'sessionData'])
     [matFileName, matPathName] = uigetfile('*.mat', 'Select AnalogDataLabels.mat file');
     load(matFileName);
     
     % Create cell array with new channel names and file names
     newNames = {'TA', 'Channel2', 'MG', 'LG', 'Channel5', 'BF', 'VM', 'VL',...
          'RF', 'Sol', 'MH'};
     
     % Load names of trials to change channel names
     load('trialsName.mat');
     
     % Replace analog labels in sessionData first
     AnalogDataLabels(1:11) = newNames(1:11);
     save('AnalogDataLabels.mat', 'AnalogDataLabels');
     
     % Then individual trials
     for files = 1:length(trialsName)
          cd(trialsName{files});
          load('AnalogData.mat')
          AnalogData.Labels(1:11) = newNames(1:11);
          % Save file
          save('AnalogData.mat', 'AnalogData');
          cd ..\
     end
     
     %% LOAD AND PROCESS C3D FILES IN THE ACQUISITION SESSION 
     
     % Prompt to choose the files for EMG normalisation
     % SQUAT TRIAL = KFJC1/2
     
     % UNCOMMENT TO ANALYSE MAX TRIALS FOR NORMALISATION
     prompt = {'Enter squat trial for EMG max normalisation:'};
     dlg_title = 'Trial for EMG normalisation'; num_lines = 1;
     def = {'KneeFJC'};
     maxName = inputdlg(prompt, dlg_title, num_lines, def);
     
     % Specify name of max file
     maxc3dFile_name = [maxName{1,1}, '.c3d'];
     
     % Loop through c3d files
     for t_trial = 1:length(c3dFiles)
          
          % Load the c3d file using btk
          c3dFile_name = c3dFiles(t_trial,1).name;
          acqLS = btkReadAcquisition([c3dFile_folder '\' c3dFile_name]);
          
          % Load data from c3d using Glenn's function
          data = btk_loadc3d([c3dFile_folder '\' c3dFile_name]);
          
          % Create directory to eventually store elaboratedData
          newPathName = [strrep(pname, 'InputData', 'ElaboratedData'), filesep, 'dynamicElaborations'];
          mkdir(newPathName, c3dFile_name(1:end-4));
          
          % Specify sessionData and dynamic folders
          sessionData = [newPathName(1:end-19), 'sessionData'];
          sessionDirs = dir(sessionData);
          isub=[sessionDirs(:).isdir];
          dynamicFolders={sessionDirs(isub).name}';
          dynamicCropFolders={sessionDirs(isub).name}';
          dynamicFolders(ismember(dynamicFolders,{'.','..'}))=[]; % dynamic subject folders
          dynamicCropFolders(ismember(dynamicCropFolders,{'.','..', 'KneeFJC2', 'KneeFJC1', 'static1'}))=[];
          
          % Navigate to InputData folder on C: because I cannot write
          % new c3dfiles to Google Drive folder.
          physFolder = uigetdir('C:\Users\s2921887\Documents\', 'Select the input data folder on your physical drive');
          physFolderName = [physFolder, filesep, subjectName, filesep, subjectFolders{i}];
          cd(physFolderName);
          
          % Insert Events and crop trials into conscutive gait cycles
          % that start on right heel-strike
          % Only run this for walking trials, not kneeFJC or static trials.
          walkingTrial = strcmp(dynamicCropFolders, c3dFile_name(1:end-4));
          
          % Function to crop
          if any(walkingTrial) == 1
               [rightHS, rightTO] = cropTrials(acqLS, c3dFile_name,data);
          end
          
          % Now have to copy the cropped c3d files from C: to Google Drive.
          % Loop through c3d files
          for newFiles = 1:length(rightHS)-1
               fileSource = [physFolderName, filesep, [c3dFile_name(1:end-4), ...
                    num2str(newFiles), '.c3d']];
               try
                    % Copy file to Google Drive folder
                    copyfile(fileSource, c3dFile_folder)
               catch error
                    
                    % Switch to the next file
                    fileSource = [physFolderName, filesep, [c3dFile_name(1:end-4), ...
                         num2str(newFiles+1), '.c3d']];
                    
                    % Try again, with modified file
                    try
                         copyfile(fileSource, c3dFile_folder)
                         
                    catch error2
                         % Re-throw original error
                         uiwait
                    end
               end
          end
          
          % Create a times variable with right hee-strike and right toe-off
          % to use in emg analysis.
          times = [rightHS, rightTO];
          
          %% EMG PROCESSING 
          %Check to see if  trial will be used as maximum for normalisation.
          isMaxExist = [];
          
          for trialName = 1:length(maxName)
               
               % Inline function to determine if string exists
               cellfind = @(string)(@(cell_contents)(strcmp(string, cell_contents)));
               cell_array = dynamicFolders;
               string = maxName{trialName};
               logicalCells = cellfun(cellfind(string), cell_array);
               
               isMaxExist = [isMaxExist, logicalCells];
          end
          
          % Set isMax based on the trial being a max trial (or not)
          A = ismember(dynamicFolders(any(isMaxExist,2)), c3dFile_name(1:end-4));
          
          if any(A(:) == 1)
               isMax = 1;
          else
               isMax = 0;
          end
          
          % --Check to see if EMG data is from txt file or from .c3d to
          % know if we need to apply a notch filter--
          % Initialise
          isASCII = [];
          
          % Loop through subject names known have txt files
          for subjectName = 1:length(asciiNames)
               k = strfind(pname, asciiNames{subjectName});
               isASCII = [isASCII, k];
          end
          tf = isempty(k);
          
          % Run EMG processing function. 
          % First check if notch filter should be applied. Only applied to data
          % collected directly into c3d
          if tf == 0
               
               %Processing for txt data.
               emgProcessingLS('no', sessionData, times, c3dFile_name(1:end-4), isMax, maxName{1,1});
               
          else
               % Processing for EMG data collected directly in nexus, this includes a
               % notch filter
               emgProcessingLS('yes', sessionData, times, c3dFile_name(1:end-4), isMax, maxName{1,1});
               
          end
          
     end
     
     %% PROCESSING CROPPED TRIALS FOR USE IN OPENSIM
     
     % Re-set folder as that chosen above to include new files
     c3dFile_folderCropped = pname;
     c3dFilesCropped=dir([c3dFile_folderCropped,'\*.c3d']);
     
     % Delete files I don't want to analyse
     c3dFilesCropped(strncmp({c3dFilesCropped.name}, 'static1', 5)) = [];
     c3dFilesCropped(strncmp({c3dFilesCropped.name}, 'fast_Processed.c3d', 18)) = [];
     c3dFilesCropped(strncmp({c3dFilesCropped.name}, 'slow_Processed.c3d', 18)) = [];
     c3dFilesCropped(strncmp({c3dFilesCropped.name}, 'KneeFJC1.c3d', 10)) = [];
     c3dFilesCropped(strncmp({c3dFilesCropped.name}, 'KneeFJC2.c3d', 10)) = [];
     
     % Loop through gait cycles
     for croppedTrials = 1:length(c3dFilesCropped)
          
          fileName = c3dFilesCropped(croppedTrials,1).name;
          
          %Load the new acquisition
          data1 = btk_loadc3d([c3dFile_folderCropped, filesep, fileName], 50);
          
          % Assign force to feet, stitch forces together, and output .trc
          % and .mot files for further analysis.
          dataFinal = assignForceOutputTrcMot(data1);
     end
end
