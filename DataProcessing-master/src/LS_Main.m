 %% -------------------------------------------------------------%%

% This script will:
% - Generate the appropriate files to run an OpenSim simulation sequence for the Load Sharing Data
% - Process EMG data for eventual CEINMS analysis
% - Process ROM trials and output joint angles

% Please acknowledge Glen Lichtwark from the University of Queensland
% for the use of his OpenSim pipeline tools

% Written by Gavin Lenton June 2016
% gavin.lenton@griffithuni.edu.au

 %% ---------------------------------------------------------------%
clear;
clc;

% Specify folder directories - may need to modify default directories as I
% set this up to quickly navigate to my folders.

if ismac
    
    % Choose subject folder here
    fName = uigetdir('/Users/s2921887/Google Drive/Load Sharing Main Data Collection/InputData', 'Select the Subject for analysis');
    % Select physical folder directory
    physFolder = uigetdir('/Users/s2921887/Documents/', 'Select the input data folder on your physical drive');
    % Select MOtoNMS directory
    motoDir = uigetdir('/Users/s2921887/Google Drive/Load Sharing Main Data Collection/',...
          'Select the folder corresponding to your MOtoNMS directory');
    
elseif ispc
    
    % Choose subject folder here
    fName = uigetdir('Z:\s2921887\Google Drive\Load Sharing Main Data Collection\InputData', 'Select the Subject for analysis');
    % Select physical folder directory
    physFolder = uigetdir('C:\Users\s2921887\Documents\', 'Select the input data folder on your physical drive');
    % Select MOtoNMS directory
    motoDir = uigetdir('Z:\s2921887\Google Drive\Load Sharing Main Data Collection\',...
          'Select the folder corresponding to your MOtoNMS directory');
    
end

% Create directory cell array of session dates for subject chosen
subjectDirs = dir(fName);
isub=[subjectDirs(:).isdir];
subjectFolders={subjectDirs(isub).name}';
subjectFolders(ismember(subjectFolders,{'.','..'}))=[]; % dynamic subject folders

subjectName = regexp(fName, 'Subject\s\d*', 'match');

%% Loop through all sessions for chosen subject.

for i = 1:length(subjectFolders)
     
     % Assign session name
     pname = [fName, filesep, subjectFolders{i}];
     
     % Initialise directories
     c3dFile_folder = pname;
     c3dFiles=dir([c3dFile_folder,'\*.c3d']);
     txtFiles=dir([c3dFile_folder,'\*.txt']);
     
     %Name of the conditions in the session
     isubby = [c3dFiles(:).bytes]';
     bigFiles = find(isubby > 1500000);
     nameOfConditions = c3dFiles(bigFiles);
     nameOfConditions(strncmp({nameOfConditions.name}, 'Static1', 7)) = [];
     nameOfConditions(strncmp({nameOfConditions.name}, 'KneeFJC', 7)) = [];
     
     sessionConditions = {};
     for n = 1:length(nameOfConditions)
          sessionName = nameOfConditions(n).name(1:end-19);
          sessionConditions = [sessionConditions, sessionName];
     end
     
     %% RUN ACQUISITION INTERFACE IF .XML DOES NOT EXIST THEN GENERATE ONE.
     
     if ~exist(fullfile(pname, 'acquisition.xml'), 'file')
          disp('acquisition.xml file does not exist, running AcquisitionInterface...');
          
          % Run Acquisition Interface for each session to determine subject name,
          % weight, and height for use in LinScale
          
          % Nav to file directory
          cd([motoDir, filesep, 'src', filesep, 'AcquisitionInterface', filesep]);
          AcquisitionInterface_LS(subjectName, pname, subjectFolders{i,1});
          
     else
          message = [sprintf('acquisition.xml already exist in folder: %s', pname), ',\n continuing with analysis...'];
          disp(message)
          
     end
     
     %% MERGE EMG IF COLLECTED INTO A TXT FILE
     
     % Prompt asking if EMG data needs merging.
     mergeEMG = questdlg('Do you need to merge the EMG data with the c3d?',...
          'EMG Merging', 'Yes', 'No', 'No');
     
     c3dFilesMerge=dir([c3dFile_folder,'\*.c3d']);
     
     physFolderName = [physFolder, filesep, subjectName, filesep, subjectFolders{i}];
     
     % Delete files I don't want to analyse
     c3dFilesMerge(strncmp({c3dFilesMerge.name}, 'HF', 2)) = [];
     c3dFilesMerge(strncmp({c3dFilesMerge.name}, 'Shoulder', 8)) = [];
     c3dFilesMerge(strncmp({c3dFilesMerge.name}, 'TF', 2)) = [];
     c3dFilesMerge(strncmp({c3dFilesMerge.name}, 'Static1', 7)) = [];
     c3dFilesMerge(strncmp({c3dFilesMerge.name}, 'UUA', 3)) = [];
     
     if strcmp(mergeEMG, 'Yes') == 1
          % Loop through all trials
          for trial = 1:length(txtFiles)
               
               % Define input file/s
               inputc3d = c3dFiles(trial,1).name(1:end-4);
               txtFile = txtFiles(trial,1).name(1:end-4);
               
               if isempty(txtFile) == 1;
                    
                    % If no text file exists
                    disp('no txt files exist for this trial');
                    
               else
                    
                    cd(pname);
                    % Ensure processing only occurs when .txt file
                    % corresponds with .c3d file.
                    if strcmp(inputc3d(1:end-10), txtFile) == 0
                         
                         % Find the correct c3d file to match txt file
                         t = struct2cell(c3dFiles)';
                         index = strfind(t(:,1), txtFile);
                         inputc3d = c3dFiles((find(not(cellfun('isempty', index)))),1).name(1:end-14);
                         
                         % Because trial names have _Processed at the end
                         newInputc3d = inputc3d;
                         
                         % Check to see if they are now the same
                         if strcmp(newInputc3d, txtFile) == 0
                              disp(['Text file (', txtFile, ') still does not match c3d File (', inputc3d, ')']);
                              disp('Please ensure ALL txt files are in the same folder as c3d files');
                              
                         else
                              % Run EMG analysis to insert .txt file EMG into c3d
                              emgAsciiAnalysis([newInputc3d, '_Processed.c3d'], [txtFile, '.txt'],...
                                   physFolderName, c3dFile_folder);
                         end
                    else
                         % If they are the same just run processing
                         emgAsciiAnalysis([inputc3d, '.c3d'], [txtFile, '.txt'],...
                              physFolderName, c3dFile_folder);
                    end
               end
          end
     else
          % If not merging then continue with analysis
          disp('EMG not merged with c3d file, continuing with analysis...');
     end
     
     %% RUN MOtoNMS C3D2BTK AND RENAME EMG CHANNELS
     
     % Navigate to directory where function is
     cd([motoDir, filesep, 'src' filesep, 'C3D2MAT_btk']);
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
          % Only run this if analogData exists
          if exist('AnalogData.mat', 'var')
               load('AnalogData.mat')
               AnalogData.Labels(1:11) = newNames(1:11);
               % Save file
               save('AnalogData.mat', 'AnalogData');
          end
          cd ..
     end
     
     %% LOAD AND PROCESS C3D FILES IN THE ACQUISITION SESSION
     
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
     
     % Check if EMG was captured in the session
     % Pick first condition in the session
     dynamicTrialsName = dynamicCropFolders{1,1};
     emgCaptured = checkSessionEMG(subjectName, dynamicTrialsName(1:end-15));
     
     w = waitbar(0,'Processing your data, be patient!');
     
     % Loop through c3d files that aren't ROM trials
     for t_trial = 1:length(dynamicFolders)
          
          % Load the c3d file using btk
          c3dFile_name = [dynamicFolders{t_trial,1}, '.c3d'];
          acqLS = btkReadAcquisition([c3dFile_folder, filesep, c3dFile_name]);
          
          % Load data from c3d using Glenn's function
          data = btk_loadc3d([c3dFile_folder, filesep, c3dFile_name]);
          
          %Create directory to store elaboratedData
          mkdir(newPathName, c3dFile_name(1:end-4));
 
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
          %% EMG PROCESSING
          
          % Only run EMG processing if EMG was collected
          if  emgCaptured == 1
               
               % Check to see if  trial will be used as maximum for normalisation.
               % Inline function to determine if string exists
               cellfind = @(string)(@(cell_contents)(strcmp(string, cell_contents)));
               cell_array = dynamicFolders;
               string = maxName;
               isMaxExist = cellfun(cellfind(string), cell_array);
              
               % Set isMax based on the trial being a max trial (or not)
               A = ismember(dynamicFolders(any(isMaxExist,2)), c3dFile_name(1:end-4));
               
               emgMaxCondition = [sessionData, filesep, maxc3dFile_name(1:end-4)];
               emgMaxFileLoc = [sessionData, filesep, maxc3dFile_name(1:end-4), filesep, 'maxEMG'];
               
               % Create folder to put max EMG if it doesn't exist already
               if ~isdir(emgMaxFileLoc)
               mkdir(emgMaxCondition, 'maxEMG');
               end
               
               % File to load emg max(if it exists)
               [Muscle,emgMax,Trial,Times,Frame] = importMaxEMGFile([emgMaxFileLoc, filesep, 'maxEmg.txt']);
               
               % --Check to see if EMG data is from txt file or from .c3d to
               % know if we need to apply a notch filter--
               
               % Initialise
               asciiNames = {'Subject 6', 'Subject 8', 'Subject 13', 'Subject 14',...
                    'Subject 15', 'Subject 16', 'Subject 17'};
               
               % Loop through subject names known have txt files
               for ii = 1:length(asciiNames)
                    k = regexp(pname, asciiNames);
               end
               
               tf = isempty(k);
               
               
               % Use isMax comparison to see if EMG processing will be performed
               if any(A(:) == 1)
                    isMax = 1;
                    
                    % Check to see if max file already exists
                    if ~exist(emgMaxFile, 'file')
                         % Load data for max trial so we can process this first
                         emgProcessingMaxLS('no', sessionData, maxc3dFile_name(1:end-4), motoDir);
                         disp('maximum trial finished processing');
                    else
                         disp('Maximum trial has already been analysed for this session, continuing with analysis...');
                    end
                    
                    
               else
                    isMax = 0;
                    
                    % Check to see if max file already exists
                    if ~exist(emgMaxFile, 'file')
                         
                         disp('EMG max does not exist, processing this max trial first...');
                         
                         if tf == 0
                              
                              % Load data for max trial so we can process this first
                              % without notch filter
                              emgProcessingMaxLS('no', sessionData, maxName, motoDir);
                              disp('Max trial done, loading for EMG processing...');
                               [Muscle,emgMax,Trial,Times,Frame] = importMaxEMGFile([emgMaxFileLoc, filesep, 'maxEmg.txt']);
                              emgProcessingLS('no', sessionData, times, c3dFile_name(1:end-4), emgMax, motoDir);
                              
                         else
                              
                              % Processing for EMG data collected directly in nexus, this includes a
                              % notch filter
                              emgProcessingMaxLS('yes', sessionData, maxName, motoDir);
                              disp('Max trial done, loading for EMG processing...');
                              [Muscle,emgMax,Trial,Times,Frame] = importMaxEMGFile([emgMaxFileLoc, filesep, 'maxEmg.txt']);
                              emgProcessingLS('yes', sessionData, times, c3dFile_name(1:end-4), emgMax, motoDir);
                              
                         end
                         
                         % If it's the other KneeFJC trial then skip this
                    elseif strcmp(c3dFile_name, maxc3dFileOther)
                         disp([c3dFile_name, ' is the other max trial, but we''re using ', maxc3dFile_name]);
                         
                    else
                         disp('Maximum trial exists, running EMG processing...');
                         % Load max trial data
                         [Muscle,emgMax,Trial,Times,Frame] = importMaxEMGFile([emgMaxFileLoc, filesep, 'maxEmg.txt']);
                         
                         if tf == 0
                              % Run EMG processing for .txt data
                              emgProcessingLS('no', sessionData, times, c3dFile_name(1:end-4), emgMax, motoDir);
                              
                         else
                              % Processing for EMG data collected directly in nexus, this includes a
                              % notch filter
                              emgProcessingLS('yes', sessionData, times, c3dFile_name(1:end-4), emgMax, motoDir);
                         end
                         
                         disp(' EMG Processing complete');
                    end
               end
               
               % If no EMG captured in the session
          else
               disp('No EMG captured in this armour condition, continuing with analysis...');
          end
          waitbar(t_trial/length(dynamicFolders));
     end
     close(w);
     
     %% PROCESSING CROPPED TRIALS FOR USE IN OPENSIM
     
     % Re-set folder as that chosen above to include new files
     croppedSessionDirs = dir([pname, '\*.c3d']);
     isub2=[croppedSessionDirs(:).bytes]';
     % Only include files above 2000000 bytes as these are walking trials
     a = find(isub2 < 2000000);
     
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
          f = figure('Name', fileName);
          plot(dataFinal.fp_data.Time(:), force_data2(:,2),...
               dataFinal.fp_data.Time(:), force_data2(:,8))
          xlabel('Time (s)'); ylabel('Force (N)'); title('Vertical GRF');
          legend('Right foot', 'Left foot');
          legend boxoff;
          
          h = uicontrol('Position',[0 0 200 40],'String','Continue',...
               'Callback','uiresume(gcbf)');
          fprintf('\nMot file printed, click continue if you''re happy with the output\n');
          uiwait(gcf, 5);
          close(f);
 
          % Save output for future use
          save([pname, filesep, fileName(1:end-4), 'Data.mat'], 'dataFinal');

          % Close vars and figure to save memory
          close(gcf);
          clearvars dataFinal force_data2
     end
     
     %% ROM TRIALS PROCESSING
     
     % Re-set folder to include the ROM trials only
     c3dFile_ROM = pname;
     c3dFilesROM=dir([c3dFile_ROM,'\*.c3d']);
     
     % Delete files I don't want to analyse
     c3dFilesForROM = {c3dFilesROM.name}';
     c3dFilesForROM = selectROMTrials(c3dFilesForROM, sessionConditions);
     
     % Output max, min, and range of joint angles
     [anglesJoint] = determineJointAngles(c3dFilesForROM, pname);
     
     % Save to .xml
     xmlFileName = [fname, '.xml'];
     
     % Save each session in separate tabs
     xmlwrite(xmlFileName, angleJoint, ['Session', num2str(i)]);
end
