% This script will create the appropriate files to run an OpenSim
% simulation sequence for the Load Sharing Data
clear;
clc;

% Choose subject folder here
fName = uigetdir('C:\Users\s2921887\Documents\Load Sharing Main Data Collection\InputData', 'Select the Subject for analysis');
subjectDirs = dir(fName);
isub=[subjectDirs(:).isdir];
subjectFolders={subjectDirs(isub).name}';
subjectFolders(ismember(subjectFolders,{'.','..'}))=[]; % dynamic subject folders

%% Loop through all sessions for that subject.

for i = 1:length(subjectFolders)
     
     % Assign session name
     pname = [fName, filesep, subjectFolders{i}];
     
     % Set folder as that chosen above
     c3dFile_folder = pname;
     c3dFiles=dir([c3dFile_folder,'\*.c3d']);
     txtFiles=dir([c3dFile_folder,'\*.txt']);
     
     % Make current folder same as above
     cd(pname);
     
     prompt = {'Do you need to merge the .txt file EMG with the c3d file?'};
     dlg_title = 'Yes or No'; num_lines = 1;
     def = {''};
     mergeEMG = inputdlg(prompt, dlg_title, num_lines, def);
     
     if strcmp(mergeEMG{1}, 'Yes') == 1
          
          for trial = 1:length(c3dFiles)
               
               % Define input file/s
               inputc3d = c3dFiles(trial,1).name(1:end-4);
               txtFile = txtFiles(trial,1).name(1:end-4);
               
               % If statement to ensure processing only occurs when .txt file
               % corresponds with .c3d file.
               if inputc3d(1:end-10) ~= txtFile
                    
                    % Find the correct c3d file to match .txt file
                    t = struct2cell(c3dFiles)';
                    index = strfind(t(:,1), txtFile);
                    inputc3d = c3dFiles((find(not(cellfun('isempty', index)))),1).name(1:end-14);
                    
                    % Run EMG analysis to insert .txt file EMG into c3d
                    try
                         emgAsciiAnalysis(inputc3d, txtFile);
                    catch me
                         error(['Text file (', txtFile, ') still does not match c3d File (', inputc3d(1:end-10), ')']);
                    end
                    
               else
                    emgAsciiAnalysis(inputc3d, txtFile);
               end
               
          end
     else
          disp('EMG not merged with c3d file, continuing with analysis');
     end
     
     % Run MOtoNMS c3d2btk function on data
     
     % Navigate to directory where function is
     cd('C:\Users\s2921887\Documents\Load Sharing Main Data Collection\MOtoNMS v2.2\src\C3D2MAT_btk');
     C3D2MAT(fName);
     % Navigate to sessionData folder and select AnalogDataLabels file
     [matFileName, matPathName] = uigetfile('*.mat', 'Select AnalogDataLabels.mat file');
     cd(matPathName);
     load(matFileName);
     
     % Create cell array with new channel names and file names
     newNames = {'TA', 'Channel2', 'MG', 'LG', 'Channel5', 'BF', 'VM', 'VL',...
          'RF', 'Sol', 'MH'};
     % Select names of trials to change channel names in (e.g., walking and KJC
     % trials)
     load('trialsName.mat');
     
     % Replace analog labels in sessionData and individual trials
     % sessionData first
     AnalogDataLabels(1:11) = newNames(1:11);
     save('AnalogDataLabels.mat', 'AnalogDataLabels');
     
     % Then individual trials
     for files = 1:length(trialsName)
          cd([matPathName,trialsName{files}]);
          load('AnalogData.mat')
          AnalogData.Labels(1:11) = newNames(1:11);
          % Save file
          save('AnalogData.mat', 'AnalogData');
          cd ..\
     end
     
     % Make current folder same as above
     cd(pname);
     
     % Initialise some variables first
     
     % Cell array containing subjects with ASCII data
     asciiNames = {'Subject 6', 'Subject 8', 'Subject 13', 'Subject 14',...
          'Subject 15', 'Subject 16', 'Subject 17', 'Subject 18',...
          'Subject 19'}; % Careful because Subjects 18 and 19 contain EMG data from both sources
     
     % Prompt to choose the files for EMG normalisation
     % SQUAT TRIAL = KFJC1/2
     
%      UNCOMMENT TO ANALYSE MAX TRIALS FOR NORMALISATION
     prompt = {'Enter squat trial for EMG max normalisation:',...
          'Enter name of sexiest girl in all of world:'};
     dlg_title = 'Trials for EMG normalisation'; num_lines = 1;
     def = {'KneeFJC', 'Salame'};
     maxNames = inputdlg(prompt, dlg_title, num_lines, def);
     
     % Load the c3dfiles and crop into multiple gait cycles - potential to expand this to analyse trials in a for loop.
     
     for t_trial = 1:length(c3dFiles)
          
          % Load the c3d file using btk
          c3dFile_name = c3dFiles(t_trial,1).name;
          acqLS = btkReadAcquisition([c3dFile_folder '\' c3dFile_name]);
          data = btk_loadc3d([c3dFile_folder '\' c3dFile_name]);
          
          % Create directory and folder to eventually store elaboratedData
          newPathName = [strrep(pname, 'InputData', 'ElaboratedData'), 'dynamicElaborations'];
          mkdir(newPathName, c3dFile_name(1:end-4));
          
          % Specify sessionData and dynamic folders
          sessionData = [newPathName(1:end-19), 'sessionData'];
          sessionDirs = dir(sessionData);
          isub=[sessionDirs(:).isdir];
          dynamicFolders={sessionDirs(isub).name}';
          dynamicCropFolders={sessionDirs(isub).name}';
          dynamicFolders(ismember(dynamicFolders,{'.','..'}))=[]; % dynamic subject folders
          dynamicCropFolders(ismember(dynamicCropFolders,{'.','..', 'KneeFJC2', 'KneeFJC1'}))=[];
          
          %      % Check to see if  trial will be used as maximum for normalisation.
          %      % UNCOMMENT TO ANALYSE WITH NORMALISATION
          
               isMaxExist = [];
          
               for trialName = 1:length(maxNames)
          
                    % Inline function to determine if string exists
                    cellfind = @(string)(@(cell_contents)(strcmp(string, cell_contents)));
                    cell_array = dynamicFolders;
                    string = maxNames{trialName};
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
          
          % Insert Events and crop trials into gait cycles beginning with right
          % HS. % Use Times as index in emgProcessing
          % Only run this for walking trials, not kneeFJC trials.
          walkingTrial = strcmp(dynamicCropFolders, c3dFile_name(1:end-4));
          
          if any(walkingTrial) == 1
               [rightHS, rightTO] = cropTrials(acqLS, c3dFile_name,data);
          else
          end
          
          times = [rightHS, rightTO];
          % Check to see if EMG data is from ASCII file or from .c3d
          % Initialise
          isASCII = [];
          for subjectName = 1:length(asciiNames)
               k = strfind(pname, asciiNames{subjectName});
               isASCII = [isASCII, k];
          end
          tf = isempty(k);
          
          % EMG PROCESSING
          
          % Check if notch filter should be applied.
          if tf == 0
               
               %Processing for ASCII data. isMax is always 0 until we decide to
               %process and analyse max data
               emgProcessingLS('no', sessionData, times, c3dFile_name(1:end-4), isMax);
               
          else
               % Processing for EMG data collected directly in nexus, this includes a
               % notch filter
               emgProcessingLS('yes', sessionData, times, c3dFile_name(1:end-4), isMax);
               
          end
          
     end
     
     % Re-set folder as that chosen above to include new files
     c3dFile_folderCropped = pname;
     c3dFilesCropped=dir([c3dFile_folderCropped,'\*.c3d']);
     
     % Delete files I don't want to analyse
     c3dFilesCropped(strncmp({c3dFilesCropped.name}, 'static1', 5)) = [];
     c3dFilesCropped(strncmp({c3dFilesCropped.name}, 'fast_Processed.c3d', 18)) = [];
     c3dFilesCropped(strncmp({c3dFilesCropped.name}, 'slow_Processed.c3d', 18)) = [];
     
     % C3D TO TRC and modify GRF data
     
     for croppedTrials = 1:length(c3dFiles)
          fileName = c3dFiles(croppedTrials,1).name;
          dataFinal = assignForceOutputTrcMot(fileName, pname);
     end
     
end
%% Run BOPS for IK
% Navigate to folder where BOPS is;
cd('C:\Users\s2921887\Documents\Load Sharing Main Data Collection\BOPS-master\src');
% Run BOPS gui
run BOPSgui.m;

%% Run BOPS for ID
