% This script will create the appropriate files to run an OpenSim
% simulation sequence for the Load Sharing Data
clear;
clc;

%% MOVE EMG DATA FROM TXT FILES INTO C3D

%Select one of the c3d files in the folder where data need to be inserted
[fname, pname] = uigetfile('*.c3d', 'Select c3d file');

% Set folder as that chosen above
c3dFile_folder = pname;
c3dFiles=dir([c3dFile_folder,'\*.c3d']);
txtFiles=dir([c3dFile_folder,'\*.txt']);


% Make current folder same as above
cd(pname);

for trial = 1:length(c3dFiles)
     
     % Define input file/s
     inputc3d = c3dFiles(trial,1).name;
     txtFile = txtFiles(trial,1).name;
     % If statement to ensure processing only occurs when .txt file
     % corresponds with .c3d file.
     if inputc3d(1:4) == txtFile(8:11)
     % Run EMG analysis to insert .txt file EMG into c3d
     emgAsciiAnalysis(inputc3d, txtFile);
     else 
          disp('no way jose');
     end
     
end

%% Run MOtoNMS c3d2btk function on data
% Navigate to directory where function is
cd('C:\Users\s2921887\Documents\Load Sharing Main Data Collection\MOtoNMS v2.2\src\C3D2MAT_btk');
run C3D2MAT;
% Navigate to sessionData folder and select AnalogDataLabels file
[matFileName, matPathName] = uigetfile('*.mat', 'Select mat file');
cd(matPathName);
load(matFileName);

% Create cell array with new channel names and file names
newNames = {'TA', 'Channel2', 'MG', 'LG', 'Channel5', 'BF', 'VM', 'VL',...
     'RF', 'Sol', 'MH'};
trialNames = {'slow', 'fast'};

% Replace analog labels in sessionData and individual trials
for files = 1:length(trialNames)
     cd([matPathName,trialNames(files)]);
     load('AnalogData.mat')
     AnalogData.Labels(1:11) = newNames(1:11);
     % Save file
     save('AnalogData.mat', 'AnalogData.mat');
     cd ..\
end

%% Select one of the c3d files in the inputData folder
[fname, pname] = uigetfile('*.c3d', 'Select c3d file');

% Set folder as that chosen above
c3dFile_folder = pname;
c3dFiles=dir([c3dFile_folder,'\*.c3d']);

% Make current folder same as above
cd(pname);

%% Load the c3dfiles and crop into multiple gait cycles
for t_trial = 1:length(c3dFiles)
     
     % Load the c3d file using btk
     c3dFile_name = c3dFiles(t_trial,1).name;
     acqLS = btkReadAcquisition([c3dFile_folder '\' c3dFile_name]);
     data = btk_loadc3d([c3dFile_folder '\' c3dFile_name]);
     
     % Create directory and folder to eventually store elaboratedData
     newpathname = [strrep(pname, 'InputData', 'ElaboratedData'), 'dynamicElaborations'];
     mkdir(newpathname, c3dFile_name(1:end-4));
     
     % Insert Events and crop trials into gait cycles beginning with right
     % HS
     times = cropTrials(acqLS, c3dFile_name,data);
     
     % Use these times to analyse EMG data
     
     % Check to see if EMG data is from ASCII file or from .c3d
     % Cell array containing subject with ASCII data
     
     asciiNames = {'Subject 6', 'Subject 8', 'Subject 13', 'Subject 14',...
          'Subject 15', 'Subject 16', 'Subject 17', 'Subject 18',...
          'Subject 19'}; % Careful because Subjects 18 and 19 contain EMG data from both sources
     isASCII = [];
     
     for subjectName = 1:length(asciiNames)
          k = strfind(pname, asciiNames{subjectName});
          isASCII = [isASCII, k];
     end
     tf = isempty(k);
     
     % Processing for ASCII data (e.g., includes notch filter).
     if tf ==0
          
          
     else
          % Processing for EMG data collected directly in c3d file. Want to
          % extract EMG and put into .mat structure
          
% Replace analogLabels to be consistent with other code

          
     end
     
end

%% Re-set folder as that chosen above to include new files
c3dFile_folder = pname;
c3dFiles=dir([c3dFile_folder,'\*.c3d']);

% Delete files I don't want to analyse
c3dFiles(strncmp({c3dFiles.name}, 'static1', 5)) = [];
c3dFiles(strncmp({c3dFiles.name}, 'fast_Processed.c3d', 18)) = [];
c3dFiles(strncmp({c3dFiles.name}, 'slow_Processed.c3d', 18)) = [];

%% C3D TO TRC and modify GRF data

for croppedTrials = 1:length(c3dFiles)
     fileName = c3dFiles(croppedTrials,1).name;
     dataFinal = assignForceOutputTrcMot(fileName, pname);
     % MOtoNMS function to print the EMG file
     % Need to modify this so it only analyses from trial two, and adds 5ms
     % of padding to the start of the crop (i.e., begins before previous
     % gait phase ends
     printEMGmot(newpathname,data.time,EMGsData,EMGsLabels, tag);
end

%% Run BOPS for IK

% Navigate to folder where BOPS is
cd('C:\Users\s2921887\Documents\Load Sharing Main Data Collection\BOPS-master\src');

% Test ID

%% Run BOPS gui
run BOPSgui.m;

%% Run BOPS for ID
