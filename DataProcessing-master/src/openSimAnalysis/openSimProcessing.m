% Script to run IK, ID, and any further analysis on load sharing data

% IK, ID, muscle analysis, and SO taken from Batch OpenSim Processing Scripts (BOPS)
% Copyright (C) 2015 Alice Mantoan, Monica Reggiani
% <https://simtk.org/home/bops>
% Please acknowledge the authors of BOPS when using this script.

%% Define BasePath with dynamicElaboration outputs and BOPS folder

addpath('Z:\s2921887\Google Drive\LS_main_data_collection\DataProcessing-master\src\openSimAnalysis');
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

% Select dynElab folder
BasePath=uigetdir('..\..\..\', 'Select Elaborated Data Folder');
% Select BOPS-master folder
folderBOPS = uigetdir('.\', 'Select the BOPS processing folder');

% Define subject names
subsdir=dir(BasePath);
for xx=1:length(subsdir)
     isub = [subsdir(:).isdir]; %# returns logical vector of subdirectories
     subjectNames = {subsdir(isub).name}';
     subjectNames(ismember(subjectNames,{'.','..'})) = [];
end

% Loop through subjects

for nS = 1:length(subjectNames)
     
     % Choose subject folder here
     fName = uigetdir(BasePath, 'Select the Subject for analysis');
     
     % Then create var for session names
     SessionDirs = dir(fName);
     isub=[SessionDirs(:).isdir];
     sessionFolders={SessionDirs(isub).name}';
     sessionFolders(ismember(sessionFolders,{'.','..', 'ROM'}))=[]; % dynamic subject folders
     
     clearvars subsdir isub
     
     % Loop through sessions for subject
     for sD = 1:length(sessionFolders)
          
          sessionName = [fName, filesep, sessionFolders{sD}];
          
          % Create var for trials names
          trialsDirs = dir([sessionName, filesep, 'dynamicElaborations']);
          isub=[trialsDirs(:).isdir];
          trialsFolders={trialsDirs(isub).name}';
          trialsFolders(ismember(trialsFolders,{'.','..', 'KneeFJC1_Processed', 'KneeFJC2_Processed'}))=[]; % dynamic trials folders
          
          % Need to select a different model for each conditions, so this
          % needs to loop only for fast and slow walking trial
          % Loop through trials (should be fast and slow walking)
          for tD = 1:length(trialsFolders)
               
               % Define input directory
               inputDir = [fName, filesep, sessionFolders{sD}, filesep,...
                    'dynamicElaborations', filesep, trialsFolders{tD}];
               
               % Select the .osim model file for that condition
               [genericModelID,genericModelPath]=uigetfile([sessionName, filesep,...
                    'staticElaborations', filesep, '*.osim'],...
                    'Select scaled OpenSim model for analysis');
               model_file=[genericModelPath genericModelID];
               
               % Navigate to folder where BOPS is;
               cd([folderBOPS, filesep, 'src']);
               
               %% INVERSE KINEMATICS and INVERSE DYNAMICS
               
               % If statement to check files you're running are for ID as well (this is
               % because Matlab keeps crashing...)
               
               idOption = questdlg('Are you running ID on same trials as IK?',...
                    'ID trials', 'Yes', 'No', 'No IK trials before, or matlab crashed', 'Yes');
               
               if strcmp(idOption, 'Yes') == 1
               
                    % Run IK as per normal
                    [IKoutputDir, IKtrialsOutputDir, IKprocessedTrials]=InverseKinematics(inputDir, model_file);

                    % Plot IK results
                    IKmotFilename='ik.mot'; %our default name
                    [coordinates,Xaxislabel]=plotResults('IK', IKoutputDir, IKtrialsOutputDir, model_file, IKprocessedTrials, IKmotFilename);
                    
                    % Run ID on same trials as IK
                    [IDoutputDir, IDtrialsOutputDir, IDprocessedTrials]=InverseDynamics(inputDir, model_file, IKoutputDir, IKprocessedTrials);
          
               elseif strcmp(idOption, 'No') == 1
                    
                    % Run IK as per normal
                    [IKoutputDir, IKtrialsOutputDir, IKprocessedTrials]=InverseKinematics(inputDir, model_file);
                    
                     % Plot IK results
                    IKmotFilename='ik.mot'; %our default name
                    [coordinates,Xaxislabel]=plotResults('IK', IKoutputDir, IKtrialsOutputDir, model_file, IKprocessedTrials, IKmotFilename);
                    
                    % Run ID on different trials as IK
                    [IDoutputDir, IDtrialsOutputDir, IDprocessedTrials]=InverseDynamics(inputDir, model_file, IKoutputDir);
                    
               else
                    
                    % Matlab crashed or you just want to run only ID
                    [IDoutputDir, IDtrialsOutputDir, IDprocessedTrials]=InverseDynamics(inputDir, model_file);
                    
               end
               
               %% Plotting ID results
               
               IDfilename='inverse_dynamics.sto';
               
               if exist('coordinates','var') && exist('Xaxislabel','var')  %same X-axis label
                    plotResults('ID', IDoutputDir, IDtrialsOutputDir, model_file, IDprocessedTrials, IDfilename, coordinates, Xaxislabel);
               else if exist('coordinates','var')       %same coordinates, different X axis
                         plotResults('ID', IDoutputDir, IDtrialsOutputDir, model_file, IDprocessedTrials, IDfilename, coordinates); %add "_moment" to coordinates
                    else %if no IK before
                         plotResults('ID', IDoutputDir, IDtrialsOutputDir, model_file, IDprocessedTrials, IDfilename);
                    end
               end
               
               
               %% Muscle Analysis
%                
%                %when run IK before & want to process the same trials
%                [MAoutputDir,MAtrialsOutputDir, MAprocessedTrials]=MuscleAnalysis(inputDir, model_file, IKoutputDir, IKprocessedTrials);
%                
%                %when run IK before, but MA on different trials
%                %[MAoutputDir, MAtrialsOutputDir, MAprocessedTrials]=MuscleAnalysis(inputDir, model_file, IKoutputDir);
%                
%                %if no IK before:
%                %[MAoutputDir, MAtrialsOutputDir, MAprocessedTrials]=MuscleAnalysis(inputDir, model_file);
%                
%                
%                %% Plot Storage (MA)
%                if exist('Xaxislabel','var')
%                     plotStorage(Xaxislabel)
%                else
%                     plotStorage()
%                end
%                
%                
%                %% STATIC OPTIMIZATION
%                
%                %when run IK and ID before & want to process the same trials (of ID)
%                [SOoutputDir,SOtrialsOutputDir, SOprocessedTrials]=StaticOptimization(inputDir, model_file, IKoutputDir, IDoutputDir, IDprocessedTrials);
%                
%                %when run IK and ID before, but SO on different trials
%                %[SOoutputDir,SOtrialsOutputDir, SOprocessedTrials]=StaticOptimization(inputDir, model_file, IKoutputDir, IDoutputDir);
%                
%                %if no IK and ID before:
%                %[SOoutputDir,SOtrialsOutputDir, SOprocessedTrials]=StaticOptimization(inputDir, model_file);
%                
%                
%                %% Plot Storage (SO)
%                
%                if exist('Xaxislabel','var')
%                     plotStorage(Xaxislabel)
%                else
%                     plotStorage()
%                end
               
          end
     end
end