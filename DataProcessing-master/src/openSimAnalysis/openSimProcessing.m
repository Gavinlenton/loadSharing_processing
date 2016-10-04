%% Script to run IK, ID, and any further analysis on load sharing data %%

% IK, ID, muscle analysis, and SO taken from Batch OpenSim Processing Scripts (BOPS)
% Copyright (C) 2015 Alice Mantoan, Monica Reggiani
% <https://simtk.org/home/bops>
% Please acknowledge the authors of BOPS when using this script.

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath('..\');

%% Define BasePath with dynamicElaboration outputs and BOPS folder

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

%% Loop through subjects

for nS = 1:length(subjectNames)
     
     % Subject folder here
     fName = [BasePath, filesep, subjectNames{nS}];
     
     % Then create var for session names
     SessionDirs = dir(fName);
     isub=[SessionDirs(:).isdir];
     sessionFolders={SessionDirs(isub).name}';
     sessionFolders(ismember(sessionFolders,{'.','..', 'ROM'}))=[]; % dynamic subject folders
     
     clearvars subsdir isub
     
     %% Loop through sessions for subject
     for sD = 1:length(sessionFolders)
          
          sessionName = [fName, filesep, sessionFolders{sD}];
          
          % Create var for trials names
          trialsDirs = dir([sessionName, filesep, 'dynamicElaborations']);
          isub=[trialsDirs(:).isdir];
          trialsFolders={trialsDirs(isub).name}';
          trialsFolders(ismember(trialsFolders,{'.','..', 'KneeFJC1_Processed', 'KneeFJC2_Processed'}))=[]; % dynamic trials folders
          
          % Need to select a different model for each condition
          
          %% Loop through trials (should be fast and slow walking)
          for tD = 1:2:length(trialsFolders)
               
               % PUT SOMETHING HERE IF THERE IS ONLY A SLOW OR FAST WALKING
               % TRIAL FOR THE CONDITION
               
               % Define input directory
               inputDir = [fName, filesep, sessionFolders{sD}, filesep,...
                    'dynamicElaborations', filesep, trialsFolders{tD}];
               
               % Define model name for IK/ID and automatically find it
               modelName = [subjectNames{nS}, '_Scaled_Pose_OptMusc.osim'];
               
               % If it's NA trial then model name is in 'Processed' folder
               if regexp(trialsFolders{tD}, 'NA')
                    trialName = 'Processed';
                    
               else
                    % Otherwise it's in the folder with the conditions name
                    trialNameIndex = regexp(trialsFolders{tD}, '\d\w');
                    trialName = trialsFolders{tD}(1:trialNameIndex+1);
               end
               
               % Define model ID
               model_file = [sessionName, filesep, 'staticElaborations', filesep,...
                    ['Static1_' trialName], filesep, 'StaticCal', filesep, modelName];
               
               % Check if model file exist
               
               if exist(model_file, 'file')
                    
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
						 
						 %                          % Plot IK results
						 %                          IKmotFilename='ik.mot'; %our default name
						 %                          [coordinates,Xaxislabel]=plotResults('IK', IKoutputDir, IKtrialsOutputDir, model_file, IKprocessedTrials, IKmotFilename);
						 
						 % Run ID on same trials as IK
						 [IDoutputDir, IDtrialsOutputDir, IDprocessedTrials]=InverseDynamics(inputDir, model_file, sessionName, IKoutputDir, IKprocessedTrials);
						 
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
						 [IDoutputDir, IDtrialsOutputDir, IDprocessedTrials]=InverseDynamics(inputDir, model_file, sessionName);
						 
					 end
					
					%% ANALYZE TOOL
					
% 					filename='C:\OpenSim 3.1\Models\PendulumExample\analyzePointKinematicsPendulum.xml'
% analyzeTool = AnalyzeTool(filename) % construct an analyzeTool object from file
% pointKinematicsObject=analyzeTool ().getAnalysisSet().get('PointKinematics') % get a handle to PK object
% ....
% Note: You will need to edit all variables to be specific for your data and models. You may want to do this in the xml file before hand.
% 
% Once you have that all setup correctly you execute by using 'analyzeTool .run()' command
				
					%% MUSCLE ANALYSIS
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
                    %% STATIC OPTIMIZATION
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
               else
                    disp('Model file does not exist');
               end
          end
     end
end

%% PLOTTING ID AND EMG RESULTS

% Select folders to get to session folder
BasePath=uigetdir('..\..\..', 'Select Elaborated Data Folder');
fName = uigetdir(BasePath, 'Select the Subject for analysis');
SessionDirs = dir(fName);
isub=[SessionDirs(:).isdir];
sessionFolders={SessionDirs(isub).name}';
sessionFolders(ismember(sessionFolders,{'.','..', 'ROM'}))=[]; % dynamic subject folders

for folder = 1:length(sessionFolders)
     sessionName = [fName, filesep, sessionFolders{folder}];
     
     % Select the .osim model file for that condition
     [genericModelID,genericModelPath]=uigetfile([sessionName, filesep,...
          'staticElaborations', filesep, '*.osim'],...
          'Select scaled OpenSim model for analysis');
     model_file=[genericModelPath genericModelID];
     
     % Plot the joint moments
     plotMoments(sessionName, model_file);
     
     % Plot the EMG data and save as .mat files
     % 	loadPlotEMG(sessionName, fName);
     
end