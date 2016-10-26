%% Script to run IK, ID, and any further analysis on load sharing data %%

% IK, ID, muscle analysis, and SO taken from Batch OpenSim Processing Scripts (BOPS)
% Copyright (C) 2015 Alice Mantoan, Monica Reggiani
% <https://simtk.org/home/bops>
% Please acknowledge the authors of BOPS when using this script.
clc; clear;

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
               modelName = [regexprep(subjectNames{nS}, ' ', ''), '_FBM_optmkrpos.osim'];
               
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
                    
                    % Make sure the hip DOFs are unlocked and the subtalar DOF is locked
                    joints = {'hip', 'subtalar'};
                    islocked = [0 1];
                    lockDOFS([sessionName, filesep, 'staticElaborations', filesep,...
                         ['Static1_' trialName], filesep, 'StaticCal', filesep], modelName, joints, islocked);
                    
                    % Navigate to folder where BOPS is;
                    cd([folderBOPS, filesep, 'src']);
                    
                    %% INVERSE KINEMATICS and INVERSE DYNAMICS
                    
                    procIdIndex = regexp(trialsFolders{tD}, '_P');
                    IK_id = trialsFolders{tD}(1:procIdIndex-1);
                    
                    % Check to see if IK has already been run
                    if ~exist([sessionName, filesep, 'inverseKinematics', filesep, IK_id, filesep, [IK_id, '_Processed1'], filesep, 'FBM_ik.mot'], 'file')
                         
                         % Run IK 
                         [IKoutputDir, IKtrialsOutputDir, IKprocessedTrials]=InverseKinematics(inputDir, model_file);
                         
                         % Plot IK results
                         IKmotFilename='FBM_ik.mot'; %our default name
                         [coordinates,Xaxislabel]=plotResults('IK', IKoutputDir, IKtrialsOutputDir, model_file, IKprocessedTrials, IKmotFilename);
                         
                         % Delete crappy trials if they exist
                         uiwait
                         
                         % Then prescribe motion to other DOFS at the knee
                         % in new model
                         if ~exist([sessionName, filesep, 'inverseKinematics', filesep, IK_id, filesep, [IK_id, '_Processed1'], filesep, 'openKnee_ik.mot'], 'file')
                              
                              % Update to the knee unlocked model file
                              model_file_kneeUnlocked = regexprep(model_file, 'FBM_optmkrpos', 'openKnee_optmkrpos');
                              
                              % Determine other knee DOF angles using
                              % simmspline values
                              prescribeKneeAngles(model_file, sessionName, IKoutputDir);
                              
                              % Check to see if ID has been run
                              if ~exist([sessionName, filesep, 'inverseDynamics', filesep, IK_id, filesep, [IK_id, '_Processed1'], filesep, 'inverse_dynamics.sto'], 'file')
                                   
                                   % Run ID
                                   [IDoutputDir, IDtrialsOutputDir, IDprocessedTrials]=InverseDynamics(inputDir, model_file_kneeUnlocked, IKoutputDir);
                                   
                                   % Plot ID results
                                   IDmotFilename='inverse_dynamics.sto'; %our default name
                                   
                                   [coordinates,Xaxislabel]=plotResults('ID', IDoutputDir, IDtrialsOutputDir, model_file_kneeUnlocked, IDprocessedTrials, IDmotFilename);
                                        
                                   % Pause UI and delete the bad trials
                                   uiwait
                              end
                         end
                         
                         % Check to see if it has been run with older model
                    elseif ~exist([sessionName, filesep, 'inverseKinematics', filesep, IK_id, filesep, [IK_id, '_Processed1'], filesep, 'openKnee_ik.mot'], 'file')
                         
                         % Update to the knee unlocked model file
                         model_file_kneeUnlocked = regexprep(model_file, 'FBM_optmkrpos', 'openKnee_optmkrpos');
                         
                         % Determine other knee DOF angles using
                         % simmspline values
                         prescribeKneeAngles(model_file, sessionName);
                           
                         % Check to see if ID has been run
                         if ~exist([sessionName, filesep, 'inverseDynamics', filesep, IK_id, filesep, [IK_id, '_Processed1'], filesep, 'inverse_dynamics.sto'], 'file')
                              
                              % Run ID
                              [IDoutputDir, IDtrialsOutputDir, IDprocessedTrials]=InverseDynamics(inputDir, model_file_kneeUnlocked, IKoutputDir);

                              % Plot ID results
                               IDmotFilename='inverse_dynamics.sto'; %our default name

                               [coordinates,Xaxislabel]=plotResults('ID', IDoutputDir, IDtrialsOutputDir, model_file, IDprocessedTrials, IDmotFilename);
                              % Pause UI and delete the bad trials
                              uiwait
                         end
                         
                         % Check to see if ID has been run
                    elseif ~exist([sessionName, filesep, 'inverseDynamics', filesep, IK_id, filesep, [IK_id, '_Processed1'], filesep, 'inverse_dynamics.sto'], 'file')
                         
                         % Update to the knee unlocked model file
                         model_file_kneeUnlocked = regexprep(model_file, 'FBM_Scaled_Pose_OptMusc', 'openKnee_optmkrpos');
                         
                         % Run ID on different trials as IK
                         [IDoutputDir, IDtrialsOutputDir, IDprocessedTrials]=InverseDynamics(inputDir, model_file_kneeUnlocked, IKoutputDir);
                         
                         % Plot ID results
                         IDmotFilename='inverse_dynamics.sto'; %our default name

                         [coordinates,Xaxislabel]=plotResults('ID', IDoutputDir, IDtrialsOutputDir, model_file, IDprocessedTrials, IDmotFilename);
    
                         % Pause UI and delete the bad trials
                         uiwait
                         
                    else
                         disp('Already processed IK and ID for this condition')
                         
                    end
                    
                    
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

%% Data analysis for ID, IK, and raw EMG results

% Select folders to get to session folder
BasePath=uigetdir('..\..\..', 'Select Elaborated Data Folder');
fName = uigetdir(BasePath, 'Select the Subject for analysis');
SessionDirs = dir(fName);
isub=[SessionDirs(:).isdir];
sessionFolders={SessionDirs(isub).name}';
sessionFolders(ismember(sessionFolders,{'.','..', 'ROM', 'AnalysedData', 'Figures'}))=[]; % dynamic subject folders

for folder = 1:length(sessionFolders)
     sessionName = [fName, filesep, sessionFolders{folder}];
     
     % Select the .osim model file for that condition
     [genericModelID,genericModelPath]=uigetfile([sessionName, filesep,...
          'staticElaborations', filesep, '*.osim'],...
          'Select scaled OpenSim model for analysis');
     model_file_kneeUnlocked=[genericModelPath genericModelID];
     
     % Plot the kinematics result
     plotKinematics(sessionName, model_file_kneeUnlocked);
     
     % Plot the joint moments
     plotMoments(sessionName, model_file_kneeUnlocked);
     
     % Plot the EMG data and save as .mat files
     % 	loadPlotEMG(sessionName, fName);
end