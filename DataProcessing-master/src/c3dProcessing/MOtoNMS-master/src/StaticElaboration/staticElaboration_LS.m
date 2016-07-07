function  staticElabPath = staticElaboration_LS(subjectNames, motoDir, BasePath)
%Runs the MOtoNMS staticElaboration for the load sharing data
%   Loops through all subjects and sessions within the subject to process
%   the static trials for further analysis. This includes processing joint
%   centres, and generating a static.trc file for scaling and mkrplacer functions
%   to be run on an OpenSim model.

%%% INPUTS %%%

% 1: Cell array containing all subjects in analysis
% 2: Path to where MOtoNMS is stored on your computer
% 3: Path to the ElaboratedData folder on your computer

%%% OUTPUTS %%%

% 1: Path to the staticElaboration folder for each subject and session

%% Structure for path names creation

for i = 1:length(subjectNames)
     subjectName = subjectNames{i}(~isspace(subjectNames{i}))
     staticElabPath.(subjectName) = struct()
end

for nS = length(subjectNames)
     
     % Remove space from subject name
     subjectName = subjectNames{nS}(~isspace(subjectNames{nS}));
     
     % Create cell array containing session folders for chosen subject
     SessionDirs = dir([BasePath, filesep, subjectNames{nS}]);
     isub=[SessionDirs(:).isdir];
     sessionFolders={SessionDirs(isub).name}';
     sessionFolders(ismember(sessionFolders,{'.','..', 'ROM'}))=[]; % dynamic subject folders
     
     % Need to use session number because we cannot use numerics as a valid
     % field name in staticElabPath structure below.
     sessionNumber = {'Session1', 'Session2', 'Session3', 'Session4'};
     
     
     for nSess = 1:length(sessionFolders) % CHECK TO MAKE SURE THAT FIRST SESSION ALWAYS HAS THE Static1_Processed trial
          
          % Define staticElabPath and folder with static c3d files
          staticElabPath.(subjectName).(sessionNumber{nSess}) = [BasePath,...
               filesep, subjectNames{nS}, filesep,  sessionFolders{nSess}, filesep, 'staticElaborations'];
          c3dFolderDirs = dir(regexprep(staticElabPath.(subjectName).(sessionNumber{nSess}), 'staticElaborations', 'sessionData'));
          isubby=[c3dFolderDirs(:).isdir];
          c3dFolder={c3dFolderDirs(isubby).name}';
          c3dFolder(~strncmp(c3dFolder, 'Static1', 7)) = [];
          c3dFolder(strncmp(c3dFolder, 'Static1_Processed', 10)) = [];
          
          % If not first session then I won't run the default static
          % Elaboration, instead skip to modified down below.
          if nSess == 1
               
               % Check to see if staticElaborations has already been performed
               if exist(staticElabPath.(subjectName).(sessionNumber{nSess}), 'dir') == 7
                    % If so print to screen and make dir of path
                    fprintf('\nStatic elaboration exists for %s in %s\n', subjectNames{nS}, sessionNumber{nSess});
                    staticFiledir=dir([staticElabPath.(subjectName).(sessionNumber{nSess}) filesep ]);
                    staticFileName='Static1_Processed';
                    staticFileFullPath=[staticElabPath.(subjectName).(sessionNumber{nSess}) filesep staticFileName filesep 'StaticCal'];
                    
                    % Run staticElab on remaining conditions in session

                    for nc3d = 1:length(c3dFolder)
                         
                         % For remaining static trials use the first static.xml to
                         % run elaborations.
                         
                         if exist([staticElabPath.(subjectName).(sessionNumber{nSess}) filesep c3dFolder{nc3d}], 'dir') == 7
                              fprintf('\nStatic elaboration exists already for %s in %s for condition %s\n',...
                                   subjectNames{nS}, sessionNumber{nSess}, c3dFolder{nc3d});
                         else
                              [foldersPaths,parameters] = StaticElaborationSettingsLS(staticFileFullPath,...
                                   sessionFolders{nSess}, c3dFolder{nc3d});
                              
                              % Run static elaboration with updated info
                              cd([motoDir, filesep, 'src', filesep, 'StaticElaboration']);
                              runStaticElaborationLS(staticFileFullPath, foldersPaths,...
                                   parameters);
                         end
                    end
               else
                    
                    % Otherwise navigate to staticElaboration src in MOtoNMS and
                    % run staticElaboration main function
                    cd([motoDir, filesep, 'src', filesep, 'StaticElaboration']);
                    run StaticInterface.m;
                    staticFiledir=dir([staticElabPath.(subjectName).(sessionNumber{nSess}) filesep ]);
                    staticFileName='Static1_Processed';
                    staticFileFullPath=[staticElabPath.(subjectName).(sessionNumber{nSess}) filesep staticFileName filesep 'StaticCal'];
                    
                    % Then run staticElab on remaining conditions in
                    % session
                    for nc3d = 1:length(c3dFolder)
                         
                         % For remaining static trials use the first static.xml to
                         % run elaborations.
                         
                         if exist([staticElabPath.(subjectName).(sessionNumber{nSess}) filesep c3dFolder{nc3d}], 'dir') == 7
                              fprintf('\nStatic elaboration exists already for %s in %s for condition %s\n',...
                                   subjectNames{nS}, sessionNumber{nSess}, c3dFolder{nc3d});
                         else
                              [foldersPaths,parameters] = StaticElaborationSettingsLS(staticFileFullPath,...
                                   sessionFolders{nSess}, c3dFolder{nc3d});
                              
                              % Run static elaboration with updated info
                              cd([motoDir, filesep, 'src', filesep, 'StaticElaboration']);
                              runStaticElaborationLS(staticFileFullPath, foldersPaths,...
                                   parameters);
                         end
                    end
                    
               end
               
               % Skip to modified elaboration for sessions 2+
          else
               
               % If other static trials for armour systems haven't been
               % processed then loop through them here and generate static.trc
               % files
               for nc3d = 1:length(c3dFolder)
                    
                    % For remaining static trials use the first static.xml to
                    % run elaborations.
                    
                    if exist([staticElabPath.(subjectName).(sessionNumber{nSess}) filesep c3dFolder{nc3d}], 'dir') == 7
                         fprintf('\nStatic elaboration exists already for %s in %s for condition %s\n',...
                              subjectNames{nS}, sessionNumber{nSess}, c3dFolder{nc3d});
                    else
                         [foldersPaths,parameters] = StaticElaborationSettingsLS(staticFileFullPath,...
                              sessionFolders{nSess}, c3dFolder{nc3d});
                         
                         % Run static elaboration with updated info
                         cd([motoDir, filesep, 'src', filesep, 'StaticElaboration']);
                         runStaticElaborationLS(staticFileFullPath, foldersPaths,...
                              parameters);
                    end
               end
               close all
          end
     end
end

end

