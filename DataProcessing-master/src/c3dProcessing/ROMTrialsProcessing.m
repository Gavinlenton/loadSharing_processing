function  ROMTrialsProcessing(pname, sessionConditions, fName)
%Evaluates the range of motion trials and saves output angles to an xml
%file
%   Input the name of directory containing the ROM c3d files and evaluate
%   the files to generate maximum, minimum, and range of joint angles for
%   each trial.

     % Re-set folder to include the ROM trials only
     c3dFile_ROM = pname;
     c3dFilesROM=dir([c3dFile_ROM, filesep, '*.c3d']);
     
     % Delete files I don't want to analyse
     c3dFilesForROM = {c3dFilesROM.name}';
     c3dFilesForROM = selectROMTrials(c3dFilesForROM, sessionConditions);
    
     % Define mat file name
     matFileDir = [regexprep(fName, 'Input', 'Elaborated'), filesep, 'ROM'];
     fileName = 'romData.mat';
     
     if ~isdir(matFileDir)
          mkdir(matFileDir)
     end
     
     % Output max, min, and range of joint angles
     [anglesJoints] = determineJointAngles(c3dFilesForROM, pname);
     
     if exist([matFileDir, filesep, fileName], 'file');
         load([matFileDir, filesep, fileName]);
     
     % Function to average the three trials
     [anglesJointMean] = findMeanOfROMTrials(anglesJoints, sessionConditions, matFileDir);
    
          % Save each session in separate tabs
          for tt = 1:length(sessionConditions)
               session = sessionConditions{tt};
			   % Only if the condition had ROM trials
			   if isfield(anglesJointMeans, session)
				   anglesJointMeans.(session) = anglesJointMean.(session);
			   end
          end

          save([matFileDir, filesep, fileName], 'anglesJointMeans');
     else
         
          [anglesJointMeans] = findMeanOfROMTrials(anglesJoints, sessionConditions, matFileDir);
          save([matFileDir, filesep, fileName], 'anglesJointMeans');
     end
     clearvars anglesJointMean anglesJoint matFileDir
end

