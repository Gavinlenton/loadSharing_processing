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
     
     % Function to average the three trials
     [anglesJointMeans] = findMeanOfROMTrials(anglesJoints, sessionConditions, matFileDir);
     
     if exist([matFileDir, filesep, fileName], 'file');
         ROM =  load([matFileDir, filesep, fileName]);
          
          % Save each session in separate tabs
          for tt = 1:length(sessionConditions)
               session = sessionConditions{tt};
               ROM.anglesJointMeans.(session) = anglesJointMeans.(session);
          end
          
          save([matFileDir, filesep, fileName], 'anglesJointMean', 'anglesJoints');
     else
          
          save([matFileDir, filesep, fileName], 'anglesJointMeans', 'anglesJoints');
     end
     clearvars anglesJointMean anglesJoint matFileDir
end

