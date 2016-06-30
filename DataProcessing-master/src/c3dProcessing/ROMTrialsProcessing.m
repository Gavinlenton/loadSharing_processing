function  ROMTrialsProcessing(pname, subjectNumber, sessionConditions)
%Evaluates the range of motion trials and saves output angles to an xml
%file
%   Input the name of directory containing the ROM c3d files and evaluate
%   the files to generate maximum, minimum, and range of joint angles for
%   each trial.

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
     xmlwrite(xmlFileName, anglesJoint, ['Session', num2str(subjectNumber)]);

end

