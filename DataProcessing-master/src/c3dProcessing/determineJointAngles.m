function [anglesJoint] = determineJointAngles(c3dFiles, pname)
%% JOINT ANGLES DETERMINATION
%  Imports markers from a vicon nexus c3d file for load shraing ROM trials.
%  Use these markers to define the angle between specific segments
%  angles of interest -  hip trunk shoulder

%          - - - - - - - %% - - - - - - - -
%
% Based on code from Nathan at Auckland Bio-engineering Institute
%
% (computeCoordinateFramesAndEulerRotations)
%         - - - - - - - %% - - - - - - - -
%% ----- VICON Marker positions

% Loop through the ROM trials
for c3dAcq = 1:length(c3dFiles)
     
     % Acquire list of marker data
     [markersList] = getMarkersFromC3D(c3dFiles{c3dAcq}, pname);
     
     % Ground FRAME
     ground=[1 0 0 ; 0 1 0 ; 0 0 1];
     
     % Change the input angle cell array depending on c3d trial name
     romTrials = {'HF', 'UUA', 'ShoulderFF', 'TF'};
     
     % Inline function checking to see which ROM trial the c3dFile
     % corresponds with
     cellfind = @(string)(@(cell_contents)(strncmp(string, cell_contents, 2)));
     cell_array = romTrials;
     string = c3dFiles{c3dAcq};
     
     % Logical output with a 1 corresponding to trial of interest.
     logicalCells = cellfun(cellfind(string), cell_array);
     
     % logicalCells(1,1) = 1 is hip flexion trial
     % logicalCells(1,2) = 1 is upper arm abduction trial
     % logicalCells(1,3) = 1 is shoulder forward flexion trial
     % logicalCells(1,4) = 1 is trunk flexion trial
     
     if logicalCells(1,1) == 1
          
          % Cell array containing angles of interest - angles for hip flexion
          angles = {'PelTilt', 'PelList', 'PelRot', 'RHip_Flex'};
          
          % Determine angles
          [anglesList] = outputJointAngles(markersList, ground, angles);
          
     elseif logicalCells(1,2) == 1
          % Angles for Shoulder upper arm abduction
          angles = {'LShld_Flex', 'RShld_Flex', 'LShld_Add', 'RShld_Add'};

          % Determine angles
          [anglesList] = outputJointAngles(markersList, ground, angles);
          
     elseif logicalCells(1,3) == 1
          % Angles for shoulder FF
          angles = {'LShld_Flex', 'RShld_Flex', 'LShld_Add', 'RShld_Add'};
          
          % Determine angles
          [anglesList] = outputJointAngles(markersList, ground, angles);

     elseif logicalCells(1,4) == 1
          % Angles for trunk flexion
          angles = {'PelTilt', 'PelList', 'PelRot', 'LumbarExtension', 'LumbarBend',...
               'LumbarRotation'};
          
          % Determine angles
          [anglesList] = outputJointAngles(markersList, ground, angles);

     else
          disp([c3dFiles{c3dAcq}, ' is not a ROM trial, please select only ROM trials'])
     end
     
     % CALCULATE THE MAXIMUM, MINIMUM, AND RANGE OF EACH ANGLE
     % THESE VALUES ARE IN DEGREES
     
     % Empty structure for outputs
     anglesJoint = struct();
     outputs = {'Max', 'Min', 'Range'};
     
     for k = 1:length(anglesList)
          
          % Define angles for analysis
          anglesData = anglesList{k};
          
          % Calculate max, min, and ROM
          maximum = max(anglesData);
          minimum = min(anglesData);
          ROM = range(anglesData);
          
          % Save data to structure
          anglesJoint.(c3dFiles{c3dAcq}).(anglesList{k}).(outputs{1}) = maximum;
          anglesJoint.(c3dFiles{c3dAcq}).(anglesList{k}).(outputs{2}) = minimum;
          anglesJoint.(c3dFiles{c3dAcq}).(anglesList{k}).(outputs{3}) = ROM;
          
     end
end
