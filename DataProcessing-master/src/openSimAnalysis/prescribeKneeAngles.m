function prescribeKneeAngles(model_file, sessionName, IKoutputDir)
%Assign knee abd/add and int/ext rotation values as a function of knee
%flexion angles using simmspline values from a musculoskeletal model
%   Input the model file, and directory to the session and IK files

% 07/10/2016 - Gavin Lenton

% Check argument inputs
narginchk(0, 3);

% If there aren't enough arguments passed in system will ask user to
% manually select file(s)
if nargin < 1
     [Model_In, modelpath] = uigetfile('.osim', 'Please select a model file');
     sessionName = uigetdir(cd,  'Select folder with INVERSE KINEMATICS results to use');
     IKoutputDir = uigetdir(sessionName, 'Select folder with INVERSE KINEMATICS results to use');
     model_file = [modelpath, Model_In,];
elseif nargin < 2
     sessionName = uigetdir(cd,  'Select folder with INVERSE KINEMATICS results to use');
     IKoutputDir = uigetdir(sessionName, 'Select folder with INVERSE KINEMATICS results to use');
elseif nargin < 3 % If no IK output or session folder name
     IKoutputDir = uigetdir(sessionName, 'Select folder with INVERSE KINEMATICS results to use');
end

% Generate list of trials
trials=dir(IKoutputDir);
j=1;

for k = 3:length(trials)
     trialsList{j}=trials(k).name;
     j = j + 1;
end
trialsList(ismember(trialsList,{'Figures','IDMetrics.mat', 'out.log', 'error.log', 'AnalysedData'}))=[];

% Be selective because I want to exclude bad trials
[trialsIndex,~] = listdlg('PromptString','Select trials to plot:',...
     'SelectionMode','multiple',...
     'ListString',trialsList);

inputTrials=trialsList(trialsIndex);

%% Get the simmspline values and use them to compute knee adduction and rotation value
import org.opensim.modeling.*

% Get model
current_model = Model(model_file);
kneeJoint = CustomJoint.safeDownCast(current_model.getJointSet().get('knee_r'));
transformKnee = kneeJoint.getSpatialTransform();

% Get SimmSpline values
% Abd/add
jointFuncR2 = transformKnee.get_rotation2().getFunction();
simmsplineR2 = SimmSpline.safeDownCast(jointFuncR2);

% Int/ext rotation
jointFuncR3 = transformKnee.get_rotation3().getFunction();
simmsplineR3 = SimmSpline.safeDownCast(jointFuncR3);

% Get X and Y values - only need one X value because they are both knee
% flexion
simmsplineX2 = simmsplineR2.getX;
simmsplineY2 = simmsplineR2.getY;
simmsplineRY3 = simmsplineR3.getY;

% Determine number of values in spline
splineLength = simmsplineX2.getSize();
% Get the values from the Vec3
splineValuesAbd = zeros(1, splineLength);
splineValuesRot = zeros(1, splineLength);
splineValuesFlex = zeros(1, splineLength);

% Manually defined the spline values because I couldn't get length
% -1 when accessing OpenSim values because it starts at 0
for splineValue = 1:splineLength
     splineValuesFlex(1, splineValue) = simmsplineX2.get(splineValue-1);
     splineValuesAbd(1, splineValue) = simmsplineY2.get(splineValue-1);
     splineValuesRot(1, splineValue) = simmsplineRY3.get(splineValue-1);
end

% Rotation for left leg is opposite sign.
splineValuesRot2 = splineValuesRot .* -1;

% Initialize model
osimModel= current_model;

% Get the model coordinate name from the new model for knee ab/add and
% int/ext rot
modelCoordSet = osimModel.getCoordinateSet();

%% AUTOMATE THIS
% Get the knee_angle_r and knee_angle_l coordinate set from the model
currentcoord = modelCoordSet.get('knee_angle_r');
currentcoord2 = modelCoordSet.get('knee_angle_l');

% Loop through the trials
for j = 1:length(inputTrials)
     
     % Define path to mot file and obtain data values
     motfilepath = [IKoutputDir, filesep, inputTrials{j}, filesep, 'FBM_ik.mot'];
     
     if exist(motfilepath, 'file')
          data = load_sto_file(motfilepath);
          
          % Create the coordinate storage object from the input .sto file
          coordinateSto=Storage(motfilepath);
          
          % construct ArrayDouble objects for time and coordinate values
          Time = ArrayDouble();
          coordvalue = ArrayDouble();
          coordvalue2 = ArrayDouble();
          
          % Get the Time stamps and Coordinate values
          coordinateSto.getTimeColumn(Time);
          coordinateSto.getDataColumn(currentcoord.getName(),coordvalue);
          coordinateSto.getDataColumn(currentcoord2.getName(),coordvalue2);
          
          % Initiate knee flexion angles
          kneeFlexR = zeros(1,coordvalue.getSize());
          kneeFlexL = zeros(1,coordvalue.getSize());
          
          % Loop through each point in the trial and get knee flexion value
          for k = 0:coordvalue.getSize()-1
               kneeFlexR(k+1) = coordvalue.getitem(k)/(180/pi);
               kneeFlexL(k+1) = coordvalue2.getitem(k)/(180/pi);
          end
          
          % Run polyfit and polyval to determine the other knee DOF values from
          % the knee flexion angle and the spline. Fits the 3rd degree (cubic)
          % polynomial to the points as per the simmspline
          [x1, ~, mu1] = polyfit(splineValuesFlex, splineValuesAbd, 3);
          [x2, ~, mu2] = polyfit(splineValuesFlex, splineValuesRot, 3);
          [x3, ~, mu3] = polyfit(splineValuesFlex, splineValuesRot2, 3);
          
          y2 = rad2deg(polyval(x1, kneeFlexR, [], mu1));
          y2L = rad2deg(polyval(x1, kneeFlexL, [], mu1));
          y3 = rad2deg(polyval(x2, kneeFlexR, [], mu2));
          y3L = rad2deg(polyval(x3, kneeFlexL, [], mu3));
          
          % Plot to check polynomial results
          %     plot(splineValuesFlex, splineValuesRot, 'o');
          %     hold on
          %     plot(kneeFlex, y3);
          
          % Add the new data to columns in .mot file  
          data.('knee_adduction_r') = y2';
          data.('knee_rotation_r') = y3';
          data.('knee_adduction_l') = y2L';
          data.('knee_rotation_l') = y3L';
          
          % Reverse polarity of knee angles because the fullbodymodel uses
          % left hand conventino instead of right hand
%           data.knee_angle_r = data.knee_angle_r.*-1;
%           data.knee_angle_l = data.knee_angle_l.*-1;
          
          %  Save the Modified data to a file
          fileoutpath = regexprep(motfilepath, 'FBM', 'openKnee');
          write_sto_file(data,fileoutpath);
          
     else
          disp(['trial ' inputTrials{j} 'cannot be found at ' IKoutputDir]);
     end
end
     
     

