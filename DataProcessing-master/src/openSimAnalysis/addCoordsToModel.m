function addCoordsToModel(model_file)
%Assign knee abd/add and int/ext rotation DOFS to an input model
%   Input the directory to the model file.

%   Get the coordinates for both knee joints and add adduction and int/ext rotation to the coordinate
%   set. Then change the simmspline functions in the knee joint spatial transform 
%   to linear functions and set the coefficients to 1 and 0

% 11/10/2016 - Gavin Lenton

% Check argument inputs
narginchk(0, 1);

% If there aren't enough arguments passed in system will ask user to
% manually select file(s)
if nargin < 1
     [Model_In, modelpath] = uigetfile('.osim', 'Please select a model file');
     model_file = [modelpath, Model_In,];
end


%% Add coordinates and delete the simmsplines
import org.opensim.modeling.*

% Get model
current_model = Model(model_file);
% Create new model
new_model = Model(current_model);
new_model.setName('FBM_openKnee');
new_model.initSystem;

% Get the model coordinate set
modelJointSet = new_model.getJointSet();

sides = {'r', 'l'};

% Loop through both sides
for i = 1:2
     
     % Get knee joint
     kneeJointR = CustomJoint.safeDownCast(new_model.getJointSet().get(['knee_', sides{i}]));
     kneeCoordSet = kneeJointR.getCoordinateSet();
     
     if i == 1
          % Specify range for joints
          rotRangeAbd = [-0.349065, 0.349065];
          rotRangeIntRot = [-0.087 0.349065];
          
     else
          rotRangeAbd = [-0.349065, 0.349065];
          rotRangeIntRot = [-0.349065 0.087];
     end
     
     % Update coordinate sets
     flexCoord = kneeCoordSet.get(0);
     
     % Create clone of knee_angle_r for other DOFs
     adductionCoord = flexCoord.clone();
     intRotationCoord = flexCoord.clone();
     
     % Update name and range for other DOF
     adductionCoord.setName(['knee_adduction_', sides{i}]); adductionCoord.setRange(rotRangeAbd);
     intRotationCoord.setName(['knee_rotation_', sides{i}]); intRotationCoord.setRange(rotRangeIntRot);
     
     % Append coordinates to coordinate set
     kneeCoordSet.cloneAndAppend(adductionCoord);
     kneeCoordSet.cloneAndAppend(intRotationCoord);
     
     % Add coordinate set to joint
     kneeJointR.set_CoordinateSet(kneeCoordSet.clone());
     
     % Add joint to joint set
     modelJointSet.cloneAndAppend(kneeJointR);
     
     %% Update spatial transforms
     transformKnee = kneeJointR.getSpatialTransform();
     
     % Get transform for rotation axes 2 and 3
     knee_add = transformKnee.getTransformAxis(1);
     knee_rot = transformKnee.getTransformAxis(2);
     
     % Set coordinate names to match coordinates
     knee_add.set_coordinates(0, ['knee_adduction_', sides{i}]);
     knee_rot.set_coordinates(0, ['knee_rotation_', sides{i}]);
     
     % Create the linear function with coefficients [1,0]
     linFunction = LinearFunction();
     functionCoeefficients = ArrayDouble();
     functionCoeefficients.setValues([1,0], 2)
     linFunction.setCoefficients(functionCoeefficients);
     
     % Set function to the transform
     knee_add.setFunction(linFunction);
     knee_rot.setFunction(linFunction);
     
     % Add transforms back to model for rotations 2 and 3 (non flexion
     % rotations)
     transformKnee.set_rotation2(knee_add.clone());
     transformKnee.set_rotation3(knee_rot.clone());
     
end

% Print the new model
new_model.print(regexprep(model_file, 'FBM','openKnee'));

