function addCoordsToModel(model_file)
%Assign knee abd/add and int/ext rotation DOFS to an input model
%   Input the model file, and directory to the session

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
state = current_model.initSystem();

% updCoord = current_model.updCoordinateSet();

% Create new coord set
newCoordSet = CoordinateSet();
axisValues = Vec3();

% Get the model coordinate set
modelCoordSet = current_model.getCoordinateSet();
kneeJoint = CustomJoint.safeDownCast(current_model.getJointSet().get('knee_r'));
kneeCoordinateSet = kneeJoint.getCoordinateSet();

kneeJoint.append_CoordinateSet(coord2);

kneeCoordinateSet.get(0).setRange(rotRangeAbd);
kneeCoordinateSet.set(1).setName;
kneeCoordinateSet.set(2).setRange(rotRangeAbd);


kneeCoord = modelCoordSet.get('knee_angle_r');
joint = kneeCoord.getJoint;
eval(['concreteJoint = ' char(joint.getConcreteClassName) '.safeDownCast(joint);'])
sptr = concreteJoint.getSpatialTransform; 

for ip = 0 : 5
     if strcmp(char(sptr.getCoordinateNames().get(ip)), char(kneeCoord.getName))
          sptr.getTransformAxis(ip).getAxis(axisValues);
          break
     end
end


% Get knee joint
kneeJoint = CustomJoint.safeDownCast(current_model.getJointSet().get('knee_r'));

% Get spatial transform
transformKnee = kneeJoint.getSpatialTransform();

% Get transform for rotation axes 2 and 3
A = transformKnee.getTransformAxis(1);
B = transformKnee.getTransformAxis(2);

% Set coordinate names to match coordinates
A.set_coordinates(0, 'knee_adduction_r');
B.set_coordinates(0, 'knee_rotation_r');

% Create the linear function with coefficients [1,0]
linFunction = LinearFunction();
functionCoeefficients = ArrayDouble();
functionCoeefficients.setValues([1,0], 2)
linFunction.setCoefficients(functionCoeefficients);

% Set function to the transform
A.setFunction(linFunction);
B.setFunction(linFunction);

% Specify range for joints
rotRangeAbd = [-1.0943951, 1.0943951];
rotRangeIntRot = [-0.087 0.349065];

% Update coordinate sets
kneeJoint = CustomJoint.safeDownCast(current_model.getJointSet().get('knee_r'));
kneeCoordSet = kneeJoint.getCoordinateSet();
flexCoord = kneeCoordSet.get(0);
% Create clone of knee_angle_r for other DOFs
adductionCoord = flexCoord.clone();
intRotationCoord = flexCoord.clone();

% Update name and range for other DOF
adductionCoord.setName('knee_adduction_r'); adductionCoord.setRange(rotRangeAbd);
intRotationCoord.setName('knee_rotation_r'); intRotationCoord.setRange(rotRangeIntRot);

kneeJoint.set_CoordinateSet(adductionCoord)
kneeJoint.upd_CoordinateSet(adductionCoord);
kneeJoint.upd_CoordinateSet().set(adductionCoord);

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
