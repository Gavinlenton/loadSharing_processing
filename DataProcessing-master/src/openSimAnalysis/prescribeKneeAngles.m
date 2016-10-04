function prescribeKneeAngles(IKoutputDir, model_file)
%Assign knee abd/add and int/ext rotation values as a function of knee
%flexion angles using simmspline values from a musculoskeletal model
%   Input the model file and directory to IK files

import org.opensim.modeling.* 

model = Model(model_file);

% Load the model and initialize
state = model.initSystem();

mname = model.getName();

model.setName(mname);

functions = model.getFunctionClassNames();

value = functions.equals();

joints = model.getJointSet();

joints.something

% Generate list of trials
trials=dir(IKoutputDir);
j=1;

for k = 3:length(trials)
     trialsList{j}=trials(k).name;
     j = j + 1;
end
trialsList(ismember(trialsList,{'Figures','IDMetrics.mat', 'out.log', 'error.log', 'AnalysedData'}))=[];

end

