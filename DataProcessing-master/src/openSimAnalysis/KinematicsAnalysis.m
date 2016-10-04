function [KAoutputDir,KAtrialsOutputDir,inputTrials, varargout] = KinematicsAnalysis(inputDir, model_file, IKoutputDir, trialsList)
% Set and run kinematics analysis for multiple trials

%% Setting KA

switch nargin
    
    case 2  %no IK before
        
        trialsList = trialsListGeneration(inputDir);
        [KAid, KATemplateXml, fcut_coordinates, inputTrials, IKmotDir] = KAinput(trialsList);
        
    case 3 %IK before, but KA on different trials (for which it has been performed IK)
        
        trialsList = trialsListGeneration(IKoutputDir);
        [KAid, KATemplateXml, fcut_coordinates, inputTrials] = KAinput(trialsList);
        IKmotDir=IKoutputDir;
        
    case 4 %IK before and ID on the same trials
        
        [KAid, KATemplateXml, fcut_coordinates] = MAinput();
        IKmotDir=IKoutputDir;
        inputTrials=trialsList;     
end


%% Running KA

[KAoutputDir, KAtrialsOutputDir]=runMuscleAnalysis(inputDir,inputTrials, model_file, IKmotDir, KAid, KATemplateXml, fcut_coordinates);