function data = LS_pipeline_treadmill(file)

% This function will create the appropriate files to run an OpenSim
% simulation sequence for the stumble recovery data
% Input - file - c3d file to process
% Originally written by Glen Lichtwark (The University of Queensland)
% Copyright 2011
% Contributions by David Graham and Chris Carty (Griffith University) 2012
% 
% This version is for use in Balance Recovery trials where data was collected
% on the instrumented treadmil. Many alterations to the initial code have
% been implemented such as migration to BTK server functions and migration to
% OpenSim 3.2 for API control functionality. 
%
% This version utlises static optimisation to estimate muscle force production
% as such the satatic optimisation outputs are required as inputs to IAA.
% This implements IAA as a plug-in which can be found on Tim Dorns' OpenSim
% page at simTK.org. Please cite Tim if you use this analysis.


%% Define path and model

global subdir_current;
global drive;



osim_path = [drive '\Opensim_Models\FullBodyModel\'];

model = 'FullBodyModel';

if nargin < 1
    [fname, pname] = uigetfile('*.c3d', 'Select C3D file');
else 
if isempty(fileparts(file))
    pname = cd;
    pname = [pname '\'];
    fname = file;
else [pname, name, ext] = fileparts(file);
    fname = [name ext];
end
end

cd(pname);

%% Load the c3dfile

% Load the c3d file using BTK
% Threshold set to 100N because of noise on force plates
[data.marker_data, data.analog_data, data.fp_data, data.sub_info] = btk_loadc3d([pname, fname], 100); 

% BTK uses First/Last so convert to fit existing routine so Start/End
data.marker_data.Last_Frame = data.marker_data.Last_Frame - data.marker_data.First_Frame;
data.marker_data.First_Frame = 1;

data.Start_Frame = data.marker_data.First_Frame;
data.End_Frame = data.marker_data.Last_Frame;

%Create a log folder to place log files for each tool
Log_folder = [cd '\Logs'];
newfile = mkdir([Log_folder]);

if isempty(strfind(lower(fname),'cal'))
%% C3D TO TRC
%'v5' of this script is modified from previous ones since person is now
%walking in the 'forward' direction of the treadmill (towards the wall) and
%so coordinate system transformations have changed.

data = btk_c3d2trc_v5(data,'off'); % animation off

% Add the mass of the harness
%data.Mass = data.Mass + 3.5;

%% Define file name prefix structure
% Added a nameing structure to avoid confusion throughout the script,
% participant code refers to the actual participant number (eg '066')
% where as trial code refers to the trial (eg '20_4')
    
participant_code    = data.Name;

trial_code          = fname(1:end-4);

data.c3d_filename   = trial_code;

%point_kinematics;

%% INVERSE KINEMATIC ANALYSIS

% DEFINE STANDARD OPENSIM FILES AND PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ModelFile           = [participant_code '_SCALED.osim'];
ResultsDirectory    = [pname 'IK_Results\'];
OutputFile          = [trial_code '_IK.sto'];
IKTasksFile         = [osim_path 'FullBodyModel_IK_Tasks.xml'];

    % Set up the XML file   
    setup_InverseKinematics('data',data,...
                            'ModelFile',ModelFile,...
                            'IKTasksFile',IKTasksFile,...
                            'OutputFile', OutputFile,...
                            'ResultsDirectory',ResultsDirectory,...
                            'Accuracy',0.000005);

% CALL THE IK TOOL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,log_mes]=dos(['ik -S ', data.Name '_Setup_InverseKinematics.xml'],'-echo');

% SAVE THE WORKSPACE AND PRINT A LOG FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen([Log_folder,'\IK_',data.Name,'.log'],'w+');
fprintf(fid,'%s\n', log_mes);
fclose(fid);

save([(subdir_current) '.mat']');


%% INVERSE DYNAMIC ANALYSIS

% DEFINE STANDARD OPENSIM FILES AND PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
MOTFile             = [pname trial_code '_IK.sto'];
GRFFile             = [pname trial_code '_grf.mot'];
GRFFileXML          = [pname trial_code '_grf.xml'];
OutputFile          = [trial_code '_ID.sto'];
ResultsDirectory    = [pname 'ID\'];

    % Write GRF XML file    
    grf2xml(data,...
            'ExternalLoadNames',{'ExternalForce_1','ExternalForce_2'},...
            'AppliedToBodies',{'calcn_r','calcn_l'},...
            'GRFFile',GRFFile,...
            'MOTFile',MOTFile,...
            'LowPassFilterForKinematics',6,...
            'OutputFile',GRFFileXML);

    % Set up the XML file 
    setup_InverseDynamics('data',data,...
                          'ModelFile',ModelFile,...
                          'MOTFile',MOTFile,...
                          'GRFFile',GRFFileXML,...
                          'LowPassFilterFreq',6,...
                          'OutputFile', OutputFile,...
                          'ResultsDirectory', ResultsDirectory); 


% CALL THE ID TOOL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,log_mes]=dos(['id -S ', participant_code '_Setup_InverseDynamics.xml'],'-echo');

% SAVE THE WORKSPACE AND PRINT A LOG FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen([Log_folder,'\ID_',data.Name,'.log'],'w+');
fprintf(fid,'%s\n', log_mes);
fclose(fid);

clear D

% Load the data from ID MOT file and tie kinetic data to data structure for
% later use in plotting for residual check at RRA
D = load_sto_file([pname 'ID\' trial_code '_ID.sto']);
data.Kinetics = D;

save([(subdir_current) '.mat']');

else % Scaling process for the static trial

 data = btk_c3d2trc_scale_v5(data, 'off');
    
% DEFINE STANDARD OPENSIM FILES AND PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ModelFile              = [osim_path model '.osim'];
 MeasurementSetFile     = [osim_path model '_Scale_MeasurementSet_PBT.xml'];
 ScaleTasksFile         = [osim_path model '_Scale_Tasks_PBT.xml'];

% Set up the XML file 
 setup_scale('data',data,...
             'ModelFile',ModelFile,...
             'ScaleTasksFile',ScaleTasksFile,...
             'MeasurementSetFile',MeasurementSetFile,...
             'PreserveMass','true');


% CALL THE SCALE TOOL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,log_mes]=dos(['scale -S ' data.Name '_Setup_Scale.xml'],'-echo');

% SAVE THE WORKSPACE AND PRINT A LOG FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen([Log_folder,'\SCALE_',data.Name,'.log'],'w+');
fprintf(fid,'%s\n', log_mes);
fclose(fid);

save([(subdir_current) '.mat']');

    
end



