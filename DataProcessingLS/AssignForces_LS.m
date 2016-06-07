%% Load walking trial and add events based on FP threshold and marker location

% Load the trial and assign data to variable 'data'.
% select the c3d file
[fname, pname] = uigetfile('*.c3d', 'Select C3D file');

%Open BTK acquisition to get access to orig force data

% load the c3d file 
data = btk_loadc3d([pname, fname], 10);

% input the details of the participant and save to data structure  (usually
% this can be passed from another function)
data.Name = 'Name';
data.Mass = '75';

% % % marker list for the static trial
% marker_list = {'LASI'; 'RASI'; 'LPSI'; 'RPSI'; ...
%     'RTHI1'; 'RTHI2'; 'RTHI3'; 'RLKN'; ...
%     'RSHA1'; 'RSHA2'; 'RSHA3'; 'RLMAL'; ...
%     'RHEEL'; 'RMET5'; 'RMET1'; ...
%     'LTHI1'; 'LTHI2'; 'LTHI3'; 'LLKN'; ...
%     'LSHA1'; 'LSHA2'; 'LSHA3'; 'LLMAL'; ...
%     'LHEEL'; 'LMET5'; 'LMET1'};

% Assign right and left foot markers
markers = {'RCAL';'LCAL'};

% Run function to assign events
[data] = grf_assign_treadmill_Lenton(data,markers);

% define the filter frequency that will be used throughout
data.FilterFreq = 30;

% define start and end frames (modify appropriately - assume all to start)
data.Start_Frame = data.marker_data.First_Frame;
data.End_Frame = data.marker_data.Last_Frame;

% assign the forces - note that an individual force vector is assigned for
% each plate and there can therefore be two forces (fore/aft plates)
% assigned to the same segment at the same time). Therefore forces are not
% 'combined' they are treated separately, but will still get the same
% inverse dynamics result. 
data = assign_forces(data,{'RCAL','LCAL'},{'calcn_r','calcn_l'},[20 0.2],data.FilterFreq);

% create the trc and grf.mot files
data = btk_c3d2trc(data,'off');

%%
% Determine units of marker data
if strcmp(data.marker_data.Info.units.ALLMARKERS,'mm')
     p_sc = 1000;
     data.marker_data.Info.units.ALLMARKERS = 'm';
else p_sc = 1;
end


%% Write motion file containing GRFs



disp('Writing c3d file...')

if isfield(data,'fp_data')
     
     F = data.fp_data.Info(1).frequency/data.marker_data.Info.frequency; % assume that all force plates are collected at the same frequency!!!
     
     fp_time = 1/data.marker_data.Info.frequency:1/data.fp_data.Info(1).frequency:(F*(data.End_Frame-data.Start_Frame+1))/data.fp_data.Info(1).frequency;
     
     % initialise force data matrix with the time array and column header
     force_data_out = fp_time';
     force_header = 'time\t';
     force_format = '%20.6f\t';
     
     % if the assign_forces function has not been run then just make a force
     % file that contains the data from each force plate
     if ~isfield(data,'GRF')
          
          % go through each marker field and re-order from X Y Z to Y Z X and place
          % into data array and add data to the force data matrix --> also need to
          % divide by 1000 to convert to mm from m if necessary and Nmm to Nm
          % these are the conversions usually used in most motion analysis systems
          % and if they are different, just change the scale factor below to p_sc
          % value for the M data to 1. It should however get this from the file.
          
          for i = 1:length(data.fp_data.FP_data)
               
               % rescale the GRF COP to meters if necessary
               data.fp_data.GRF_data(i).P =  data.fp_data.GRF_data(i).P/p_sc;
               % rescale the moments to meters as well, if necessary
               data.fp_data.GRF_data(i).M =  data.fp_data.GRF_data(i).M/p_sc;
               
               % do some cleaning of the COP before and after contact
               b = find(abs(diff(data.fp_data.GRF_data(i).P(:,3)))>0);
               if ~isempty(b)
                    for j = 1:3
                         data.fp_data.GRF_data(i).P(1:b(1),j) = data.fp_data.GRF_data(i).P(b(1)+1,j);
                         data.fp_data.GRF_data(i).P(b(end):end,j) = data.fp_data.GRF_data(i).P(b(end)-1,j);
                    end
               end
               
               % define the period which we are analysing
               K = (F*data.Start_Frame):1:(F*data.End_Frame);
               
               % add the force, COP and moment data for current plate to the force matrix
               force_data_out = [force_data_out data.fp_data.GRF_data(i).F(K,:) data.fp_data.GRF_data(i).P(K,:) data.fp_data.GRF_data(i).M(K,:)];
               % define the header and formats
               force_header = [force_header num2str(i) '_ground_force_vx\t' num2str(i) '_ground_force_vy\t' num2str(i) '_ground_force_vz\t'...
                    num2str(i) '_ground_force_px\t' num2str(i) '_ground_force_py\t' num2str(i) '_ground_force_pz\t' ...
                    num2str(i) '_ground_torque_x\t' num2str(i) '_ground_torque_y\t' num2str(i) '_ground_torque_z\t'];
               force_format = [force_format '%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t'];
               
          end
          
          force_header = [force_header(1:end-2) '\n'];
          force_format = [force_format(1:end-2) '\n'];
          
          % assign a value of zero to any NaNs
          force_data_out(logical(isnan(force_data_out))) = 0;
          
     else
          
          % go through each marker field and re-order from X Y Z to Y Z X and place
          % into data array and add data to the force data matrix --> also need to
          % divide by 1000 to convert to mm from m if necessary and Nmm to Nm
          % these are the conversions usually used in most motion analysis systems
          % and if they are different, just change the scale factor below to p_sc
          % value for the M data to 1. It should however get this from the file.
          data.AssignForce.ExForce = [];
          data.AssignForce.ApBodies = [];
          
          for i = 1:length(data.GRF.FP)
               
               fnames = fieldnames(data.GRF.FP(i));
               
               for j = 1:length(fnames)
                    
                    if isfield(data.GRF.FP(i).(fnames{j}),'F')
                         % rescale the GRF COP and moments to meters if necessary
                         data.GRF.FP(i).(fnames{j}).P =  data.GRF.FP(i).(fnames{j}).P/p_sc;
                         data.GRF.FP(i).(fnames{j}).M =  data.GRF.FP(i).(fnames{j}).M/p_sc;
                         
                         % do some cleaning of the COP before and after contact
                         b = find(abs(diff(data.GRF.FP(i).(fnames{j}).P(:,3)))>0);
                         if ~isempty(b)
                              for k = 1:3
                                   data.GRF.FP(i).(fnames{j}).P(1:b(1),k) = data.GRF.FP(i).(fnames{j}).P(b(1)+1,k);
                                   data.GRF.FP(i).(fnames{j}).P(b(end):end,k) = data.GRF.FP(i).(fnames{j}).P(b(end)-1,k);
                              end
                         end
                         
                         % define the period which we are analysing
                         K = (F*data.Start_Frame):1:(F*data.End_Frame);
                        
                         % add the force, COP and moment data for current plate to the force matrix
                         force_data_out = [force_data_out data.GRF.FP(i).(fnames{j}).F(K,:) data.GRF.FP(i).(fnames{j}).P(K,:) data.GRF.FP(i).(fnames{j}).M(K,:)];
                         % define the header and formats
                         force_header = [force_header fnames{j} num2str(i) '_ground_force_vx\t' fnames{j} num2str(i) '_ground_force_vy\t' fnames{j} num2str(i) '_ground_force_vz\t'...
                              fnames{j} num2str(i) '_ground_force_px\t' fnames{j} num2str(i) '_ground_force_py\t' fnames{j} num2str(i) '_ground_force_pz\t' ...
                              fnames{j} num2str(i) '_ground_torque_x\t' fnames{j} num2str(i) '_ground_torque_y\t' fnames{j} num2str(i) '_ground_torque_z\t'];
                         force_format = [force_format '%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t'];
                         
                         data.AssignForce.ExForce = [data.AssignForce.ExForce;{[fnames{j} num2str(i)]}];
                         data.AssignForce.ApBodies = [data.AssignForce.ApBodies;{[fnames{j}]}];
                    end
               end
          end
          
          force_header = [force_header(1:end-2) '\n'];
          force_format = [force_format(1:end-2) '\n'];
          
          % assign a value of zero to any NaNs
          force_data_out(logical(isnan(force_data_out))) = 0;
     end
     
end

%% Edit the force plate data in Nexus

% REMEMBER TO LOAD TRIAL IN NEXUS FIRST
vicon = ViconNexus();

% Device ID:  1 = FP1, 2 = FP2, 3 = EMG
% Device Output ID: 1 = Force, 2 = Moment, 3 = CoP
% Channel ID: 1 = x, 2 = y, 3 = z

% Start with front plate
% Force
vicon.SetDeviceChannel(1,1,1,FPFront(:,1))
vicon.SetDeviceChannel(1,1,2,FPFront(:,2))
vicon.SetDeviceChannel(1,1,3,FPFront(:,3))

% COP
vicon.SetDeviceChannel(1,3,1,COPFront(:,1))
vicon.SetDeviceChannel(1,3,2,COPFront(:,2))
vicon.SetDeviceChannel(1,3,3,COPFront(:,3))

% Moment
vicon.SetDeviceChannel(1,2,1,MomentFront(:,1))
vicon.SetDeviceChannel(1,2,2,MomentFront(:,2))
vicon.SetDeviceChannel(1,2,3,MomentFront(:,3))

% Now back plate