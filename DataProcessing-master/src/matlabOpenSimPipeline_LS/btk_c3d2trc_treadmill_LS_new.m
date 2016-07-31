function [data, force_data2] = btk_c3d2trc_treadmill_LS_new(varargin)
% function btk_c3d2trc_treadmill(file) OR
% function btk_c3d2trc_treadmill(data)
%
% Function to convert data from a C3D file into the TRC and MOT file
% formats for OpenSim when using an AMTI treadmill, where the force
% assignments need to be adjusted to create the GRF mot file.
%
% INPUT -   file - the C3D file path that you wish to load (leave blank to
%               choose from a dialog box) OR
%           data - structure containing fields from from previously loaded
%               C3D file using btk_loadc3d.m
%           anim - animate 'on' or 'off' (default - 'on')
%           croppedTrialNum - number of loop used to name the output file
%           folder
%
% OUTPUT -  data - structure containing the relevant data from the c3dfile
%                  Creates the TRC file and _grf.MOT file for OpenSim
%
% example - data = btk_c3dtrc('filein.c3d','off');
%           data = btk_c3dtrc(data,'on');
%
% Written by Glen Lichtwark (University of Queensland)
% Updated September 2012

% Modified by Gavin Lenton for DST Group load sharing project
% Updated July 2016

%% load data
if nargin > 0
	if ~isstruct(varargin{1})
		% load C3d file
		file = varargin{1};
		if isempty(fileparts(file))
			pname = cd;
			if ispc
				pname = [pname '\'];
			else pname = [pname '/'];
			end
			fname = file;
		else [pname, name, ext] = fileparts(file);
			fname = [name ext];
		end
		% load the c3dfile
		data = btk_loadc3d([pname, fname], 10);
		
	else
		data = varargin{1};
		if ~isfield(data,'marker_data')
			error('Please ensure that the following field is included in structure - marker_data. Please use btk_loadc3d for correct outputs');
		end
		if isfield(data,'marker_data')
			[pname, name, ext] = fileparts(data.marker_data.Filename);
			if ispc
				pname = [pname '\'];
			else pname = [pname '/'];
			end
			fname = [name ext];
		else fname = data.marker_data.Filename;
		end
		
	end
	if length(varargin) < 2
		anim = 'on';
	else anim = varargin{2};
	end
else [fname, pname] = uigetfile('*.c3d', 'Select C3D file');
	% load the c3dfile
	data = btk_loadc3d([pname, fname], 10);
	anim = 'on';
end

%% define the start and end frame for analysis as first and last frame unless
% this has already been done to change the analysed frames
if ~isfield(data,'Start_Frame')
	data.Start_Frame = 1;
	data.End_Frame = data.marker_data.Info.NumFrames;
end

%%
% define some parameters
nrows = data.End_Frame-data.Start_Frame+1;
nmarkers = length(fieldnames(data.marker_data.Markers));

data.time = (1/data.marker_data.Info.frequency:1/data.marker_data.Info.frequency:(data.End_Frame-data.Start_Frame+1)/data.marker_data.Info.frequency)';

nframe = 1:nrows;

% anim the trial if animation = on
if strcmp(anim,'on')
	data.marker_data.First_Frame = data.Start_Frame;
	data.marker_data.Last_Frame = data.End_Frame;
	if isfield(data,'fp_data')
		btk_animate_markers(data.marker_data, data.fp_data, 5)
	else btk_animate_markers(data.marker_data)
	end
end

%%
% we need to reorder the lab coordinate system to match that of the OpenSim
% system --> SKIP THIS STEP IF LAB COORDINATE SYSTEM IS SAME AS MODEL
% SYSTEM
markers = fieldnames(data.marker_data.Markers); % get markers names

if strcmp(data.marker_data.Info.units.ALLMARKERS,'mm')
	p_sc = 1000;
	data.marker_data.Info.units.ALLMARKERS = 'm';
else p_sc = 1;
end

% go through each marker field and re-order from X Y Z to Y Z X
for i = 1:nmarkers
	data.marker_data.Markers.(markers{i}) =  [data.marker_data.Markers.(markers{i})(:,2)...
		data.marker_data.Markers.(markers{i})(:,3) data.marker_data.Markers.(markers{i})(:,1)]/p_sc;
end

%% Write trc file containing marker data

% initialise the matrix that contains the data as a frame number and time row
data_out = [nframe; data.time'];

% each of the data columns (3 per marker) will be in floating format with a
% tab delimiter - also add to the data matrix
for i = 1:nmarkers
	
	% add 3 rows of data for the X Y Z coordinates of the current marker
	% first check for NaN's and fill with a linear interpolant - warn the
	% user of the gaps
	clear m
	m = find(isnan(data.marker_data.Markers.(markers{i})((data.Start_Frame:data.End_Frame),1))>0);
	if ~isempty(m);
		clear t d
		disp(['Warning -' markers{i} ' data missing in parts. Frames ' num2str(m(1)) '-'  num2str(m(end))])
		t = time;
		t(m) = [];
		d = data.marker_data.Markers.(markers{i})((data.Start_Frame:data.End_Frame),:);
		d(m,:) = [];
		data.marker_data.Markers.(markers{i})((data.Start_Frame:data.End_Frame),:) = interp1(t,d,time,'linear','extrap');
	end
	data_out = [data_out; data.marker_data.Markers.(markers{i})((data.Start_Frame:data.End_Frame),:)'];
end

% Assign marker info to markersData
markersData = data_out';

% Apply filter here
data_out_filtered = zeros(size(markersData));

% Sampling freq
Fs = 1/data.marker_data.Info.frequency;

% Filter freq
Fc = 8;

% Apply 2nd order Butterworth filt with two passes.
for col = 3:length(data_out)
	data_out_filtered(:,col) = lpfilter(markersData(:,col),Fc,Fs, 'butter', 2);
	markersData(:,col) = data_out_filtered(:,col);
end

%% Define path and file names
indexName = regexp(fname(1:end-4), 'd\d*');
fileNameTRC = [fname(1:end-4), '.trc'];
fileNameGRF = regexprep(fileNameTRC, '.trc', '_grf.mot');

% Define path name to new folder
newpathname = [strrep(pname, 'InputData', 'ElaboratedData'),...
	'dynamicElaborations', filesep, fileNameTRC(1:indexName)];



% Final path name to store .trc and .mot files
finalpathname = [newpathname, filesep, fileNameTRC(1:end-4)];

% Add folder for the condition and walking speed in session
if ~exist(finalpathname, 'dir')
	mkdir(finalpathname);
end

fullFileNameTRC = [finalpathname filesep fileNameTRC];

% Marker labels names
MLabels = markers;

% Write the trc file
writetrc_LS(markersData,MLabels,data.marker_data.Info.frequency,fullFileNameTRC);

%% Write motion file containing GRFs

if isfield(data,'fp_data')
	
	Fp_change = data.fp_data.Info(1).frequency/data.marker_data.Info.frequency; % assume that all force plates are collected at the same frequency!!!
	
	% Changed this to start from frame 1, not from frame 10.
	%      fp_time = 1/data.marker_data.Info.frequency:1/data.fp_data.Info(1).frequency:(F*(data.End_Frame-data.Start_Frame+1))/data.fp_data.Info(1).frequency;
	
	fp_time1 = 0.001:1/data.fp_data.Info(1).frequency:(Fp_change*(data.End_Frame-data.Start_Frame+1))/data.fp_data.Info(1).frequency;
	
	% initialise force data matrix with the time array
	force_data_out = fp_time1';
	
	% add the force, COP and moment data for current plate to the force matrix
	% Needs to loop through both force plates
	for i = 1:2
		
		% rescale the GRF COP to meters if necessary
		data.fp_data.GRF_data(i).P =  data.fp_data.GRF_data(i).P/p_sc;
		% rescale the moments to meters as well, if necessary
		data.fp_data.GRF_data(i).M =  data.fp_data.GRF_data(i).M/p_sc;
		
		% do some cleaning of the COP before and after contact
		b = find(abs(diff(data.fp_data.GRF_data(i).P(:,2)))>0);
		if ~isempty(b)
			for j = 1:3
				data.fp_data.GRF_data(i).P(1:b(1),j) = data.fp_data.GRF_data(i).P(b(1)+1,j);
				data.fp_data.GRF_data(i).P(b(end):end,j) = data.fp_data.GRF_data(i).P(b(end)-1,j);
			end
		end
		
		% Define the period which we are analysing
		
		% Modified K to start from first frame and not start from first
		% frame * 10.
		K = (data.Start_Frame):1:(Fp_change*data.End_Frame);
		
		% reorder data so lab coordinate system to match that of the OpenSim
		% system
		if ~isempty(fieldnames(data.GRF.FP(i)))
			data.GRF.FP(i).P =  [data.GRF.FP(i).P(:,2)/p_sc ...
				data.GRF.FP(i).P(:,3)/p_sc data.GRF.FP(i).P(:,1)/p_sc];
			data.GRF.FP(i).F =  [data.GRF.FP(i).F(:,2) ...
				data.GRF.FP(i).F(:,3) data.GRF.FP(i).F(:,1)];
			data.GRF.FP(i).M =  [data.GRF.FP(i).M(:,2) ...
				data.GRF.FP(i).M(:,3) data.GRF.FP(i).M(:,1)]/p_sc;
			
			% If plate one then it's the right foot data
			
			if i == 1
				force_data_out = [force_data_out, data.GRF.FP(i).F(K,:)...
					data.GRF.FP(i).P(K,:) data.GRF.FP(i).M(K,:)];
				
				% If plate two it's the left foot data
			elseif i == 2
				
				force_data_out = [force_data_out, data.GRF.FP(i).F(K,:)...
					data.GRF.FP(i).P(K,:) data.GRF.FP(i).M(K,:)];
			else
				fprintf('GRF data does not exist for trial %s', fname);
			end
			
		end
		
	end
	
	% Apply 10 Hz filter to smooth stiching out - don't want to filter the
	% COP here.
	
	
	% Clean up data at beginning and end
	a = find(force_data_out(:,5) > 0);
	a2 = find(force_data_out(300:end,14) > 0) + 299;
	a3 = find(force_data_out(:,3) > 0);
	a4 = find(force_data_out(300:end,12) > 0) + 299;
	
	% Define filter parameters
	filt_freq = 20;
	dt = 1/data.fp_data.Info(1).frequency;
	
	force_data_filtered = zeros(size(force_data_out));
	
	% Apply 20th order (5 passes of a second order) critically-damped zero-lag filter
	for col1 = 1:3
		% Front plate forces
		force_data_filtered(a3(1):a3(end)+5,col1+1) = lpfilter(force_data_out(a3(1):a3(end)+5,col1+1),filt_freq,dt, 'damped', 10);
		% Front plate moments
		force_data_filtered(a3(1):a3(end)+5,col1+7) = lpfilter(force_data_out(a3(1):a3(end)+5,col1+7),filt_freq,dt, 'damped', 10);
		% Rear plate forces
		force_data_filtered(a4,col1+10) = lpfilter(force_data_out(a4,col1+10),filt_freq,dt, 'damped', 10);
		% Rear plate moments
		force_data_filtered(a4,col1+16) = lpfilter(force_data_out(a4,col1+16),filt_freq,dt, 'damped', 10);
		
		% COP front
		if a(end) ~= length(fp_time1)
		force_data_filtered(a(1)+1:a(end)+4, col1+4) = lpfilter(force_data_out(a(1)+1:a(end)+4, col1+4), 6,dt, 'butter', 2);
		else
			force_data_filtered(a(1)+1:a(end), col1+4) = lpfilter(force_data_out(a(1)+1:a(end), col1+4), 6,dt, 'butter', 2);
		end
		% COP rear
		force_data_filtered(a2, col1+13) = lpfilter(force_data_out(a2, col1+13), 6,dt, 'butter', 2);
	end
	
	% Butterworth filter to smooth out FP cropping
	b_filt_freq = 8;
	[pks, loc1] = findpeaks(force_data_filtered(:,3));
	force_data_filtered(a3(1):a3(end)+5, 2:4) = matfiltfilt(dt,b_filt_freq, 2, force_data_filtered(a3(1):a3(end)+5, 2:4));
	force_data_filtered(a3(1):a3(end)+5, 8:10) = matfiltfilt(dt,b_filt_freq,2, force_data_filtered(a3(1):a3(end)+5, 8:10));
	force_data_filtered(a4, 11:13) = matfiltfilt(dt,b_filt_freq, 2, force_data_filtered(a4, 11:13));
	force_data_filtered(a4, 17:19) = matfiltfilt(dt,b_filt_freq, 2, force_data_filtered(a4, 17:19));
	
	% Clean up COP
	if a(end) ~= length(fp_time1)
	force_data_filtered(a(1)+1:a(end)+5, 5:7) = matfiltfilt(dt,8,2, force_data_filtered(a(1)+1:a(end)+5, 5:7));
	else
		force_data_filtered(a(1)+1:a(end), 5:7) = matfiltfilt(dt,8,2, force_data_filtered(a(1)+1:a(end), 5:7));
	end
	force_data_filtered(a(1):a(1)+1, 7) = 0;
	force_data_filtered(a2, 14:16) = matfiltfilt(dt,8,2, force_data_filtered(a2, 14:16));
	
	% assign a value of zero to any NaNs
	force_data_filtered(logical(isnan(force_data_filtered))) = 0;
	
	% Re-arrange so data matches MOtoNMS convention
	force_data2 = force_data_filtered(:, 2:19);
	
	force_data2(:,7:12) = force_data_filtered(:,11:16);
	force_data2(:,13:15) = force_data_filtered(:,8:10);
	force_data2(:,16:18) = force_data_filtered(:,17:19);
	
	%% Find where force on left leg is zeroing and fix
	salame = find(force_data2(600:end, 8) == 0) + 599;
	
	% Only fix if there's a zero value
	if ~isempty(salame)
		if length(salame) > 1
			% If there's more than 1 frame zeroed then apply interp
			for columnM = 1:3
				force_data2(salame(1):salame(end), columnM+6) = (force_data2(salame(1)-1, columnM+6)...
					+ force_data2(salame(end)+1, columnM+6))/2;
				force_data2(salame(1):salame(end), columnM+15) = (force_data2(salame(1)-1, columnM+15)...
					+ force_data2(salame(end)+1, columnM+15))/2;
			end
			% Loop through columns and fix by taking mean of previous and
			% following frame
		else
			for columnF = 1:3
				force_data2(salame, columnF+6) = (force_data2(salame-1, columnF+6)...
					+ force_data2(salame+1, columnF+6))/2;
				force_data2(salame, columnF+15) = (force_data2(salame-1, columnF+15)...
					+ force_data2(salame+1, columnF+15))/2;
			end
		end
	end
	
	%% Print MOT
	
	if any(isempty(fieldnames(data.GRF.FP(i))))
		% Specify new file name if there is missing data name so I know to
		% check data
		disp('Trial is missing data, GRFs not printed')
		
	elseif any(force_data2(loc1(1):loc1(end), 2) < 50)
		% If there is issue with force assignment then print with modified
		% name
		disp('Trial has dodgy data, printing with modified filename');
		fullFileNameGRF = [finalpathname, filesep, fileNameGRF(1:end-4), '_NFU.mot'];
		% Write the MOT file using MOtoNMS function
		writeMot_LS(force_data2 ,force_data_out(:,1), fullFileNameGRF);
		
	else
		
		% Otherwise name and print normally
		fullFileNameGRF = [finalpathname, filesep, fileNameGRF];
		
		% Write the MOT file using MOtoNMS function
		writeMot_LS(force_data2 ,force_data_out(:,1), fullFileNameGRF);
	end
	
	cd ..
	
else disp('No force plate information available.')
end