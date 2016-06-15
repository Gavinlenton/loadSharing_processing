function [] = emgProcessingLS(isNotch, sessionData, timePoints, sessionName, isMax, maxc3dFileName)
%Input if Notch is required (Yes, No), session Data, timePoints from
%cropping for processing of EMG data, and sessionName (e.g., fast or slow
%walking). isMax specifies whether the EMG is used as max for normalisation
%   The Notch filter will only be applied for data collected directly
%   through Nexus. The sessionData is from the MOtoNMS c3d2mat function

%% EMG scaling code
% Designed to work with Load sharing data set
% Code searchers through user provided EMG mat files, process the signal
% and search for it global max.
% Code then stores all max emg values to file using CEINMS standard format.
% Code then asks for trials to scale the EMGs against max, and print to
% CEINMS standard format.

% By David Saxby, October 18th, 2015, d.saxby@griffith.edu.au
% Edited by Gavin Lenton, May 2016

%% EMG processing settings
videoSamplingRate = 100;

% Sampling rate
emgSamplingRate = 1000;
a2v = emgSamplingRate/videoSamplingRate;
emgSamplingTimeStep = 1/emgSamplingRate;
emgNyquistLimit = emgSamplingRate/2;
nominalPassBand = [30, 450]; % Need to double check what the upper pass band should be
passBandLowerEdge = nominalPassBand(1)/emgNyquistLimit;
passBandUpperEdge = nominalPassBand(2)/emgNyquistLimit;
normalizedPassBand = [passBandLowerEdge, passBandUpperEdge];
% high- and low-pass cut-off frequencies for band-pass filter
% N.B. In pass-band designs, cutoff frequencies are specified
% as a proportion of the Nyquist limit, i.e. rad/sample, not as a scalar.
passBandFilterOrder = 10; % this is double in matlab's conversion
lowPassCutOff = 6; % low-pass cut-off frequency for linear envelope
filterType = 'damped'; % Critically damped filter
emgPrintTag = '.mot'; % file format for printing EMG files for CEINMS

% Notch filter settings to remove 50 Hz (+harmonics) electromagnetic
% interference in signal. *ONLY APPLY TO DATA COLLECTED THROUGH NEXUS*
% Design
fo = 50; % nominal tone
q = 50; % filter "quality factor". Related to bandwidth of notch filter by q = fo/bw
bw = (fo/emgNyquistLimit)/q;
notchFilterOrder = (emgNyquistLimit/fo)*2; % specifies a filter with n+1 coefficients (numerator and demoninator).
% The filter will have n notches distributed from -1 to 1 normalized
% frequency domain

%% Generate pass band filter coefficents - change passBand for 1000Hz data
%  collection otherwise highPass will be = 1 and error thrown
[b_bp, a_bp] = butter(passBandFilterOrder,normalizedPassBand,'bandpass');

%% Optional, uncomment to inspect the band pass filter design
% plot filter magnitude and phase response
% freqz(b_bp,a_bp)

%% Generate notch filter
notchDesign = fdesign.comb('notch', 'N,BW', notchFilterOrder, bw);
notchFilter = design(notchDesign);

% Alternative filter formulation
% [b_notch, a_notch] = iircomb(emgSamplingRate/fo, bw, 'notch');

%% Optional, uncomment to inspect the notch filter design
% fvtool(notchFilter);

% Alternative notch visualization
% fvtool(b_notch, a_notch);

%% EMG labels in the .mat file - MAY BE SOURCE OF ERROR BECAUSE CHANNEL2 AND CHANNEL5 WERE SKIPPED
emgLabelsInMatFile = {'TA', 'Channel2',  'MG', 'LG', 'Channel5', 'BF', 'VM', 'VL',...
     'RF', 'Sol', 'MH'};
emgLabelsCEINMS = {'tibant_r', 'Channel2', 'gasmed_r', 'gaslat_r', 'Channel5' 'bicfemlh_r', 'vasmed_r', ...
     'vaslat_r', 'recfem_r', 'sol_r', 'semimem_r'};

%% Muscle groups
quadsLabels = {'rf', 'vaslat', 'vasmed'};
hamsLabels = {'bifem', 'semimem'};
medialLabels = {'medgas', 'semimem', 'vasmed'};
lateralLabels = {'latgas', 'bifem', 'vaslat'};

%% Gait sub-phases as defined in Sturnieks et al 2011
preContactTime = 0.05; % 50 ms prior to heelstrike
preContactInEMGFrames = 0.05*emgSamplingRate;
gaitPhases = {'preContact', 'loadingStance', 'earlyStance', 'midStance', 'terminalStance'};

%% EMG measures
emgMeasures = {'directedCoContraction', 'totalActivation', 'individualMuscles'};

%% Pairing
musclePairs = {'quadAndHams' , 'flexorsVsExtensors', 'medialToLateralwGastrocs', 'medialToLateralNoGastrocs'};

%% Directories and files
motoNmsSrcDir = 'C:\Users\s2921887\Documents\Load Sharing Main Data Collection\MOtoNMS v2.2\src';
% trialsFile = 'C:\Users\s2790936\Desktop\SCOPEX\SCOPEX Trials DJS.xlsm';
% maxDir = 'C:\Users\s2790936\Desktop\SCOPEXTemp\Baseline'; % Max's

originalPath=pwd;
cd(motoNmsSrcDir)
cd('shared')
sharedFunctionsPath=pwd;
addpath(sharedFunctionsPath)
cd(originalPath)

%% Load session data
analogData = [sessionData, filesep, sessionName]; % Dynamic subject dir
maxData = [sessionData, filesep, maxc3dFileName]; % Max trial dir

%% Make empty summary structure for end of processing
for p = 1:length(gaitPhases)
     for eM = 1:length(emgMeasures)
          if eM ~= 3
               for mP = 1:length(musclePairs)
                    emgMetrics.('means').('FastWalking').(gaitPhases{p}).(emgMeasures{eM}).(musclePairs{mP}) = [];
                    emgMetrics.('means').('NormalWalking').(gaitPhases{p}).(emgMeasures{eM}).(musclePairs{mP}) = [];
                    emgMetrics.('stats').('FastWalking').(gaitPhases{p}).(emgMeasures{eM}).(musclePairs{mP}) = [];
                    emgMetrics.('stats').('NormalWalking').(gaitPhases{p}).(emgMeasures{eM}).(musclePairs{mP}) = [];
                    emgMetrics.('stack').('FastWalking').(gaitPhases{p}).(emgMeasures{eM}).(musclePairs{mP}) = [];
                    emgMetrics.('stack').('NormalWalking').(gaitPhases{p}).(emgMeasures{eM}).(musclePairs{mP}) = [];
               end
          else
               for iM = 1:length(emgLabelsCEINMS)
                    emgMetrics.('means').('FastWalking').(gaitPhases{p}).(emgMeasures{eM}).(emgLabelsCEINMS{iM}) = [];
                    emgMetrics.('means').('NormalWalking').(gaitPhases{p}).(emgMeasures{eM}).(emgLabelsCEINMS{iM}) = [];
                    emgMetrics.('stats').('FastWalking').(gaitPhases{p}).(emgMeasures{eM}).(emgLabelsCEINMS{iM}) = [];
                    emgMetrics.('stats').('NormalWalking').(gaitPhases{p}).(emgMeasures{eM}).(emgLabelsCEINMS{iM}) = [];
                    emgMetrics.('stack').('FastWalking').(gaitPhases{p}).(emgMeasures{eM}).(emgLabelsCEINMS{iM}) = [];
                    emgMetrics.('stack').('NormalWalking').(gaitPhases{p}).(emgMeasures{eM}).(emgLabelsCEINMS{iM}) = [];
               end
          end
     end
end

%% Specify emg.mat file

% Specify for dynamic trial
suffix = 'AnalogData.mat';
emgFile1 = fullfile(analogData, suffix);

% Specify for max trial
filenameMax = fullfile(maxData,  suffix);

%% Load emg mat files

% File 1
try
     emgData1 = load(emgFile1);
catch me1
     disp(['No such file ', emgFile1, ' , continuing with analysis'])
end

%% Initialize some values
try
     channels = fieldnames(emgData1);
catch me
     disp(['No such file ', emgFile1, ' , continuing with analysis'])
end
emgMax = zeros(2, length(emgLabelsInMatFile));

%% Condition raw emg signal
if exist('emgData1', 'var') && exist('emgData2', 'var') && exist('emgData3', 'var')
     nDataPoints = 1;
     filesUsedForMax = [1,2,3];
     for xx = 1:length(channels)
          lengthOfData = size([emgData1.(channels{xx}).values ; emgData2.(channels{xx}).values ; emgData3.(channels{xx}).values],1);
          if lengthOfData > nDataPoints
               nDataPoints = lengthOfData;
          end
     end
elseif exist('emgData1', 'var') && exist('emgData2', 'var')
     nDataPoints = 1;
     filesUsedForMax = [1,2];
     for xx = 1:length(channels)
          lengthOfData = size([emgData1.(channels{xx}).values ; emgData2.(channels{xx}).values],1);
          if lengthOfData > nDataPoints
               nDataPoints = lengthOfData;
          end
     end
elseif exist('emgData1', 'var') && exist('emgData3', 'var')
     nDataPoints = 1;
     filesUsedForMax = [1,3];
     for xx = 1:length(channels)
          lengthOfData = size([emgData1.(channels{xx}).values ; emgData3.(channels{xx}).values],1);
          if lengthOfData > nDataPoints
               nDataPoints = lengthOfData;
          end
     end
elseif exist('emgData2', 'var') && exist('emgData3', 'var')
     nDataPoints = 1;
     filesUsedForMax = [2,3];
     for xx = 1:length(channels)
          lengthOfData = size([emgData2.(channels{xx}).values ; emgData3.(channels{xx}).values],1);
          if lengthOfData > nDataPoints
               nDataPoints = lengthOfData;
          end
     end
elseif exist('emgData1', 'var')
     nDataPoints = 1;
     filesUsedForMax = 1;
     for xx = 1:length(channels)
          lengthOfData = size(emgData1.(channels{xx}).RawData, 1);
          if lengthOfData > nDataPoints
               nDataPoints = lengthOfData;
          end
     end
elseif exist('emgData2', 'var')
     nDataPoints = 1;
     filesUsedForMax = 2;
     for xx = 1:length(channels)
          lengthOfData = size(emgData2.(channels{xx}).values, 1);
          if lengthOfData > nDataPoints
               nDataPoints = lengthOfData;
          end
     end
elseif exist('emgData3', 'var')
     nDataPoints = 1;
     filesUsedForMax = 3;
     for xx = 1:length(channels)
          lengthOfData = size(emgData3.(channels{xx}).values, 1);
          if lengthOfData > nDataPoints
               nDataPoints = lengthOfData;
          end
     end
end

%% Initialize matrix
emgData = zeros(nDataPoints, size(emgLabelsInMatFile,2));


for n = 1:length(emgLabelsInMatFile)
     for x = 1:size(channels,1)
          if exist('emgData1', 'var') && exist('emgData2', 'var') && exist('emgData3', 'var')
               if strcmp(emgData1.(channels{x}).title, emgLabelsInMatFile{n})
                    % Concatenate the data
                    try
                         emgData(:,n) = [emgData1.(channels{x}).values ; emgData2.(channels{x}).values ; emgData3.(channels{x}).values];
                         break;
                    catch me
                         tempData = [emgData1.(channels{x}).values ; emgData2.(channels{x}).values ; emgData3.(channels{x}).values];
                         if length(tempData) > size(emgData,1)
                              emgData(:,n) = tempData(1:size(emgData,1),:);
                         else
                              emgData(:,n) = [tempData ; zeros(size(emgData,1) - length(tempData), 1)];
                         end
                         break;
                    end
               end
          elseif exist('emgData1', 'var') && exist('emgData2', 'var')
               if strcmp(emgData1.(channels{x}).title, emgLabelsInMatFile{n})
                    % Concatenate the data
                    try
                         emgData(:,n) = [emgData1.(channels{x}).values ; emgData2.(channels{x}).values];
                         break;
                    catch me
                         tempData = [emgData1.(channels{x}).values ; emgData2.(channels{x}).values];
                         if length(tempData) > size(emgData,1)
                              emgData(:,n) = tempData(1:size(emgData,1),:);
                         else
                              emgData(:,n) = [tempData ; zeros(size(emgData,1) - length(tempData), 1)];
                         end
                         break;
                    end
               end
          elseif exist('emgData1', 'var') && exist('emgData3', 'var')
               if strcmp(emgData1.(channels{x}).title, emgLabelsInMatFile{n})
                    % Concatenate the data
                    try
                         emgData(:,n) = [emgData1.(channels{x}).values ; emgData3.(channels{x}).values];
                         break;
                    catch me
                         tempData = [emgData1.(channels{x}).values ; emgData3.(channels{x}).values];
                         if length(tempData) > size(emgData,1)
                              emgData(:,n) = tempData(1:size(emgData,1),:);
                         else
                              emgData(:,n) = [tempData ; zeros(size(emgData,1) - length(tempData), 1)];
                         end
                         break;
                    end
               end
          elseif exist('emgData2', 'var') && exist('emgData3', 'var')
               if strcmp(emgData2.(channels{x}).title, emgLabelsInMatFile{n})
                    % Concatenate the data
                    try
                         emgData(:,n) = [emgData2.(channels{x}).values ; emgData3.(channels{x}).values];
                         break;
                    catch me
                         tempData = [emgData2.(channels{x}).values ; emgData3.(channels{x}).values];
                         if length(tempData) > size(emgData,1)
                              emgData(:,n) = tempData(1:size(emgData,1),:);
                         else
                              emgData(:,n) = [tempData ; zeros(size(emgData,1) - length(tempData), 1)];
                         end
                         break;
                    end
               end
          elseif exist('emgData1', 'var')
               if strcmp(emgData1.(channels{x}).Labels{n}, emgLabelsInMatFile{n})
                    try
                         emgData(:,n) = emgData1.(channels{x}).RawData(:,n);
                         break;
                    catch me
                         tempData = emgData1.(channels{x}).RawData(:,n);
                         if length(tempData) > size(emgData,1)
                              emgData(:,n) = tempData(1:size(emgData,1),:);
                         else
                              emgData(:,n) = [tempData ; zeros(size(emgData,1) - length(tempData), 1)];
                         end
                         break;
                    end
               end
          elseif exist('emgData2', 'var')
               if strcmp(emgData2.(channels{x}).title, emgLabelsInMatFile{n})
                    try
                         emgData(:,n) = emgData2.(channels{x}).values;
                         break;
                    catch me
                         tempData = emgData2.(channels{x}).values;
                         if length(tempData) > size(emgData,1)
                              emgData(:,n) = tempData(1:size(emgData,1),:);
                         else
                              emgData(:,n) = [tempData ; zeros(size(emgData,1) - length(tempData), 1)];
                         end
                         break;
                    end
               end
          elseif exist('emgData3', 'var')
               if strcmp(emgData3.(channels{x}).title, emgLabelsInMatFile{n})
                    try
                         emgData(:,n) = emgData3.(channels{x}).values;
                         break;
                    catch me
                         tempData = emgData3.(channels{x}).values;
                         if length(tempData) > size(emgData,1)
                              emgData(:,n) = tempData(1:size(emgData,1),:);
                         else
                              emgData(:,n) = [tempData ; zeros(size(emgData,1) - length(tempData), 1)];
                         end
                         break;
                    end
               end
          end
     end
end

%% Find max signal and process only if trial is a Max trial

if isMax == 1
     
     MaxEMG_trials = cell(size(emgLabelsInMatFile,2),1);
     % signal dimensions
     lengthOfSignal = size(emgData,1);
     % define frequency domain
     f = emgSamplingRate*(0:(lengthOfSignal/2))/lengthOfSignal;
     
     for l = 1:size(emgLabelsInMatFile,2)
          
          % Explore trend/DC offsets in signal
          detrendedSignal = detrend(emgData(:,l), 'constant');
          
          % Band pass filter
          bandPassOut = filter(b_bp, a_bp, detrendedSignal);
          
          % Apply notch filter
          notchedSignal = filter(notchFilter, bandPassOut);
          
          % Fast fourier transform on the band-passed signal
          fastFourierBandPassSignal = fft(bandPassOut);
          
          % FFT on the notched signal as well
          fastFourierNotchedSignal = fft(notchedSignal);
          
          % Two sided spectrum of FFT non-notched
          P2 = abs(fastFourierBandPassSignal/lengthOfSignal);
          P1 = P2(1:round(lengthOfSignal/2)+1);
          P1(2:end-1) = 2*P1(2:end-1);
          
          % Two sided spectrum of FFT notched signal
          P2n = abs(fastFourierNotchedSignal/lengthOfSignal);
          P1n = P2n(1:round(lengthOfSignal/2)+1);
          P1n(2:end-1) = 2*P1n(2:end-1);
          
          if ~isequal(length(f), length(P1))
               f = emgSamplingRate*(0:((lengthOfSignal+1)/2))/(lengthOfSignal+1);
          end
          
          % Full-wave rectify
          if isNotch(strcmp('yes', 'yes')) == 1
               fullWaveRectifiedSignal = abs(notchedSignal);
          else
               fullWaveRectifiedSignal = abs(bandPassOut);
          end
          
          % Low-pass filter
          lowPassOut = lpfilter(fullWaveRectifiedSignal, lowPassCutOff, emgSamplingRate, filterType);
          
          % Find max
          [emgMax(1, l), emgMax(2, l)] = max(lowPassOut);
          
          % Trial (a or b) from which max was from
          if filesUsedForMax == 3
               
               if emgMax(2, l) > size(emgData1.(channels{1}).values, 1)
                    MaxEMG_trials{l} = filenameMax;
               else
                    MaxEMG_trials{l} = filename1;
               end
               
          elseif filesUsedForMax == 2
               
               MaxEMG_trials{l} = filenameMax;
               
          elseif filesUsedForMax == 1
               
               MaxEMG_trials{l} = filename1;
               
          end
          
          % Optional plotting of signal processing steps, uncomment to
          % view
          fig(1) = figure('Name', [analogData, ' MaxEfforts ', char(emgLabelsInMatFile{l})]);
          
          subplot(4,2,1);
          plot(emgData(:,l), 'k');
          box off
          legend('Raw Signal');
          legend boxoff;
          subplot(4,2,2);
          plot(detrendedSignal, 'g');
          box off
          legend('Detrended Signal');
          legend boxoff;
          subplot(4,2,3);
          plot(bandPassOut, 'r');
          box off
          legend('Band-passed Signal');
          legend boxoff;
          subplot(4,2,4);
          plot(f,P1, 'y');
          box off
          xlabel('frequency (Hz)');
          ylabel('|P1(f)|');
          legend('Single-sided amplitude spectrum');
          legend boxoff;
          subplot(4,2,5);
          plot(f,P1n, 'm');
          box off
          xlabel('frequency (Hz)');
          ylabel('|P1Notched (f)|');
          legend('Single-sided amplitude spectrum (notched)');
          legend boxoff;
          subplot(4,2,6)
          plot(fullWaveRectifiedSignal, 'b');
          box off
          legend('Full-wave Rectified Signal');
          legend boxoff;
          subplot(4,2,7)
          plot(lowPassOut, 'c');
          box off
          legend('Linear Envelope');
          legend boxoff;
          
     end
     
     sessionDataMax = [sessionData, filesep, sessionName];
     % Print max emg to MotoNMS directory
     printDir = [sessionDataMax, filesep, 'maxEMG'];
     if ~isdir(printDir)
          mkdir(printDir)
     end
     nEMGChannels=length(emgLabelsInMatFile);
     fid = fopen([printDir filesep 'maxEmg.txt'], 'w');
     fprintf(fid,'Muscle\tMaxEMGvalue\tTrial\tTime(s)\tFrame#\n');
     for i=1:nEMGChannels
          MaxEMGLabel = char(emgLabelsInMatFile{i});
          fprintf(fid,'%s\t%6.4e\t%s\t%6.4f\t%6.4f\n', MaxEMGLabel, emgMax(1,i), MaxEMG_trials{i}, (emgMax(2,i)/emgSamplingRate), emgMax(2,i));
     end
     fclose(fid);
     clearvars me
     
else
     % This is where I should include normal EMG processing steps for
     % dynamic trials - also load max data here
     
     
end
%% Dynamic signal conditioning

% If not a max trial, then load the max data and use for processing

if isMax == 0
     
     emgMaxDir = [sessionDataMax, filesep, 'maxEMG'];
     emgMaxFile = fullfile(emgMaxDir, 'maxEmg.txt');
     
else
end

analogData = sessionName;

%% Create empty structure for outputs
for p = 1:length(gaitPhases)
     for eM = 1:length(emgMeasures)
          if eM ~= 3
               for mP = 1:length(musclePairs)
                    emgMetrics.(analogData).(gaitPhases{p}).(emgMeasures{eM}).(musclePairs{mP}) = [];
                    subjectEmgMetrics.(analogData).(gaitPhases{p}).(emgMeasures{eM}).(musclePairs{mP}) = [];
               end
          else
               for iM = 1:length(emgLabelsCEINMS)
                    emgMetrics.(analogData).(gaitPhases{p}).(emgMeasures{eM}).(emgLabelsCEINMS{iM}) = [];
                    subjectEmgMetrics.(analogData).(gaitPhases{p}).(emgMeasures{eM}).(emgLabelsCEINMS{iM}) = [];
               end
          end
     end
end

%% Loop through gait cycles

% Load session data analog file
analogFile = fullfile([sessionData, filesep, analogData, filesep, 'AnalogData.mat']);
try
     load(analogFile);
catch me
     error(['Cannot load analog data mat file for ' , analogData, ' ensure this data has been converted from c3d to mat'])
end
nFrames = size(AnalogData.RawData,1);

% Loop through gait cycles

for cycle = 1:length(timePoints)-1
          
     % Define event times
          heelStrikeTime = timePoints(cycle,1);
          toeOffTime = timePoints(cycle, 2);
          heelStrike2Time = timePoints(cycle+1,1);
          
          % Correct if first Toe-off occured after first heel-strike
          if toeOffTime < heelStrikeTime
               toeOffTime = timePoints(cycle+1,2);
          else
          end
          
          % correct for any cropping
          if   heelStrikeTime ~= 1 && heelStrikeTime > 5
               
               heelStrikeInEMGFrame = heelStrikeTime*a2v;
               toeOffInEMGFrame = toeOffTime*a2v;
               heelStrike2InEMGFrame = heelStrike2Time*a2v;
               
          elseif heelStrikeTime <= 5 
               
               try
               heelStrikeTime = timePoints(cycle+1,1);
               toeOffTime = timePoints(cycle+1, 2);
               heelStrike2Time = timePoints(cycle+2,1);
               
               % Error if timePoints index is wrong
               catch err
               
                    disp('No more heel strikes to work with')
                    break
                    
               end 
               heelStrikeInEMGFrame = heelStrikeTime*a2v;
               toeOffInEMGFrame = toeOffTime*a2v;
               heelStrike2InEMGFrame = heelStrike2Time*a2v;
               
          else  
               toeOffInEMGFrame = toeOffTime*a2v;
               heelStrike2InEMGFrame = heelStrike2Time*a2v;
               
          end
          
     
     % Define gait cycle length
     gaitCycleLength = ((heelStrikeInEMGFrame-preContactInEMGFrames):heelStrike2InEMGFrame);
     startAt = gaitCycleLength(1);
     endAt = gaitCycleLength(end);
     
     % Correct for cropping on time
     if AnalogData.FirstFrame ~= 1
          time = (AnalogData.FirstFrame*a2v)*emgSamplingTimeStep : emgSamplingTimeStep : (nFrames-1+(AnalogData.FirstFrame*a2v))*emgSamplingTimeStep;
     else
          time = 0 : emgSamplingTimeStep : (nFrames-1)*emgSamplingTimeStep;
     end
     
     % Crop time to be from heelStrikeInEMGFrame - preContactTime to
     % heelStrike2InEMGFrame
     timeGaitCycle = time(heelStrikeInEMGFrame-40):emgSamplingTimeStep:time(heelStrike2InEMGFrame+10);
     
     % Initialise
     dynamicTrialScaledEmg = zeros(length(gaitCycleLength), length(emgLabelsInMatFile));
     
     % Map analog channels to proper EMG names
     dynamicEMGMatrix = [AnalogData.RawData(startAt:endAt, 1) , AnalogData.RawData(startAt:endAt, 2),...
          AnalogData.RawData(startAt:endAt, 3), AnalogData.RawData(startAt:endAt, 4),...
          AnalogData.RawData(startAt:endAt, 5), AnalogData.RawData(startAt:endAt, 6),...
          AnalogData.RawData(startAt:endAt, 7), AnalogData.RawData(startAt:endAt, 8), ...
          AnalogData.RawData(startAt:endAt, 9), AnalogData.RawData(startAt:endAt, 10), AnalogData.RawData(startAt:endAt,11)];
     
     % signal dimensions
     lengthOfSignal = size(dynamicEMGMatrix,1);
     % define frequency domain
     f = emgSamplingRate*(0:(lengthOfSignal/2))/lengthOfSignal;
     
     % Condition the dynamic EMG
     for n = 1:length(emgLabelsInMatFile)
          
          % De-trended
          detrendedSignal = detrend(dynamicEMGMatrix(:, n), 'constant');
          
          % Band pass filter
          bandPassOut = filter(b_bp, a_bp, detrendedSignal);
          
          % Apply notch filter
          notchedSignal = filter(notchFilter, bandPassOut);
          
          % Fast fourier transform on the band-passed signal
          fastFourierBandPassSignal = fft(bandPassOut);
          
          % FFT on the notched signal as well
          fastFourierNotchedSignal = fft(notchedSignal);
          
          % Two sided spectrum of FFT non-notched
          P2 = abs(fastFourierBandPassSignal/lengthOfSignal);
          P1 = P2(1:round(lengthOfSignal/2)+1);
          P1(2:end-1) = 2*P1(2:end-1);
          
          % Two sided spectrum of FFT notched signal
          P2n = abs(fastFourierNotchedSignal/lengthOfSignal);
          P1n = P2n(1:round(lengthOfSignal/2)+1);
          P1n(2:end-1) = 2*P1n(2:end-1);
          
          if ~isequal(length(f), length(P1))
               f = emgSamplingRate*(0:((lengthOfSignal+1)/2))/(lengthOfSignal+1);
          end
          
          % Full-wave rectify
          if strcmp(isNotch, 'yes') == 1
               fullWaveRectifiedSignal = abs(notchedSignal);
          else
               fullWaveRectifiedSignal = abs(bandPassOut);
          end
          
          % Low-pass filter
          lowPassOut = lpfilter(fullWaveRectifiedSignal, lowPassCutOff, emgSamplingRate, filterType);
          
          if isMax == 1
               % Find max 
               dynamicTrialScaledEmg(:, n) = lowPassOut/emgMax(1,n);
               
          else
               dynamicTrialScaledEmg(:, n) = lowPassOut;
          end
          
          % Check if the linear envelope has scaled signal > 1
%           if max(dynamicTrialScaledEmg(:, n))>1
%                display(['The maximum scaled EMG signal for ' , char(emgLabelsInMatFile{n}), ' during ' , char(analogData), ' is greater than 1'])
%           end
%           
          % Optional plotting of signal processing steps, uncomment to
          % view
          %                fig(1) = figure('Name', char(emgLabelsInMatFile{n}));
          %
          %                subplot(4,2,1);
          %                plot(dynamicEMGMatrix(:, n), 'k');
          %                box off
          %                legend('Raw Data');
          %                legend boxoff;
          %                subplot(4,2,2);
          %                plot(detrendedSignal, 'g');
          %                box off
          %                legend('Detrended');
          %                legend boxoff;
          %                subplot(4,2,3);
          %                plot(bandPassOut, 'r');
          %                box off
          %                legend('Band-passed');
          %                legend boxoff;
          %                subplot(4,2,4);
          %                plot(f,P1,'y');
          %                box off
          %                xlabel('frequency (Hz)');
          %                ylabel('|P1(f)|');
          %                legend('Single-sided amplitude spectrum');
          %                legend boxoff;
          %                subplot(4,2,5)
          %                plot(f,P1n,'m');
          %                box off
          %                xlabel('frequency (Hz)');
          %                ylabel('|P1 Notched (f)|');
          %                legend('Single-sided amplitude spectrum (notched)');
          %                legend boxoff;
          %                subplot(4,2,6)
          %                plot(fullWaveRectifiedSignal, 'b');
          %                box off
          %                legend('Full-wave Rectified');
          %                legend boxoff;
          %                subplot(4,2,7)
          %                plot(lowPassOut, 'c');
          %                box off
          %                legend('Linear Envelope');
          %                legend boxoff;
          %                subplot(4,2,8)
          %                plot(dynamicTrialScaledEmg(:, n), 'k');
          %                box off
          %                legend('Scaled Linear Envelope');
          %                legend boxoff;
          %
          %                set(gcf,'color','w');
          %                checkHandle1 = get(gcf,'Renderer');
          %                if strcmp(checkHandle1,'opengl')
          %                     set(gcf,'Renderer','painters');
          %                end
          
          if ~exist([sessionData(1:end-11) , filesep, 'dynamicElaborations',  filesep, analogData, filesep, 'EMGs'], 'dir')
               mkdir([sessionData(1:end-11) , filesep, 'dynamicElaborations',  filesep, analogData, filesep, 'EMGs\Envelope']);
               mkdir([sessionData(1:end-11) , filesep, 'dynamicElaborations',  filesep, analogData, filesep, 'EMGs\Raw']);
          end
          
          %                savefig(figure(1), [sessionData(1:end-11) , filesep, 'dynamicElaborations',  filesep, analogData, filesep, 'EMGs', filesep, char(emgLabelsInMatFile{n}), '.fig']);
          %                close all;
          
     end
     
     % Save EMGs to dynamic elaboration folder
     save([sessionData(1:end-11) , filesep, 'dynamicElaborations',  filesep, analogData, filesep, 'EMGs', filesep, 'Raw', filesep, ['EMGsSelectedRaw', num2str(cycle), '.mat']], 'emgLabelsInMatFile', 'dynamicEMGMatrix');
     save([sessionData(1:end-11) , filesep, 'dynamicElaborations',  filesep, analogData, filesep, 'EMGs', filesep, 'Envelope', filesep, ['EMGsSelectedEnvelope', num2str(cycle), '.mat']], 'emgLabelsInMatFile', 'dynamicTrialScaledEmg');
     
     % Print the scaled EMG envelopes to .mot format for future
     % CEINMS use.
     emgMotDir = [sessionData(1:end-11) , 'dynamicElaborations', filesep, analogData];
     printEMGmot_LS(emgMotDir,timeGaitCycle',dynamicTrialScaledEmg, emgLabelsInMatFile, emgPrintTag, cycle);
     
     % GAIT CYCLE ANALYSIS
     
     % full gait cycle
     stancePhaseInEMGFrame = toeOffInEMGFrame - heelStrikeInEMGFrame;
     fifteenPercent = round(15 * (stancePhaseInEMGFrame/100));
     thirtyfivePercent = round(35 * (stancePhaseInEMGFrame/100));
     twentyPercent = round(20 * (stancePhaseInEMGFrame/100));
     
     % Break the gait cycle into sub-phases
     preContactA = 1;
     pC = preContactA : preContactInEMGFrames; % Pre contact
     lS = (pC(end)+1) : (pC(end)+fifteenPercent); % loading stance
     eS = (lS(end)+1) : (lS(end)+thirtyfivePercent); % early stance
     mS = (eS(end)+1) : (eS(end)+twentyPercent); % mid stance
     tS = (mS(end)+1) : stancePhaseInEMGFrame;
     
     %Compute the EMG metrics for each phase
     for p = 1:length(gaitPhases)
          if p == 1
               phase = pC;
          elseif p == 2
               phase = lS;
          elseif p == 3
               phase = eS;
          elseif p == 4
               phase = mS;
          else
               phase = tS;
          end
          
          % individual muscles
          
          for eM = 3:length(emgMeasures)
               for iM = 1:length(emgLabelsCEINMS)
                    emgMetrics.(analogData).(gaitPhases{p}).(emgMeasures{eM}).(emgLabelsCEINMS{iM}) = dynamicTrialScaledEmg(phase, iM);
                    subjectEmgMetrics.(analogData).(gaitPhases{p}).(emgMeasures{eM}).(emgLabelsCEINMS{iM}) = dynamicTrialScaledEmg(phase, iM);
               end
          end
          
          % groups of muscles
          rf = dynamicTrialScaledEmg(phase, 5);
          vl = dynamicTrialScaledEmg(phase, 7);
          vm = dynamicTrialScaledEmg(phase, 8);
          quad = rf + vm + vl;
          bf = dynamicTrialScaledEmg(phase, 1);
          sm = dynamicTrialScaledEmg(phase, 6);
          hams = bf + sm;
          mg = dynamicTrialScaledEmg(phase, 3);
          lg = dynamicTrialScaledEmg(phase, 2);
          hamsPlus = hams + mg + lg;
          medial = vm + sm;
          lateral = vl + bf;
          medialPlus = medial + mg;
          lateralPlus = lateral + lg;
          
          for eM = 1:length(emgMeasures)-1
               for mP = 1:length(musclePairs)
                    if mP == 1 % quad vs hams
                         emgMetric = zeros(length(quad), 1);
                         if eM == 1 % directed co-contraction
                              for y = 1:length(quad)
                                   if quad(y,:)>hams(y,:)
                                        emgMetric(y,:) = 1 - (hams(y,:)/quad(y,:));
                                   else
                                        emgMetric(y,:) = (quad(y,:)/hams(y,:)) - 1;
                                   end
                              end
                         else % total activation
                              for y = 1:length(quad)
                                   emgMetric(y,:) = hams(y,:) + quad(y,:);
                              end
                         end
                    elseif mP == 2 % flexors vs extensors
                         emgMetric = zeros(length(quad), 1);
                         if eM == 1 % directed co-contraction
                              for y = 1:length(quad)
                                   if quad(y,:)>hamsPlus(y,:)
                                        emgMetric(y,:) = 1 - (hamsPlus(y,:)/quad(y,:));
                                   else
                                        emgMetric(y,:) = (quad(y,:)/hamsPlus(y,:)) - 1;
                                   end
                              end
                         else % total activation
                              for y = 1:length(quad)
                                   emgMetric(y,:) = hamsPlus(y,:) + quad(y,:);
                              end
                         end
                    elseif mP == 3 % medial to lateral with the gastrocs
                         emgMetric = zeros(length(medialPlus), 1);
                         if eM == 1 % directed co-contraction
                              for y = 1:length(medialPlus)
                                   if medialPlus(y,:)>lateralPlus(y,:)
                                        emgMetric(y,:) = 1 - (lateralPlus(y,:)/medialPlus(y,:));
                                   else
                                        emgMetric(y,:) = (medialPlus(y,:)/lateralPlus(y,:)) - 1;
                                   end
                              end
                         else % total activation
                              emgMetric(y,:) = lateralPlus(y,:) + medialPlus(y,:);
                         end
                         
                    elseif mP == 4 % medial to lateral with no gastrocs
                         emgMetric = zeros(length(medial), 1);
                         if eM == 1 % directed co-contraction
                              for y = 1:length(medial)
                                   if medial(y,:)>lateral(y,:)
                                        emgMetric(y,:) = 1 - (lateral(y,:)/medial(y,:));
                                   else
                                        emgMetric(y,:) = (medial(y,:)/lateral(y,:)) - 1;
                                   end
                              end
                         else % total activation
                              emgMetric(y,:) = lateral(y,:) + medial(y,:);
                         end
                    end
                    
                    % Passing to structure
                    emgMetrics.(analogData).(gaitPhases{p}).(emgMeasures{eM}).(musclePairs{mP}) = emgMetric;
                    subjectEmgMetrics.(analogData).(gaitPhases{p}).(emgMeasures{eM}).(musclePairs{mP}) = emgMetric;
               end
          end
     end
     save([sessionData(1:end-11) , filesep, 'dynamicElaborations',  filesep, analogData, filesep, 'EMGs', filesep, ['EMGsMetrics', num2str(cycle), '.mat']], 'emgMetrics')
end

% %% make means from the individual's stacks
% type=fieldnames(subjectEmgMetrics.(analogData));
% for t = 1:length(type)
%      for gP = 1:length(gaitPhases)
%           for emgM = 1:length(emgMeasures)
%                theMeasures = fieldnames(subjectEmgMetrics.(analogData).(gaitPhases{gP}).(emgMeasures{emgM}));
%                for tM = 1:length(theMeasures)
%                     if emgM ~= 3
%                          subjectEmgMetrics.(analogData).('means').(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}) = mean(subjectEmgMetrics.(analogData).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}));
%                     else % time-varying
%                          lOM = size(subjectEmgMetrics.(analogData).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}),2);
%                          sOM = size(subjectEmgMetrics.(analogData).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}),1);
%                          tempVar = zeros(1,lOM);
%                          for l = 1:lOM
%                               tempVar(1,l) = mean(subjectEmgMetrics.(analogData).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM})(:,l));
%                          end
%                          subjectEmgMetrics.(analogData).('means').(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}) = tempVar;
%                     end
%                end
%           end
%      end
% end
% 
% %% Save that person's data to a structure
% save([sessionData, filesep, analogData, filesep, [analogData, '_emgMetrics.mat']], 'subjectEmgMetrics');
% clearvars subjectEmgMetrics
% 
% %% Take basic stats on the summary stack
% for t = 1:length(type)
%      for gP = 1:length(gaitPhases)
%           for emgM = 1:length(emgMeasures)
%                theMeasures = fieldnames(emgMetrics.('stack').(type{t}).(gaitPhases{gP}).(emgMeasures{emgM}));
%                for tM = 1:length(theMeasures)
%                     if emgM ~= 3
%                          emgMetrics.('means').(type{t}).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}) = mean(emgMetrics.('stack').(type{t}).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}));
%                          % std, var, n measures
%                          emgMetrics.('stats').(type{t}).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}) = std(emgMetrics.('stack').(type{t}).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}));
%                          emgMetrics.('stats').(type{t}).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}) = [emgMetrics.('stats').(type{t}).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}) , var(emgMetrics.('stack').(type{t}).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}))];
%                          emgMetrics.('stats').(type{t}).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}) = size(emgMetrics.('stack').(type{t}).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}),1);
%                     else
%                          lOM = size(emgMetrics.('stack').(type{t}).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}),2);
%                          sOM = size(emgMetrics.('stack').(type{t}).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}),1);
%                          tempVar = zeros(1, lOM);
%                          tempVar1 = zeros(1, lOM);
%                          tempVar2 = zeros(1, lOM);
%                          for l = 1:length(lOM)
%                               tempVar(1,l) = mean(emgMetrics.('stack').(type{t}).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM})(:,l));
%                               tempVar1(1,l) = var(emgMetrics.('stack').(type{t}).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM})(:,l));
%                               tempVar2(1,l) = std(emgMetrics.('stack').(type{t}).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM})(:,l));
%                          end
%                          emgMetrics.('means').(type{t}).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}) = tempVar;
%                          emgMetrics.('stats').(type{t}).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}) = tempVar1;
%                          emgMetrics.('stats').(type{t}).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}) = [emgMetrics.('stats').(type{t}).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}) , tempVar2];
%                          emgMetrics.('stats').(type{t}).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}) = [emgMetrics.('stats').(type{t}).(gaitPhases{gP}).(emgMeasures{emgM}).(theMeasures{tM}) , sOM];
%                     end
%                end
%           end
%      end
% end
% 
% % Save it
% save([sessionData, filesep, 'emgMetrics.mat'], 'emgMetrics');
% 
% 
% end
