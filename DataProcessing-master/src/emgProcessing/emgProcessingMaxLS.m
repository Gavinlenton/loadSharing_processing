function [] = emgProcessingMaxLS(isNotch, sessionData,  maxc3dFileName, motoDir)
% Input if Notch is required (Yes, No), session Data folder, name of max
% file name, and MOtoNMS directory
%   The Notch filter will only be applied for data collected directly
%   through Nexus. The sessionData is from the MOtoNMS c3d2mat function

%% EMG scaling code
% Designed to work with Load sharing data set
% Code processes the user-provided max trial
% Code then stores all max emg values to file using CEINMS standard format.

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

%% Directories and files
motoNmsSrcDir = motoDir;
% trialsFile = 'C:\Users\s2790936\Desktop\SCOPEX\SCOPEX Trials DJS.xlsm';
% maxDir = 'C:\Users\s2790936\Desktop\SCOPEXTemp\Baseline'; % Max's

originalPath=pwd;
cd([motoNmsSrcDir filesep, 'src']);
cd('shared')
sharedFunctionsPath=pwd;
addpath(sharedFunctionsPath)
cd(originalPath)

%% Load session data
maxData = [sessionData, filesep, maxc3dFileName]; % Max trial dir


%% Specify emg.mat file
suffix = 'AnalogData.mat';

% Specify for max trial
filenameMax = fullfile(maxData,  suffix);

%% Load emgMax mat file

% File 1
try
     emgData1 = load(filenameMax);
catch me1
     disp(['No such file ', filenameMax, ' , continuing with analysis'])
end

%% Initialize some values
try
     channels = fieldnames(emgData1);
catch me
     disp(['No such file ', filenameMax, ' , check the filename'])
end
emgMax = zeros(2, length(emgLabelsInMatFile));

%% Condition raw emg signal
     nDataPoints = 1;
     for xx = 1:length(channels)
          lengthOfData = size(emgData1.(channels{xx}).RawData, 1);
          if lengthOfData > nDataPoints
               nDataPoints = lengthOfData;
          end
     end

%% Initialize matrix
emgData = zeros(nDataPoints, size(emgLabelsInMatFile,2));


for n = 1:length(emgLabelsInMatFile)
     for x = 1:size(channels,1)
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
     end
end


% Load max data analog file
analogFile = fullfile([sessionData, filesep, analogData, filesep, 'AnalogData.mat']);
try
     load(analogFile);
catch me
     error(['Cannot load analog data mat file for ' , analogData, ' ensure this data has been converted from c3d to mat'])
end

%% Find max signal 
     
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
          
          % Trial from which max was from
          
          MaxEMG_trials{l} = filenameMax;
               
          
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
     
     % Print max emg to MotoNMS sessionData directory.
     printDir = [maxData, filesep, 'maxEMG'];
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
end
