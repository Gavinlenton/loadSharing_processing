function [plot1, plot2, plot3, plot4, plot5] = emgAnalysis(EMG, Fs, fcolow, fcohigh, actualFrames)
%Function to process raw EMG data and run an FFT on the data
%   Input raw EMG signal, the frequency it was collected at, the low-pass
%   filter frequency, the high-pass filter frequency, and vector with number of frames in trial. 

EMGdetrend = detrend(EMG);
Fnyq = Fs/2;

%% Plot the raw data and rectified signal
plot1 = figure(1);
subplot(2,1,1)
plot(EMGdetrend), xlabel(' (Hz)'), ylabel('Volts (V)'), title('EMG Raw')
subplot(2,1,2)
plot(EMGdetrend), xlabel(' (Hz)'), ylabel('Amplitude (V)'), title('Remove DC Offset EMG')

%High-pass filter
[b,a] = butter(2,fcohigh/Fnyq,'high');
filter_EMG = filtfilt(b,a,EMGdetrend);

%Rectify signal (Obtain absolute value)
rec_EMG = abs(filter_EMG-mean(filter_EMG)); 


plot2 = figure(2); %Create new figure with plot of rectified signal
plot(rec_EMG), xlabel(' (Hz)'), ylabel('Amplitude (V)'), title('Rectified EMG')
% axis([min(actualFramess), max(actualFramess), 0, 6e-3]);

%% Low-pass filter
[b,a] = butter(2,fcolow*1.25/Fnyq,'low');
filter_EMG1 = filtfilt(b,a,rec_EMG);

%%
plot3 = figure(3);
plot(actualFrames, EMGdetrend, 'b',actualFrames, rec_EMG, 'g', actualFrames, filter_EMG1, 'r'), xlabel(' (Hz)'), ylabel('Amplitude (V)'), title('Low Pass filtered EMG')
legend('Raw', 'Rectified', 'Linear envelope');

%% Vector of frequencies present in spectra
N = length(EMGdetrend);
freqs = 0:Fs/N:Fnyq;

%Compute fft and plot amplitude spectrum
LGastroc_fft = fft(EMGdetrend-mean(EMGdetrend));
plot4 = figure(4);
plot(freqs, abs(LGastroc_fft(1:N/2+1)));
xlabel('Sampling frequency up to Nyquist frequency (Hz)'), ylabel('Magnitude'), title(sprintf('FFT at %d Hz', Fs));
% axis([0, Fnyq, 0, 0.25])
%Compute and plot power spectrum
Px_Gastroc = LGastroc_fft.*conj(LGastroc_fft);
plot5 = figure(5);
plot(freqs,abs(Px_Gastroc(1:N/2+1)));
xlabel('Sampling frequency up to Nyquist frequency (Hz)'), ylabel('Magnitude'), title(sprintf('FFT at %d Hz', Fs));
% axis([0, Fnyq, 0, 0.07])

saveas(figure(3), sprintf('Envelope_%d', Fs), 'fig'); 
saveas(figure(4), sprintf('FFT_%d', Fs), 'fig'); 
saveas(figure(5), sprintf('powerSpec_%d', Fs), 'fig');
end

