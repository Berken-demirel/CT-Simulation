%% Berken Utku Demirel - 2166221
clearvars
close all
clc
%% Q1)
% Load the dataset
data = load('SingleLeadMeasurement.mat').SingleLeadMeasurement;
fs = 360;
t = 1/fs:1/fs:10;
% Plot the waveform
figure,
plot(t,data)
title('The plot of waveform')
xlabel('Time in seconds')
% Find the index of the peaks 
[~,peaks,delay] = pan_tompkin(data,fs,0);
peaks_index = find_peaks(data,peaks,fs);
% Find the beats per minute
[RR_interval,bpm] = calculate_RR_intervals(peaks_index,fs);

% Check the correct peak
figure,
plot(t,data)
hold on
plot(t(peaks),data(peaks),'*r')
title('ECG waveform with R-peaks')
xlabel('Time in seconds')
% Plot the beats per minute vs time
time = linspace(0,10,12);
figure,
plot(time,bpm)
title('Beats per minute')
xlabel('Time in seconds')
ylabel('Beats per minute (bpm)')
%% Q2)
data2 = load('2LeadMeasurements.mat');
fs = 257;
data2_lead1 = data2.Lead_I;
data2_lead2 = data2.Lead_II;
% Plot the Lead I signal in the specified interval
interval_begin = 23 * 60 * fs;
interval_end = interval_begin + fs * 5;

figure,
plot(data2_lead1(1,interval_begin:interval_end))
title('Plot of Lead1 for the specified interval')
% CTFT of the overall signal
T = 1/fs;
fft_signal = T * 0.25 * abs(fftshift(fft(data2_lead1)));
frequency_axis = linspace(-fs/2,fs/2,length(fft_signal));
figure,
plot(frequency_axis,fft_signal)
title('The CTFT of the overall signal')
xlabel('The frequency in Hertz')
ylabel('Magnitude')

%% Q3)
L = 2;
h_n = 1/L * ones(1,L);
filtered_signal = conv(h_n,data2_lead1);
filtered_signal = filtered_signal(:,1+L/2:end-L/2+1);

baseline = data2_lead1 - filtered_signal;
fft_signal = T * 0.25 * abs(fftshift(fft(filtered_signal)));

figure,
plot(data2_lead1(1,interval_begin:interval_end))
hold on 
plot(filtered_signal(:,interval_begin:interval_end))
hold on
plot(baseline(:,interval_begin:interval_end))
legend('Raw data','Filtered data','Baseline')
title('Filtered signal')

figure,
plot(frequency_axis,fft_signal)
title('The CTFT of the baseline-corrected signal')
xlabel('The frequency in Hertz')
ylabel('Magnitude')

%% Q4)
wo = 50/(fs/2);  
bw = wo/5;
[b,a] = iirnotch(wo,bw);
%fvtool(notchfilt,'Color','white');
filtered_signal = filter(b,a,filtered_signal);

figure,
plot(filtered_signal(:,interval_begin:interval_end))
title('Filtered signal (Elimination of Powerline)')

T = 1/fs;
fft_signal = T * 0.25 * abs(fftshift(fft(filtered_signal)));
frequency_axis = linspace(-fs/2,fs/2,length(fft_signal));
figure,
plot(frequency_axis,fft_signal)
title('Frequecy response of the signal')
%% Q5)
filtered_signal_2 = conv(h_n,data2_lead2);
filtered_signal_2 = filtered_signal_2(:,1+L/2:end-L/2+1);
filtered_signal_2 = filter(b,a,filtered_signal_2);

Lead_III = filtered_signal_2 - filtered_signal;

figure,
plot(Lead_III(:,interval_begin:interval_end))
title('Lead III')
%% Q6)



