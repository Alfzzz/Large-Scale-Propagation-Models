clear
clc

q=5;
f1=2500;
%% Record your voice for q seconds.
recObj = audiorecorder;
disp('Start speaking.')
recordblocking(recObj, q);
disp('End of Recording.');

% Play back the recording.
play(recObj);

% Store data in double-precision array.
myRecording = getaudiodata(recObj);
% Time axis
qa=recObj.TotalSamples;
t=(0:q/qa:q-q/qa)';

%% for the frrequency axis
Ts=q/qa; %sampling time
fs=1/Ts; %sampling frequency
[na,nb]=size(t(:)); % na=number of points of signal
ff=fs*[0:na-1]/na-fs/2;

% Plot the waveform.
plot(t,myRecording);
xlabel('time (secs)')
ylabel('amplitude (V)')
figure

cs=cos(2*pi*f1*t);
mods=cs.*myRecording;
plot(ff,fftshift(abs(fft(myRecording))),'r')
% hold on
% plot(ff,fftshift(abs(fft(mods))),'k')
% hold off
xlabel('frequency (Hz)')
ylabel('Magnitude (V)')






