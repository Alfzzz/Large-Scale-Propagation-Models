close all; clc; clear;

%% Signal plus noise

% Limits of time axis
ti=0; % Initial time
tf=4; % final time

% Generate a cosine signal
A=2; % Amplitude in Volts
da=A/2; % for plotting purposes
fo=2; % Frequency in Hz
phi=pi/6; % phase of signal
delta=1/(100*fo); % Sampling period
t=(ti:delta:tf)'; % Time axis
[tx,ty]=size(t(:)); % size of time vector or number of samples
% Signal
xt=A*cos(2*pi*fo*t+phi);

% % Plot noiseless signal
% plot(t,xt,'b','LineWidth',2)
% grid
% xlabel('Time (sec)')
% ylabel('Amplitude (V)')
% axis([ti tf -A-da A+da])

%% Generate Gaussian random variables with unit variance and zero mean
rn = randn(tx,1); % Generate vector of random numbers

% Generate signal plus noise
xtn = xt + rn;
da2=max(rn);
da3=min(rn);

%figure
% Plot noisy signal
plot(t,xtn,'b')
hold on
plot(t,xt,'r','LineWidth',2)
hold off
grid
xlabel('Time (sec)')
ylabel('Amplitude (V)')
axis([ti tf -A+da3 A+da2])

% Statistics Noise
mn1=mean(rn);
stdn1=std(rn);

%% Generate Gaussian random variables with sigma variance and zero mean
sigma = 2;
rn2 = sigma*randn(tx,1); % Generate vector of random numbers

% Generate signal plus noise
xtn2 = xt + rn2;
da2s=max(rn2);
da3s=min(rn2);

figure
% Plot noisy signal
plot(t,xtn2,'b')
hold on
plot(t,xt,'r','LineWidth',2)
hold off
grid
xlabel('Time (sec)')
ylabel('Amplitude (V)')
axis([ti tf -A+da3s A+da2s])

% Statistics Noise
mn2=mean(rn2);
stdn2=std(rn2);



%% SNR
pws = A^2/2;
pwn1 = stdn1^2;
pwn2 = stdn2^2;

SNR1 = pws/pwn1;
SNR1_dB = 10*log10(SNR1);
SNR2 = pws/pwn2;
SNR2_dB = 10*log10(SNR2);


% Statistics
[mn1 mn2 ; pwn1 pwn2 ; SNR1 SNR2 ; SNR1_dB SNR2_dB]







