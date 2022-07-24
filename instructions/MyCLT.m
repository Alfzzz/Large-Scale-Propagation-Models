close all; clc; clear;

%% Program to see graphically the Central Limit Theorem based on uniform
% random variables.
% Generate uniformly distributed random numbers, and obtain the density
% or distribution using histogram and relative frequency of the sum of the
% random variables.

rng('shuffle'); % seed for random number generator

N=10000; % Number of random numbers to generate
M=50; % Number of random vectors to sum, M>=2 or pdf to convolve
nbins=50; % number of bins in histogram, must be < N
npts = 20; % Number of points of PDF
a = -5; % parameters of the uniform distribution
b = 5;

%% generate vector of random uniform numbers between a and b
x = a + (b-a).*rand(N,1);
% Get histogram of vector x
HH=histogram(x,nbins,'Normalization','pdf')
qq=zeros(M,1);
qq(1)=max(HH.Values);
hold on
for i=2:M
    % Generate new uniformly distributed vector and accumulate it
    x = x + (a + (b-a).*rand(N,1));
    % Get histogram of vector x
    hh=histogram(x,nbins,'Normalization','pdf')
    qq(i)=max(hh.Values);
end
hold off

%% PDF of uniformly distributed r.v. U[a,b]
xx=(a:(b-a)/npts:b)';
[sx,sy]=size(xx(:));
fx=(1/(b-a))*ones(sx,1);
fx2=fx;
figure
plot(xx,fx,'k')
a1=a;
b1=b;
sx1=sx;
qq2=zeros(M,1);
qq2(1)=max(fx2);
hold on
% Convolution of PDFs
for k=2:M
    % Accumulate convolutions
    fx2=conv(fx2,fx);
    % Generate x axis of length(fx2) + length(fx) - 1
    a1=a1+a;
    b1=b1+b;
    nk=length(fx2);
    xx2=(a1:(b1-a1)/(nk-1):b1)';
    qq2(k)=max(fx2);
    cft=qq(k)/qq2(k);
    plot(xx2,cft*fx2)
    sx1=nk;
end
hold off
grid



