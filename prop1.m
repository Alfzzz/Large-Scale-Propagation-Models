%%
clear all
close all
clc

pt=0.5;
ple=3.5;
a=8; % Standard deviation of Gaussian rv

[X,Y] = meshgrid(1:100,100:-1:1);
X=X*0.5;
Y=Y*0.5;
[i,j]=size(X);

% Antenna is assumed to be at position X=7.5, Y=7.5
D=sqrt((X-7.5).^2+(Y-7.5).^2);

% Simple model
pr=pt*(D.^(-ple));

prdbm=10*log10(pr)+30;

surf(X,Y,prdbm)
title('Power received, measurements (dBm)')
xlabel('meters')
ylabel('meters')
zlabel('$P_{R,dBm}(d)$','Interpreter','latex')
view(0,90)

% Gaussian rv
y = a.*randn(i,j);

% Model plus rv
prdbm_r=prdbm+y;
% figure
hold on
% Plot random surface 100 dBm above smooth surface
surf(X,Y,prdbm_r+100)
title('Power received, measurements plus rv (dBm)')
xlabel('meters')
ylabel('meters')
zlabel('$P_{R,dBm}(d)$','Interpreter','latex')
view(0,90)
hold off
% xx=xlsread('test2.xlsx')
