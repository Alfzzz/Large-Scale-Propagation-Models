%%
clear all
close all
clc
%
%% Transmission power
pt=1e-3;

% Generate axes
% Generate number of samples per axis (X,Y)
%[X,Y]=meshgrid(0:7,28:-1:0);
[X,Y]=meshgrid(0:7,5.5:-0.5:-8.5);
% Each sample was taken at 0.5 meter increments
X=X*.5;
%Y=Y*.5;
[i,j]=size(X);

% Antenna was located at coordinates X=0, Y=8.5 (X=0, Y=0)
% Get Euclidean distances from every point to the antenna and store these
% in an array
%D=sqrt((X).^2+(Y-8.5).^2);
D=sqrt((X).^2+Y.^2);
[i2,j2]=size(D);

% Read Excel file (file name, Sheet, Range)
dbm=xlsread('Book2.xlsx',1,'B2:I30');

% Standard Deviation obtained from data
% change the matrix to column vector
dbm2=dbm(:);
sigma=std(dbm2);
% Generate a matrix of normal (Gaussian) random numbers the same size as
% the matrix of distances. The mean value of the random numbers is zero
y=sigma*randn(i2,j2);

% For plotting purposes, change 0 dbm to infinity
dbm(dbm==0)=inf;
%pr=pt*(d.^(-ple)
dc=(D(end:-1:1)).'; %vector columna de la matrix Distancia 'D'
dbmc=(dbm(end:-1:1)).'; %vector columna de la matrix pr (potencia recibida)
[r, c]=size(dc);
% for n=r:-1:1
%     if dc(n,1)>0
%         Dc(n,1)=log10(dc(n,1));
%     end
% end

% Array of distances and dbm
cc=[dc(:) dbmc(:)];
% Generate matrix from cc including only those numbers > 0 in the first
% column of cc
z3=cc(cc(:,1)>0,:);
Dc=log10(z3(:,1));
dbmc=z3(:,2);

%% Linear regression with polynomial of order 1 
% polyfit(logarithmic distances,dbm,order) the result is the set of
% coefficients of polynomial are given in order starting with that for the
% highest power of x
coeficientes=polyfit(Dc,dbmc,1); 
% Evaluate polynomial with coefficients just obtained at locarithmic 
% distances Dc. This will be the line in the plot for the linear regression
yy=polyval(coeficientes,Dc);
% take first coefficient (slope of line) and divide by 10 to get PLE
ple=abs(coeficientes(1,1))/10;

% Here you need to calculate your PLE using maximum likelihood

% Here you need to calculate the received power using 2-ray model and
% free-space model

% Plot in logarithmic x axis the dbm and the straight line. For the dbm use
% blue stars, for the line use color red
semilogx(z3(:,1),dbmc,'*b',z3(:,1),yy,'r')
xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
ylabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
title('Power received, data vs. regression')
grid
% You will need to add on the same plot those for 2-ray and free-space
% propagation models


% Calculate received power in dBm using simple model pr=(D.^(-ple)), get
% these from the line equation in dBm
prdbm=-10*ple*log10(D)+coeficientes(2);

%% Plot surfaces for dbm
% first plot is for data measured
figure(2)
subplot(1,3,1)
surf(X,Y,dbm)
title('Power received, measurements (dBm)')
xlabel('meters')
ylabel('meters')
zlabel('$P_{R,dBm}(d)$','Interpreter','latex')
% shading 'interp'
% colormap jet
view(0,90)
colorbar

%% Second plot is for received power in dBm obtained from straight line
%figure(3)
subplot(1,3,2)
surf(X,Y,prdbm)
% hold on
% surf(X,Y,dbm)
% hold off
title('Power received, Simple model (dBm)')
xlabel('meters')
ylabel('meters')
zlabel('$P_{R,dBm}(d)$','Interpreter','latex')
% shading 'interp'
% colormap jet
view(0,90)
colorbar

%% Get dBm received by using dBm from straight line equation plus the normal
% (Gaussian) random variable y, already obtained in line 32 of this program
prdbm_r=prdbm+y;

% Plot surface for dBm predicted by model (line plus random number)
%figure(4)
subplot(1,3,3)
surf(X,Y,prdbm_r)
% hold on
% surf(X,Y,dbm)
% hold off
title('Power received, Model and rv (dBm)')
xlabel('meters')
ylabel('meters')
zlabel('$P_{R,dBm}(d)$','Interpreter','latex')
% shading 'interp'
% colormap jet
view(0,90)
colorbar
map = [0  0  0;
       0  0.5  0.25;
       0  0.5  0.5;
       0 0.5 0.75;
       0  0.5  1;
       0  1 1;
       0.5  1 1;
       1  1  0];
