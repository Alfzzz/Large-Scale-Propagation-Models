clear all;
close all;
clc

%% Transmission power or P(do)
pt=1e-3;

% Read Excel file (file name, Sheet, Range)
dbm=xlsread('Book2.xlsx',3,'B2:I30');

%% Generate axes
% Generate number of samples per axis (X,Y)
%[X,Y]=meshgrid(0:7,28:-1:0);
[X,Y]=meshgrid(0:7,5.5:-0.5:-8.5);
% Each sample was taken at 0.5 meters increments
X=X*.5;
%Y=Y*.5;
[i,j]=size(X);

%% Antenna was located at coordinates X=0, Y=0
% Get Euclidean distances from every point to the antenna and store these
% in an array
D=sqrt((X).^2+Y.^2);
[i2,j2]=size(D);

%% Get from dBm matrix only those elements that are less than zero, the
% result will be given in a column vector
mat1=dbm(dbm<0);
% For those elements with dbm less than zero, get those distances in a
% column vector
mat2=D(dbm<0);

%% Get PLE using linear regression
% Linear regression with polynomial of order 1 
% polyfit(logarithmic distances,dbm,order) the result is the set of
% coefficients for the polynomial, which are given in order starting with 
% that coefficient for the highest power of x
coefficient=polyfit(log10(mat2),mat1,1); 

% Evaluate polynomial with coefficients just obtained at locarithmic 
% distances Dc. This will be the line in the plot for the linear regression
yy=polyval(coefficient,log10(mat2));
% take first coefficient (slope of line) and divide by 10 to get PLE
ple=abs(coefficient(1,1))/10;

%% Get standard deviation from data
sigma=std(mat1);
% Generate a matrix of normal (Gaussian) random numbers the same size as
% the matrix of distances. The mean value of the random numbers is zero
y=sigma*randn(i2,j2);

% Use of simple model. You could also use the equation with do, i.e., get
% the value using pt with P(do), recall that do=1m
dbm2=10*log10(pt*(D.^(-ple)));
for i=1:i2
    for j=1:j2
        if (dbm(i,j) == 0)
            dbm3(i,j)=dbm2(i,j)+y(i,j); %Matrix with random dBm
        else
            dbm3(i,j)=dbm(i,j);
        end
    end
end

% dbm3=dbm;
% dbm3(dbm3==0)=dbm2(dbm3==0)+y(dbm3==0);

semilogx(mat2,mat1,'*r',mat2,yy,'b')
xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
ylabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
title('Power received, data vs. regression')
grid



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

%% second plot is for simple model, no random variable
%figure(3)
subplot(1,3,2)
surf(X,Y,dbm2)
title('Power received, Simple model (dBm)')
xlabel('meters')
ylabel('meters')
zlabel('$P_{R,dBm}(d)$','Interpreter','latex')
% shading 'interp'
% colormap jet
view(0,90)
colorbar

%% third plot is for model plus random variable
%figure(4)
subplot(1,3,3)
surf(X,Y,dbm3)
title('Power received, Model and rv (dBm)')
xlabel('meters')
ylabel('meters')
zlabel('$P_{R,dBm}(d)$','Interpreter','latex')
% shading 'interp'
% colormap jet
view(0,90)
colorbar
