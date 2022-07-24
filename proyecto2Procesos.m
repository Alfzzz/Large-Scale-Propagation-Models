%%
%Alfredo Zhu Chen
%Andrés Islas Bravo
%Luis Arturo Dan Fong
%Juan Pablo Valverde López

%% Clear environment
clc
clear all
close all
%% Read data
dbm=xlsread('potenciaWifi.xlsx',1,'C3:J13');

%%  Provide a dBm received power 3D map (surface) detected, use interpolation to generate received 
%power at more points than those measured. 
%********Adaptar la meshgrid de acuerdo a las coordenadas********************
[X,Y]=meshgrid(0:0.5:3.5,0:0.5:5)   %Coordinates
[XInt,YInt]=meshgrid(0:0.1:3.5,0:0.1:5) %Coordinates interpolated

var_new = griddata(X,Y,dbm,XInt,YInt,'cubic');
figure
surface(var_new); %Interpolated 3D map

%% For the antenna, with your data obtained, provide a 2D plot of distance in meters vs. power 
%received in dBm. Consider the average of dBms when you have several measurements at the same 
%distances. This plot is distance dependent, not location dependent.
[i,j]=size(X);
%********Adaptar la ecuación a la posición del modem********************
D=sqrt((X).^2+Y.^2);
[i2,j2]=size(D);
dbm2=dbm(:);
sigma=std(dbm2);

dbm(dbm==0)=inf;
dc=(D(end:-1:1)).'; %vector columna de la matrix Distancia 'D'
dbmc=(dbm(end:-1:1)).'; %vector columna de la matrix pr (potencia recibida)
[r, c]=size(dc);

% Array of distances and dbm
[dcSorted,I] =sort(dc);
dbmcSorted = dbmc(I);
cc=[dcSorted(:) dbmcSorted(:)]; %concatenación de las distancias euclidianas con potencias recibidas
% Generate matrix from cc including only those numbers > 0 in the first
% column of cc
z3=cc(cc(:,1)>0,:); %Quitar potencias recibidas mayores a 0
d=z3(:,1) %distance
Dc=log10(d);  %pasar a log10 las distancias
dbmc=z3(:,2); %Tomar datos de la segunda columna, potencias recibidas

count=0;
prom=[];
z3New=[];
for i=2:size(z3,1)
    check=z3(i-1,1)
    if z3(i,1)==check
        count=count+1 
    else
        for j=i-1:-1:i-count-1
            prom(j)=z3(j,2)
        end
        prom=mean(nonzeros(prom))
        z3New=[z3New;[z3(i-1,1),prom]]
        prom=[]
        count=0;
    end
end
figure
semilogx(z3New(:,1),z3New(:,2),'*b')
xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
ylabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
title('Power received, data vs. regression')
%% Use the free-space Friis model to estimate power received and path loss and plot them as a function of distance.
% Free space
%Model de espacio libre
%********Adaptar Gt,Gr,Pt,f de acuerdo a los datos********************
Gt=1;%0dB
Gr=1;%0dB
Pt=.001;%0dBm
f=2.4*10^9;%2.4GHertz
lambda=299792458/f;
PRel0=Gt*Gr*Pt*lambda^2/((4*pi)^2);
PRel=PRel0.*(1./(d.^2));%watts
PR_dbm_free=10*log10(PRel.*1000);%conversion a dbms
%figure
%semilogx(d,PR_dbm_free,'g','LineWidth',2);
%xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
%ylabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
%title('Power received, data vs. regression')
%%  2-ray model
%********Adaptar las alturas de acuerdo al caso********************
ht=1.1; %Tx high in m
hr=0.75; %Rx high in m
K = Pt*Gt*Gr*((ht*hr)^2);
PR2r=(K./(d.^4));%watts
PR_dbm_2R=10*log10(PR2r.*1000);%dbms
%semilogx(d,PR_dbm_2R,'k','LineWidth',2);
%xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
%ylabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
%title('Power received, data vs. regression')

%% Linear regression with polynomial of order 1 
% polyfit(logarithmic distances,dbm,order) the result is the set of
% coefficients of polynomial are given in order starting with that for the
% highest power of x
coeficientes=polyfit(Dc,dbmc,1); 
% Evaluate polynomial with coefficients just obtained at locarithmic 
% distances Dc. This will be the line in the plot for the linear regression
PR_dbm_LR=polyval(coeficientes,Dc);
% take first coefficient (slope of line) and divide by 10 to get PLE
ple=abs(coeficientes(1,1))/10;

% Here you need to calculate your PLE using maximum likelihood

% Here you need to calculate the received power using 2-ray model and
% free-space model

% Plot in logarithmic x axis the dbm and the straight line. For the dbm use
% blue stars, for the line use color red

%semilogx(d,dbmc,'*b',d,PR_dbm_linearRegression,'r')
%xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
%ylabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
%title('Power received, data vs. regression')
%grid

%% Simple model
%********Adaptar prdo de acuerdo al caso********************
prdo=-36.5; % in dBm
MLnum = ones(1,size(z3,1));
MLden = ones(1,size(z3,1));
for N = 1:size(z3,1)    
     %MLnum(1,N)=(log10(N/1)*(PromMed(N) - PR2rdbm(N)));
     MLnum(1,N)=prdo*   log10(z3(N,1))  -   log10(z3(N,1)) *  z3(N,2);
     MLden(1,N)=(10*((log10(z3(N,1)))^2));
end
num = 0;
den = 0;
for N2 = 1:size(z3,1)    
     num = num + MLnum(1,N2);
     den = den + MLden(1,N2);
end
nML = num./den;    %Exponente de perdidas n

%Modelo sencillo with Pr(do) and ML
do=1;
%nML=-polinomio(1)/10;
PR=Pt.*(d.^(-nML));%watts
%PRdbm=-10*log10(PR.*1000);%dbms
%PRdbm=10*log10(PR)+30;
PR_dbm_simple=prdo-10*nML*log10(d);

%semilogx(d,PR_dbm_simple,'m')
xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
ylabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
title('Power received, data vs. regression')
grid
figure
semilogx(d,dbmc,'*r',d,PR_dbm_LR,'b',d,PR_dbm_2R,'k',d,PR_dbm_free,'m',d,PR_dbm_simple,'g');
grid
xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
ylabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
title('Power received, data vs. regression')
legend('Data','Linear regression','2-ray model','Free model','Simple space');
%% Propagation map of the power received in the rooms by coloring the area according 
%to the power received measured, using the model that achieved the best fit
%****************Adaptar x,y para distancias**************************** 

map_dbm=ones(size(dbm))
for x=0:0.5:3.5
    for y=0:0.5:5
        %****************Uncommnet the map for the best model*****
        %map_dbm(1+y*2,1+x*2)=10*log10((K./(sqrt(x^2. +y^2.).^4)).*1000);%2-Ray model
        map_dbm(1+y*2,1+x*2)=polyval(coeficientes,sqrt(x^2. +y^2.));%Linear Regression
    end
end  
figure
surf(map_dbm)

%% estimate the ? of the environment using statistical method ML
sigmaMLArr = ones(1,size(d,1));
for N = 1:size(d,1)    
     sigmaMLArr(1,N)=(10*nML*log10(d(N,1))+(z3(N,2) - PR_dbm_2R(N,1))).^2;
end
sigmaMLArr2 =0;
for N2 = 1:size(d,1)     
     sigmaMLArr2 = sigmaMLArr2 + sigmaMLArr(1,N2);
end
sigmaML = (1/size(d,1))*sqrt(sigmaMLArr2);

sigma2Ray = std(PR_dbm_2R);
sigmaLR = std(PR_dbm_LR);

%% estimate the standard deviation of your data by using statistical tools (averages)
sigmaMeds = std(z3(:,2));

%% superimpose the result of adding such best model and the Gaussian random numbers (in dB) that 
%simulate fading effects. You will add to each point a Gaussian random variable with zero mean 
%and standard deviation ?
%******************Use this if the 2-Ray model is the best one*************
%{
PR_dbm_2R_RV=PR_dbm_2R+normrnd(0,sigma2Ray,size(PR_dbm_2R,1),1);
figure
hold on
semilogx(d,PR_dbm_2R_RV,'k')
semilogx(d, dbmc,'r*')
title('Graficas de mediciones y modelos de Potencias recibidas + ruido')
ylabel('dBm')
xlabel('distancia')
legend('Modelo con variable aleatoria','Mediciones');
hold off
%}
%******************Use this if the Linear Regression is the best one*************
PR_dbm_LR_RV=PR_dbm_LR+normrnd(0,sigmaLR,size(PR_dbm_LR,1),1);
figure
hold on
semilogx(d,PR_dbm_LR_RV,'k')
semilogx(d, dbmc,'r*')
title('Graficas de mediciones y modelos de Potencias recibidas + ruido')
ylabel('dBm')
xlabel('distancia')
legend('Modelo con variable aleatoria','Mediciones');
hold off

%% 3D map given the meassurments and 3D map of the best model with RV
figure
surface(dbm); %Non interpolated 3D map
map_dbm_RV=map_dbm+normrnd(0,sigmaLR,size(map_dbm,1),size(map_dbm,2));
figure
surface(map_dbm_RV)

%% Estimate coverage by calculating the outage probability with your best model
gamma_db=-37;
%*****Uncommnent the line for the best model****
outage=1-qfunc((gamma_db-mean(PR_dbm_LR))/std(PR_dbm_LR)); %Linear Regression
%outage=1-qfunc((gamma_db-mean(PR_dbm_2R))/std(PR_dbm_2R)); %2-Ray model
outage_idx=find(dbm<gamma_db)
dbm_outageCol=dbm;
dbm_outageCol(outage_idx)=-37;
outage_idx=find(dbm>gamma_db)
dbm_outageBlack=dbm;
dbm_outageBlack(outage_idx)=-37;

figure
hold on
s1=surface(dbm_outageCol)
%colormap hot
s1.FaceColor='blue';
hold on
%caxis('manual')
s2=surface(dbm_outageBlack)
%colormap([0 0 0]);
s2.FaceColor='black';


