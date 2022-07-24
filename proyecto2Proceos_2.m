%%
%Alfredo Zhu Chen
%Andrés Islas Bravo
%Luis Arturo Dan Fong
%Juan Pablo Valverde López
%%%%%%%%%%%%Modelos utilizando Pr(d0)%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear environment
clc
clear all
close all
%% Read data
%********Adaptar los datos del .xlsx********************
dbm=xlsread('potenciaWifi.xlsx',1,'C3:J13');

%%  Provide a dBm received power 3D map (surface) detected, use interpolation to generate received 
%power at more points than those measured. 
%********Adaptar la meshgrid de acuerdo a las cordenadas********************
[X,Y]=meshgrid(0:0.5:3.5,0:0.5:5);   %Cordenadas
[XInt,YInt]=meshgrid(0:0.1:3.5,0:0.1:5); %Coordenadas con interpolación

mapInt=griddata(X,Y,dbm,XInt,YInt,'cubic'); %Interpolation
figure %figure 1
surface(mapInt); %mostrar mapa 3D con interpolación
%********Adaptar escalas de X y Y********************
xticks(0:5:35)
xticklabels(0:0.5:3.5)
yticks(0:5:50)
yticklabels(0:0.5:5)
xlim([1 35])
ylim([1 50])
xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
ylabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
zlabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
title("Propagation 3D map of the power received using interpolation")
colorbar
%% For the antenna, with your data obtained, provide a 2D plot of distance in meters vs. power 
%received in dBm. Consider the average of dBms when you have several measurements at the same 
%distances. This plot is distance dependent, not location dependent.
%********Adaptar la ecuación a la posición del modem********************
D=sqrt((X).^2+Y.^2); %Distancias euclidianas
sigma=std(dbm(:)); %Desviación estándard de las potencias recibidas
%dbm(dbm==0)=inf; %poner la potencia en donde se encuentra la antena en infinito
dc=(D(end:-1:1)).'; %pasar las distancias euclidianas a un vector columna
dbmc=(dbm(end:-1:1)).'; %pasar las potencias recibidas en un vector columna
[dcSorted,I] =sort(dc); %Ordenar las distancias 
dbmcSorted = dbmc(I); %ordenar las potencias recibidas
d_and_dbm=[dcSorted(:) dbmcSorted(:)]; %concatenación de las distancias euclidianas con potencias recibidas
d_and_dbm=d_and_dbm(d_and_dbm(:,1)>0,:); %Quitar potencias recibidas mayores a 0
d_and_dbm=d_and_dbm(d_and_dbm(:,2)>-90,:); %Quitar potencias donde no se pudieron hacer mediciones
d=d_and_dbm(:,1); %distancias
dbmc=d_and_dbm(:,2); %potencias recibidas
d_log=log10(d);  %pasar a log10 las distancias
count=0; %Contador de distancias repetidas
prom=[]; %promedios
d_and_dbm_meanD=[]; %promediar potencias que tienen distancias iguales
for i=2:size(d_and_dbm,1) %Interar todos los datos empezando por el segundo
    if d_and_dbm(i,1)==d_and_dbm(i-1,1) %Checar si la distancia actual es igual a la previa
        count=count+1; %Acumular si hay distancia repetida 
    else
        for j=i-1:-1:i-count-1 %Iterar con las distancias que se repiten
            prom=[prom d_and_dbm(j,2)]; %guardar(concatenar) potencias con mismas distancias
        end
        prom=mean(prom); %obtener el promedio de las potencias con las mismas distancias
        d_and_dbm_meanD=[d_and_dbm_meanD;[d_and_dbm(i-1,1),prom]]; %Concatenar la distancia y el promedio con matriz
        prom=[]; 
        count=0; %Resetear contador de distancias repetidas
    end
end
figure %figura 2
semilogx(d_and_dbm_meanD(:,1),d_and_dbm_meanD(:,2),'*b') %plot de datos, las potencias se promedian si tienen distancias iguales
xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
ylabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
title('Power received')
grid on
%% Use the free-space Friis model to estimate power received and path loss and plot them as a function of distance.
% Free space Modelo de espacio libre
%*******************Adaptar d0 y PR_d0 de acuerdo al caso*******************
d0=1; %distancia pivote
PR_d0=-36.5; %potencia a distancia pivote d0=1m en dbm con LOS
n_free=2; %exponente del modelo de espacio libre
PR_dbm_free=PR_d0+10*n_free*log10(d0./d); %Potencia recibida con modelo de espacio libre
%%  2-ray model Modelo de 2 rayos
n_2R=4; %exponente del modelo de 2 rayos
PR_dbm_2R=PR_d0+10*n_2R*log10(d0./d); %potencia recibida con el modelo de 2 rayos
%% Linear regression with polynomial of order 1 
% polyfit(logarithmic distances,dbm,order) the result is the set of 
%coefficients of polynomial are given in order starting with that for the
% highest power of x
coeficientes=polyfit(d_log,dbmc,1);  %obtener coefficientes a partir de la regresión lineal
ple=abs(coeficientes(1,1))/10; %dividir la pendiente entre 10 para obtener el Path Loss Exponent
y_LR=polyval(coeficientes,d_log); %evaluar las distancias con la regresión lineal que se obtuvo
PR_dbm_LR=PR_d0+10*ple*log10(d0./d); %evaluar ple en el modelo general
%% Now estimate the path loss exponent using the statistical method of maximum likelihood (ML)
%Simple model Modelo simple
n_ML=sum(PR_d0*log10(d)-log10(d).*dbmc)./sum((10*((log10(d)).^2))); %exponente del modelo simple 
PR_dbm_simple=PR_d0+10*n_ML*log10(d0./d); %Potencias recibidas con modelo simple
figure %figura3
semilogx(d,dbmc,'*r',d,y_LR,'b',d,PR_dbm_LR,'y',d,PR_dbm_2R,'k',d,PR_dbm_free,'m',d,PR_dbm_simple,'g');
xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
ylabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
title('Power received')
legend('Data','Linear regression(LR)','ple in Model','2-ray model','Free model','Simple model');
%% Squared errors for all models
n_all=[ple,ple,n_2R,n_free,n_ML]; %exponentes: regresión lineal, 2 rayos, espacio libre, simple
sq_errorLR1= mean((dbmc-polyval(coeficientes,d_log)).^2);
sq_errorLR2= mean((dbmc-PR_d0+10*n_all(2)*log10(d)).^2);
sq_error2R= mean((dbmc-PR_d0+10*n_all(3)*log10(d)).^2);
sq_errorFree= mean((dbmc-PR_d0+10*n_all(4)*log10(d)).^2);
sq_errorML= mean((dbmc-PR_d0+10*n_all(5)*log10(d)).^2);
sq_errors=[sq_errorLR1,sq_errorLR2,sq_error2R,sq_errorFree,sq_errorML];
[errormin,minError_idx]=min(sq_errors); %obtener el error menor y el índice
n_best=n_all(minError_idx); %exponente para el modelo con menor error
%% Propagation map of the power received in the rooms by coloring the area according 
%to the power received measured, using the model that achieved the best fit
%****************Adaptar x,y junto con el exponente del modelo mejor**************************** 
map_dbm_best=ones(size(dbm));
for x=0:0.5:3.5
    for y=0:0.5:5
        if n_best==n_all(1)
            map_dbm_best(1+y*2,1+x*2)=polyval(coeficientes,sqrt(x^2. +y^2.));
        else
            map_dbm_best(1+y*2,1+x*2)=PR_d0+10*n_best*log10(d0/sqrt(x^2. +y^2.));
        end
    end
end  
figure %figura 4
surface(map_dbm_best)
%*****************Adaptar a los datos!************************
xticks(1:1:8)
xticklabels(0:0.5:3.5)
yticks(1:1:11)
yticklabels(0:0.5:5)
xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
ylabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
zlabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
title("Propagation 3D map of the power received with the best model")
colorbar
%% estimate the sigma of the environment using statistical method ML
sigmaML=sqrt(sum((dbmc-PR_d0+10*n_ML*log10(d)).^2)/size(d,1));
sigma2R = std(PR_dbm_2R);
sigmaLR = std(PR_dbm_LR);
sigma_yLR=std(y_LR);
sigmaFree=std(PR_dbm_free);
%% estimate the standard deviation of your data by using statistical tools (averages)
sigmaMeds = std(dbmc);
%% superimpose the result of adding such best model and the Gaussian random numbers (in dB) that 
%simulate fading effects. You will add to each point a Gaussian random variable with zero mean 
%and standard deviation
sigma_all=[sigmaLR,sigma_yLR,sigma2R,sigmaFree,sigmaML];
sigma_best=sigma_all(minError_idx); %sigmal del mejor modelo
if n_best==ple
    PR_dbm_best=PR_dbm_LR; %potencia recibida de acuerdo al mejor modelo
else
    PR_dbm_best=PR_d0+10*n_best*log10(d0./d); %potencia recibida de acuerdo al mejor modelo
end
PR_dbm_best_RV=PR_dbm_best+normrnd(0,sigmaML,size(PR_dbm_best,1),1); %agregar RV al que tiene mejor modelo
figure %figura 5
hold on
semilogx(d,PR_dbm_best_RV,'k')
semilogx(d, dbmc,'r*')
title('Data and power received+noise')
ylabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
legend('Best model with random variable','Data');
hold off

%% 3D map given the meassurments and 3D map of the best model with RV
figure %figura 6
surface(dbm); %Non interpolated 3D map
%************Adaptar a los datos!***********************
xticks(1:1:8)
xticklabels(0:0.5:3.5)
yticks(1:1:11)
yticklabels(0:0.5:5)
xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
ylabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
zlabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
title("Propagation 3D map of the power received from the measurements")
colorbar
map_dbm_RV=map_dbm_best+lognrnd(0,sigma_best,size(map_dbm_best,1),size(map_dbm_best,2)); %Agregar RV al mejor modelo
figure %figura 7
surface(map_dbm_RV)
%************Adaptar a los datos!***********************
xticks(1:1:8)
xticklabels(0:0.5:3.5)
yticks(1:1:11)
yticklabels(0:0.5:5)
xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
ylabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
zlabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
title("Propagation 3D map of the power received with the best model with Random Variable")
colorbar
%% Estimate coverage by calculating the outage probability with your best model
%****************Cambiar gamma_db****************************
gamma_db=-40; %Valor umbral, solo los que son mayores a este cumplen con el criterio de outage
outage=1-qfunc((gamma_db-mean(PR_dbm_best))/std(PR_dbm_best));%probabilidad del outage
figure %figura 8
hold on
surface(dbm)
colormap default
cmap = colormap;       %current map
maximumdBm=max(dbm);
minimumdBm=min(dbm);
outColor=(gamma_db-minimumdBm)/(maximumdBm-minimumdBm)
cmap(1:floor(size(colormap,1)*outColor),:) =zeros(floor(size(colormap,1)*outColor),3);   %make first color white
colormap(cmap);
colorbar
%************Adaptar a los datos!***********************
xticks(1:1:8)
xticklabels(0:0.5:3.5)
yticks(1:1:11)
yticklabels(0:0.5:5)
xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
ylabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
zlabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
title("Outage 3D map ")
