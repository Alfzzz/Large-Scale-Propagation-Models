%% Evaluación 2
%Alfredo Zhu Chen
%% Clear environment
clc
clear all
close all
%% Read data
%********Adaptar los datos del .xlsx********************
dbm_sala=xlsread('potenciaWifi.xlsx',1,'C3:J13'); %datos del cuarto donde se encuentra la antena
dbm_cocina=xlsread('potenciaWifi.xlsx',1,'C15:J24'); %datos del cuarto adyacente a donde se encuentra la antena

%% Present heat maps of the power received of the WiFi signal in your house, for at least two
%rooms, one where the modem is and an adjacent room.
figure %figura 1
surface(dbm_sala); %mapa del cuarto donde se encuentra la antena
xticks(1:1:8)
xticklabels(0:0.5:3.5)
yticks(1:1:11)
yticklabels(0:0.5:5)
xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
ylabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
zlabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
title("Propagation map of the power received in the room with the modem")
colorbar
figure %figura 2
surface(dbm_cocina);%mapa del cuarto adayacente a donde se encuentra la antena
xticks(1:1:8)
xticklabels(0:0.5:3.5)
yticks(1:1:10)
yticklabels(6:0.5:10.5)
xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
ylabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
zlabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
title("Propagation map of the power received in the adjacent room")
colorbar

%% Present plots of the power received data measurements as distance increases, show plots of
%these data together with the linear regression, the free space and the two-ray models.
[X_sala,Y_sala]=meshgrid(0:0.5:3.5,0:0.5:5);   %Cordenadas
D_sala=sqrt((X_sala).^2+Y_sala.^2); %Distancias euclidianas del cuarto con modem
%dbm_sala(dbm==0)=inf; %poner la potencia en donde se encuentra la antena en infinito
dc_sala=(D_sala(end:-1:1)).'; %pasar las distancias euclidianas a un vector columna
dbmc_sala=(dbm_sala(end:-1:1)).'; %pasar las potencias recibidas en un vector columna
d_and_dbm_sala=[dc_sala(:) dbmc_sala(:)]; %concatenación de las distancias euclidianas con potencias recibidas
d_and_dbm_sala=d_and_dbm_sala(d_and_dbm_sala(:,1)>0,:); %Quitar potencias recibidas mayores a 0 del cuarto con modem
d_sala=d_and_dbm_sala(:,1); %distancias
dbmc_sala=d_and_dbm_sala(:,2); %potencias recibidas
d_log_sala=log10(d_sala);  %pasar a log10 las distancias

[X_cocina,Y_cocina]=meshgrid(0:0.5:3.5,6:0.5:10.5);  
D_cocina=sqrt((X_cocina).^2+Y_cocina.^2);
%dbm_cocina(dbm==0)=inf;
dc_cocina=(D_cocina(end:-1:1)).'; 
dbmc_cocina=(dbm_cocina(end:-1:1)).'; 
d_and_dbm_cocina=[dc_cocina(:) dbmc_cocina(:)];
d_and_dbm_cocina=d_and_dbm_cocina(d_and_dbm_cocina(:,2)>-80,:);%Quitar las potencias de donde se encuentra el refrigerador
d_cocina=d_and_dbm_cocina(:,1); 
dbmc_cocina=d_and_dbm_cocina(:,2); 
d_log_cocina=log10(d_cocina);  

% Free space Modelo de espacio libre
d0=1; %distancia pivote
PR_d0=-36.5; %potencia a distancia pivote d0=1m en dbm con LOS
n_free=2; %exponente del modelo de espacio libre
PR_dbm_free_sala=PR_d0+10*n_free*log10(d0./d_sala); %Potencia recibida con modelo de espacio libre en el cuarto con el modem
PR_dbm_free_cocina=PR_d0+10*n_free*log10(d0./d_cocina); %Potencia recibida con modelo de espacio libre en el cuarto adyacente

%2-ray model Modelo de 2 rayos
n_2R=4; %exponente del modelo de 2 rayos
PR_dbm_2R_sala=PR_d0+10*n_2R*log10(d0./d_sala); %potencia recibida con el modelo de 2 rayos del cuarto con modem
PR_dbm_2R_cocina=PR_d0+10*n_2R*log10(d0./d_cocina); %potencia recibida con el modelo de 2 rayos del cuarto adyacente

%Linear regression
coeficientes_sala=polyfit(d_log_sala,dbmc_sala,1);  %obtener coefficientes a partir de la regresión lineal
ple_sala=abs(coeficientes_sala(1,1))/10; %dividir la pendiente entre 10 para obtener el Path Loss Exponent
y_LR_sala=polyval(coeficientes_sala,d_log_sala); %evaluar las distancias con la regresión lineal que se obtuvo
PR_dbm_LR_sala=PR_d0+10*ple_sala*log10(d0./d_sala);
coeficientes_cocina=polyfit(d_log_cocina,dbmc_cocina,1);  
ple_cocina=abs(coeficientes_cocina(1,1))/10; 
y_LR_cocina=polyval(coeficientes_cocina,d_log_cocina);
PR_dbm_LR_cocina=PR_d0+10*ple_cocina*log10(d0./d_cocina);

figure %figura3 Sala
semilogx(d_sala,dbmc_sala,'*r',d_sala,y_LR_sala,'b',d_sala,PR_dbm_LR_sala,'c',d_sala,PR_dbm_2R_sala,'k',d_sala,PR_dbm_free_sala,'m');
xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
ylabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
title('Power received in the room with the modem')
legend('Data','Linear regression(LR)','ple in Model','2-ray model','Free model');

figure %figura4 Cocina
semilogx(d_cocina,dbmc_cocina,'*r',d_cocina,y_LR_cocina,'b',d_cocina,PR_dbm_LR_cocina,'c',d_cocina,PR_dbm_2R_cocina,'k',d_cocina,PR_dbm_free_cocina,'m');
xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
ylabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
title('Power received in adyacent room')
legend('Data','Linear regression(LR)','ple in Model','2-ray model','Free model');

%% Determine the PLE with the slope of the linear regression, and with the MMSE method,
%discuss their similarity or differences and determine which one works best for the area you
%chose by plotting the line corresponding to the simple model with these slopes.
n_ML_sala=sum(PR_d0*log10(d_sala)-log10(d_sala).*dbmc_sala)./sum((10*((log10(d_sala)).^2))); %exponente usando MMSE
PR_dbm_simple_sala=PR_d0+10*n_ML_sala*log10(d0./d_sala); %Potencias recibidas con modelo simple
figure %figura 5 Sala
semilogx(d_sala,dbmc_sala,'*r',d_sala,y_LR_sala,'b',d_sala,PR_dbm_LR_sala,'c',d_sala,PR_dbm_2R_sala,'k',d_sala,PR_dbm_free_sala,'m',d_sala,PR_dbm_simple_sala,'g');
xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
ylabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
title('Power received in the room with the modem')
legend('Data','Linear regression(LR)','ple in Model','2-ray model','Free model','Simple space');

n_ML_cocina=sum(PR_d0*log10(d_cocina)-log10(d_cocina).*dbmc_cocina)./sum((10*((log10(d_cocina)).^2))); %exponente usando MMSE
PR_dbm_simple_cocina=PR_d0+10*n_ML_cocina*log10(d0./d_cocina); %Potencias recibidas con modelo simple
figure %figura 5 Sala
semilogx(d_cocina,dbmc_cocina,'*r',d_cocina,y_LR_cocina,'b',d_cocina,PR_dbm_LR_cocina,'c',d_cocina,PR_dbm_2R_cocina,'k',d_cocina,PR_dbm_free_cocina,'m',d_cocina,PR_dbm_simple_cocina,'g');
xlabel('Distance $d$ (m)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
ylabel('$P_{R,dBm}(d)$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
title('Power received in the adjacent room')
legend('Data','Linear regression(LR)','ple in Model','2-ray model','Free model','Simple space');

%% Determine the standard deviation in dB ? and by taking the STD from the data. Discuss
%their similarities and which one works best.
sigmaML_sala=sqrt(sum((dbmc_sala-PR_d0+10*n_ML_sala*log10(d_sala)).^2)/size(d_sala,1));
sigmaML_cocina=sqrt(sum((dbmc_cocina-PR_d0+10*n_ML_cocina*log10(d_cocina)).^2)/size(d_cocina,1));
sigmaData_sala=std(dbmc_sala);
sigmaData_cocina=std(dbmc_cocina);