
%% Proyecto 2: propagación
%Redes de comuniacion Personal
%Octubre 2014

clear all;
close all;
clc

rng('shuffle'); % seed for random number generator

%% Mapa de medi;thiciones
first = [-44 -40 -48 -48 -59 -52 -57 -61 -63 -62 -71 -61 -58 -64 -72 -81 -73 -71 -73 -71 -75 -74 -72 -75 -75];
second = [-47 -44 -51 -46 -48 -52 -65 -62 -56 -62 -62 -68 -69 -64 -68 -69 -76 -69 -72 -79 -85 -82 -77 -83 -76];
third = [-37 -39 -45 -49 -51 -53 -56 -60 -55 -63 -65 -62 -66 -67 -72 -70 -68 -90 -83 -70 -75 -71 -69 -79 -75];
fourth = [-37 -40 -42 -45 -46 -55 -60 -63 -58 -57 -67 -63 -65 -69 -66 -68 -82 -71 -74 -67 -85 -73 -69 -74 -72];
PP=[first; second;third;fourth];
PPaver=mean(PP);

%Promedio de mediciones de 0 a 25m
figure()
distancia=1:0.5:25;
PromMediciones=[-42.00 -42.00 -42.67 -45.00 -47.33 -45.33 -49.67 -49.00 -56.00 -55.33 -53.33 -60.33 -56.33 -59.67 -59.67 -57.00 -60.67 -62.33 -64.33 -67.00 -67.00 -66.33 -67.00 -66.33 -64.33 -65.67 -66.67 -71.33 -72.00 -68.67 -77.33 -78.00 -71.33 -69.67 -77.33 -72.67 -76.33 -75.00 -70.67 -85.00 -75.00 -79.00 -73.00 -74.33 -74.33 -82.00 -76.00 -74.67 -78.67];
distancia2=1:1:25;

%regresion lineal
grado=1;
polinomio=polyfit(10*log10(distancia),PromMediciones,grado);
poleval=polyval(polinomio,10*log10(distancia));
% hold on;

%Modelo 2 rayos
Pt=.001;%0dBm
Gt=1;%0dB
Gr=1;%0dB
ht=0.72;
hr=.1;
K = Pt*Gt*Gr*((ht*hr)^2);
PR2r=(K./(distancia.^4));%watts
PR2rdbm=10*log10(PR2r.*1000);%dbms

% Free space
%Model de espacio libre
f=2.4*10^9;%2.1GHertz
lambda=299792458/f;
PRel0=Gt*Gr*Pt*lambda^2/((4*pi)^2);
PRel=PRel0.*(1./(distancia.^2));%watts
PReldbm=10*log10(PRel.*1000);%conversion a dbms

%Modelo de Hata
ahr=3.2*(log10(11.54*hr)^2)-4.97; %a(hr) for large city fc>300MHz
Lurban=69.55+26.16*log10(2.4e3)-13.82*log10(ht)-ahr+(44.9-6.55*log10(ht)).*log10(distancia./1000);
LurbandB=Lurban + 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exponente de perdidas n
%PromMed = [-44 -40 -48 -48 -59 -52 -57 -61 -63 -62 -71 -61 -58 -64 -72 -81 -73 -71 -73 -71 -75 -74 -72 -75 -75];
PromMed = PPaver;
prdo=-42; % in dBm
MLnum = ones(1,25);
MLden = ones(1,25);
for N = 1:25    
     %MLnum(1,N)=(log10(N/1)*(PromMed(N) - PR2rdbm(N)));
     MLnum(1,N)=(log10(N/1)*(prdo - PR2rdbm(N)));
     MLden(1,N)=(10*((log10(N/1))^2));
end
num = 0;
den = 0;
for N2 = 1:25    
     num = num + MLnum(1,N2);
     den = den + MLden(1,N2);
end
nML = num./den;    %Exponente de perdidas n

%Modelo sencillo with Pr(do) and ML
do=1;
prdo=-42; % in dBm
%n=-polinomio(1)/10;
n=nML;
PR=Pt.*(distancia.^(-n));%watts
%PRdbm=-10*log10(PR.*1000);%dbms
%PRdbm=10*log10(PR)+30;
PRdbm=prdo-10*n*log10(distancia);


distancia3=distancia2+rand;
distancia4=distancia3-rand;
distancia5=distancia4+rand;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PP2=[first(:);second(:);third(:);fourth(:)];
distancia33=[distancia2(:);distancia3(:);distancia4(:);distancia5(:)];
[ix,iy]=sort(distancia33);
ppord=PP2(iy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    HATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fc=2.5e9;
ht=1;
hr=0.75;
dkm=ix;
aht=(1.1*log10(fc)-0.7)*hr-(1.56*log10(fc)-0.8);
%aht=3.2*(log10(11.75*hr))^2-4.97;
PL_dBm=99.55-aht+26.16*log10(fc)-13.82*log10(ht)+(44.9-6.55*log10(ht))*log10(dkm);
semilogx(ix,ppord,'r*')
%semilogx(distancia33, PP2,'*r')
%semilogx(distancia, PromMediciones,'*r')
hold on;
% semilogx(distancia2, first,'*r');
% semilogx(distancia2, second,'*r');
% semilogx(distancia2, third,'*r');
% semilogx(distancia2, fourth,'*r');
%semilogx(ix,-PL_dBm,'c--')
semilogx(distancia,-LurbandB+30,'r--')
semilogx(distancia,poleval,'b','LineWidth',2);
semilogx(distancia,PR2rdbm,'k','LineWidth',2)
semilogx(distancia,PRdbm,'m','LineWidth',2);
semilogx(distancia,PReldbm,'g','LineWidth',2);
%semilogx(distancia,LurbandB,'x');
hold off
%title('Mediciones y Modelo 2 Rayos')
grid
ylabel('Received power (dBm)')
xlabel('distance [m]')
legend('Data','Hata model','Linear regression','2-ray model','Simple model','Free space');

%pause
%%
%%%%%%%%%%%%%%%%%%%%%%%%% Analisis de diferentes Modelos
figure()
%plot(distancia,poleval,'r');
semilogx(distancia,poleval,'r')
hold on;
%Modelo sencillo
Pt=.001;%0dBm
n=polinomio(1)/10;
PR=Pt.*(distancia.^(-n));%watts
PRdbm=-10*log10(PR.*1000);%dbms

%plot(distancia,PRdbm,'b');
semilogx(distancia,PRdbm,'b')
hold on;
%Model de espacio libre
Gt=1;%0dB
Gr=1;%0dB
f=2.1*10^9;%2.1GHertz
lambda=299792458/f;
PRel0=Gt*Gr*Pt*lambda^2/((4*pi)^2);
PRel=PRel0.*(1./(distancia.^2));%watts
PReldbm=10*log10(PRel.*1000);%conversion a dbms

%plot(distancia,PReldbm,'g');
semilogx(distancia,PReldbm,'g')
hold on;
%Modelo de 2 rayos
ht=0.72;
hr=.1;
PR2r=Pt*Gt*Gr*((ht*hr)^2./(distancia.^4));%watts
PR2rdbm=10*log10(PR2r.*1000);%dbms

%plot(distancia,PR2rdbm,'m')
semilogx(distancia,PR2rdbm,'m')
hold on;
title('Graficas de mediciones y modelos de Potencias recibidas')
ylabel('dBm')
xlabel('distancia')
%Modelo de Hata
ahr=3.2*(log10(1.54*hr)^2)-4.97; %a(hr) for large city fc>300MHz
Lurban=69.55+26.16*log10(2.4e9)-13.82*log10(ht)-ahr+(44.9-6.55*log10(ht)).*log10(distancia);
LurbandB=-0.16*Lurban - 30;

%plot(distancia,LurbandB,'x');
semilogx(distancia,LurbandB,'x')
legend('Mediciones','Sencillo','Espacio Libre','Modelo 2 rayos', 'Modelo de Hata');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mapa 3D de Modelo de 2 Rayos
%Crear vectores de distancia
d = ones(25,25);
for v1 = 1:25    
   for v2 = 1:25
       d(v1,v2)=sqrt(v1^2. + v2^2.);
   end
end
%Obtenemos valores de Potencia recibida mediante Modelo 2 rayos
PRcuadrante = ones(25,25);
for v1 = 1:25    
   for v2 = 1:25
       PRcuadrante(v1,v2)=10*log10((K./(d(v1,v2).^4)).*1000);
   end
end
%Segundo Cuadrante
PRcuadranteDos = ones(25,25);
for v1 = 1:25    
   for v2 = 1:25
       PRcuadranteDos(v1,v2)=10*log10((K./(d(26-v2,v1).^4)).*1000);
   end
end
%Tercer Cuadrante
PRcuadranteTres = ones(25,25);
for v1 = 1:25    
   for v2 = 1:25
       PRcuadranteTres(v1,v2)=10*log10((K./(d(26-v2,26-v1).^4)).*1000);
   end
end
%Cuarto Cuadrante
PRcuadranteCuatro = ones(25,25);
for v1 = 1:25    
   for v2 = 1:25
       PRcuadranteCuatro(v1,v2)=10*log10((K./(d(v2,26-v1).^4)).*1000);
   end
end
figure()
%Creamos Mapa 3D de nuestra region
ModeloNeg = [PRcuadranteTres,PRcuadranteCuatro];
ModeloPos = [PRcuadranteDos,PRcuadrante];
Model3D= [ModeloNeg ; ModeloPos];
%surface(Model3D);    %Plot de plano
surf(Model3D);   %Plot en 3D
title('Potencia Recibida en Modelo 2 Rayos Ideal')
zlabel('dBm')
xlabel('distancia [m]')
ylabel('distancia [m]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Obtenemos valores de Potencia recibida mediante interpolacion
PRcuadrante1 = ones(25,25);
for v1 = 1:25    
   for v2 = 1:25
       PRcuadrante1(v1,v2)=((first(1,v1)) + (second(1,v2)))/2.;
   end
end
%Obtenemos valores de Potencia recibida mediante interpolacion
PRcuadrante2 = ones(25,25);
for v1 = 1:25    
   for v2 = 1:25
       PRcuadrante2(v2,v1)=((third(1,26-v1)) + (second(1,v2)))/2.;
   end
end
%Obtenemos valores de Potencia recibida mediante interpolacion
PRcuadrante3 = ones(25,25);
for v1 = 1:25    
   for v2 = 1:25
       PRcuadrante3(v1,v2)=((third(1,26-v1)) + (fourth(1,26-v2)))/2.;
   end
end
%Obtenemos valores de Potencia recibida mediante interpolacion 
PRcuadrante4 = ones(25,25);
for v1 = 1:25    
   for v2 = 1:25
       PRcuadrante4(v1,v2)=((fourth(1,26-v1)) + (first(1,v2)))/2.;
   end
end
figure()
%Creamos Mapa 3D de nuestra region
medicionesNeg = [PRcuadrante3,PRcuadrante4];
medicionesPos = [PRcuadrante2,PRcuadrante1];
Mapa3D= [medicionesNeg ; medicionesPos];
surface(Mapa3D);    %Plot de plano
title('Potencia Recibida')
xlabel('distancia [m]')
ylabel('distancia [m]')
figure()
surf(Mapa3D);   %Plot en 3D
title('Potencia Recibida')
zlabel('dBm')
xlabel('distancia [m]')
ylabel('distancia [m]')



%% Desviacion Estandar
sigmaMLArr = ones(1,25);
for N = 1:25    
     sigmaMLArr(1,N)=(10*nML*log10(N/1)+(PromMed(N) - PR2rdbm(N))).^2;
end
sigmaMLArr2 =0;
for N2 = 1:25    
     sigmaMLArr2 = sigmaMLArr2 + sigmaMLArr(1,N2);
end
sigmaML = (1/25.)*sqrt(sigmaMLArr2);

sigma2Ray = std(PR2rdbm);
sigmaMeds = std(PromMediciones);

%Potencias recibidas mas ruido:
%Modelo de 2 rayos:
%prm2r=PR2rdbm+normrnd(0,sigmaMeds,1,49);
%prm2r=PR2rdbm+normrnd(0,sigmaML,1,49);
prm2r=PR2rdbm+normrnd(0,sigma2Ray,1,49);

figure
%plot(distancia,prm2r,'m');
semilogx(distancia,prm2r,'k')

hold on;
%plot(distancia,PR2rdbm,'r')
%plot(distancia, PromMediciones,'x');
semilogx(distancia, PromMediciones,'r*')
title('Graficas de mediciones y modelos de Potencias recibidas + ruido')
ylabel('dBm')
xlabel('distancia')
legend('Modelo con variable aleatoria','Mediciones');



%% Mapa 3D de Modelo de 2 Rayos CON Variable Aleatoria
%Primer cuadrante
PRcuadrante = ones(25,25);
for v1 = 1:25    
   for v2 = 1:25
       PRcuadrante(v1,v2)=10*log10((K./(d(v1,v2).^4)).*1000) + normrnd(0,sigmaMeds/2,1,1);
   end
end
%Segundo Cuadrante
PRcuadranteDos = ones(25,25);
for v1 = 1:25    
   for v2 = 1:25
       PRcuadranteDos(v1,v2)=10*log10((K./(d(26-v2,v1).^4)).*1000) + normrnd(0,sigmaMeds/2,1,1);
   end
end
%Tercer Cuadrante
PRcuadranteTres = ones(25,25);
for v1 = 1:25    
   for v2 = 1:25
       PRcuadranteTres(v1,v2)=10*log10((K./(d(26-v2,26-v1).^4)).*1000) + normrnd(0,sigmaMeds/2,1,1);
   end
end
%Cuarto Cuadrante
PRcuadranteCuatro = ones(25,25);
for v1 = 1:25    
   for v2 = 1:25
       PRcuadranteCuatro(v1,v2)=10*log10((K./(d(v2,26-v1).^4)).*1000) + normrnd(0,sigmaMeds/2,1,1);
   end
end

figure()
%Creamos Mapa 3D de nuestra region
ModeloNeg = [PRcuadranteTres,PRcuadranteCuatro];
ModeloPos = [PRcuadranteDos,PRcuadrante];
Model3D= [ModeloNeg ; ModeloPos];
%surface(Model3D);    %Plot de plano
surf(Model3D);   %Plot en 3D
title('Potencia Recibida en Modelo 2 Rayos')
zlabel('dBm')
xlabel('distancia [m]')
ylabel('distancia [m]')
   
   
   