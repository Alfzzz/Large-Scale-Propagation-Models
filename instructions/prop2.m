clear

pt=0.5;
ple=3.5;
sigma=8;

d=(1:500)';
d=d*0.5;
[i,j]=size(d(:));

pr=pt*(d.^(-ple));
y = sigma.*randn(i,j);

prdBm_r=10*log10(pr)+30+y;

% p=polyfit(log10(d),prdBm_r,1);
% prn=polyval(p,log10(d));
p=polyfit(log10(d),prdBm_r,1);
prn=polyval(p,log10(d));

semilogx(d,prdBm_r,'r.',d,prn,'b')
xlabel('distancia (m)')
ylabel('Potencia (dBm)')
grid






