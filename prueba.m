clc
clear all
close all
%suburbano
fsub=[55.25 61.25 67.25 77.25 175.25 187.25 471.25];
ht_sub=40;
hr_sub=3;
d_urban=36;
asuburbana=(1.1*log10(fsub)-0.7)*hr_sub-(1.56*log10(fsub)-0.8);
lurban1=69.55+26.16*log10(fsub)-13.82*log10(ht_sub)-asuburbana+(44.9-6.55*log10(ht_sub))*log10(d_urban);
atenua_sub=lurban1-2*(log10(fsub/28)).^2-5.4

%rural
frural=[175.25];
ht_rural=40;
hr_rural=3;
d_rural=99;
arural=(1.1*log10(frural)-0.7)*hr_rural-(1.56*log10(frural)-0.8);
lurban2=69.55+26.16*log10(frural)-13.82*log10(ht_rural)-arural+(44.9-6.55*log10(ht_rural))*log10(d_rural);
atenua_rural=lurban2-4.78*(log10(frural)).^2-18.33*log10(frural)-40.98