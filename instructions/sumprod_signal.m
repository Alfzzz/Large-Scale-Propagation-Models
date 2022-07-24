clear all; close all; clc;

% Final time in seconds
tf=10;
% for plotting final time
tfp=tf/4;
% Initial time in seconds
to=0;
% Low frequency in Hz
fL=1;
% High frequency in Hz
fh=20;
% Amplitude high frequency signal
ah=0.75;
% Amplitude low frequency signal
aL=2;

% Sampling
ts=1/(10*fh);

% Sum of a low frequncy signal with a high frequency signal of zero mean
% value
x=(to:ts:tf)';
yh=ah*cos(2*pi*fh*x);
yL=aL*cos(2*pi*fL*x)+(aL/2)*cos(2*pi*(fL/2)*cos(x));
z=yh+yL;
plot(x,z,'b',x,yL,'r')
xlabel('Time $t$ (sec)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
ylabel('$y_h + y_L$','Interpreter','latex','FontSize',18,'FontWeight','bold','Color','k')
grid
axis([to tfp min(z) max(z)])

pause

figure
% Product of a low frequncy signal with a high frequency signal of zero mean
% value
x=(to:ts:tf)';
yhp=ah*cos(2*pi*fh*x);
yLp=aL*cos(2*pi*fL*x)+(aL/2)*cos(2*pi*(fL/2)*cos(x));
zp=yhp.*yLp;
plot(x,zp,'b',x,ah*yLp,'r')
xlabel('Time $t$ (sec)','Interpreter','latex','FontSize',16,'FontWeight','bold','Color','k')
ylabel('$y_h y_L$','Interpreter','latex','FontSize',18,'FontWeight','bold','Color','k')
grid
axis([to tfp min(zp) max(zp)])


