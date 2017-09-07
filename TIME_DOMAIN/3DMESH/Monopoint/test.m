clear all ; close all ; clc ; 
%
% ROTOR =======================================================================
RPM = 5000 ;
W   = RPM * 2.0 * pi / 60;
P   = 60 / RPM ;
% OBSERVATEUR =================================================================
X1   = 0.1; 
X2   = 0.1;
nt   = 101 ; Tobs=zeros(nt,1);
Tobs = [2*P : P/(nt-1) : 3*P];
%
% FIGURE ======================================================================
%
xmin = -1.5*sqrt(X1^2 + X2^2);ymin=xmin;
xmax =  1.5*sqrt(X1^2 + X2^2);ymax=xmax;
%
startx=2;
starty=8;
sizex=15;
sizey=15;
nb_fig=1;
%
fig = figure(nb_fig);
set(fig,'visible','on');
orient landscape;
set(fig,'papertype','a4letter');
set(fig,'units','centimeters','paperunits','centimeters');
set(fig,'paperposition',[startx starty sizex sizey]);
set(fig,'position',[startx starty sizex sizey]);
startx=startx+0.5;starty=starty-0.5;nb_fig=nb_fig+1;
set(gcf,'color','white')
%
for nn=1:nt
%
R   = sqrt( X1^2 + X2^2 );
phi = atan( X2/X1 ); 
%
Y1 = R*cos( W*Tobs(nn) +phi );
Y2 = R*sin( W*Tobs(nn) +phi );
%Zsrc(j,k) = 0.0 ;
%end
%end
%
plot(X1,X2,'o','color','r','markerfacecolor','r');
hold on
plot(Y1,Y2,'o','color','b','markerfacecolor','b');
set(gca,'Xlim',[xmin xmax],'Ylim',[ymin ymax]);
axis off;
pause;
hold off;
%
end
