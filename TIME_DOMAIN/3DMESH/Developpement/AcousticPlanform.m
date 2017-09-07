clear all ; close all ; clc ; 
%
% BLADE GEOMETRIE / ROTATION SPEED / OBSERVER =================================
Input
% MAILLAGE ====================================================================
Rx = 10; % Facteur de raffinement
N = Rx*NCP;dR=(CP_R(NCP)-CP_R(1))/(N-1);R=zeros(N,1);
R = [CP_R(1):dR:CP_R(NCP)];
C=zeros(N,1);Fx=C;
%
C  = bezier(CP_R,CP_C,R);
Fx = bezier(CP_R,CP_Fx,R);
%
XG=zeros(NC,N);YG=XG;ZG=XG;Rsrc=XG;
% Effets Corde
for k=1:N
   XG(1,k) = R(k);
   YG(1,k) = C(k)/2;
   ZG(1,k) = 0.0 ; 
   dY     = C(k)/(NC-1);
   for j=2:NC
      XG(j,k) = R(k);
      YG(j,k) = YG(j-1,k) - dY; 
      ZG(j,k) = 0.0 ; 
   end
end
%
% Effets Fleche
for k=1:N
   for j=1:NC
      YG(j,k) = YG(j,k) - Fx(k); 
   end
end
%
% Chord MAX
Cmax = abs(YG(1,N) - YG(NC,N));
%
% PROPAGATION =================================================================
C0 = 340.0 ;                % Speed of Sound
% ROTOR =======================================================================
W    = RPM * 2.0 * pi / 60; % Rotational speed
P    = 60 / RPM ;           % Revolution time 
Mach = W*CP_R(NCP) / C0 ;      % Mach number at blade tip
strMach=sprintf('%2.2f',Mach);
% OBSERVATEUR =================================================================
Rdum       = 1.1*R(end) ;         % Fake observer Distance for plotting
strRobs=sprintf('%2.2f',Robs);
theta      = theta * pi / 180 ;   % Elevation angle (rad)
phi        = phi   * pi / 180 ;   % Azimuthal angle (rad)
X1         = Robs * sin(theta)*cos(phi) ;
X2         = Robs * sin(theta)*sin(phi) ;
X3         = Robs * cos(theta)          ;
X1dummy    = Rdum * sin(theta)*cos(phi) ;
X2dummy    = Rdum * sin(theta)*sin(phi) ;
X3dummy    = Rdum * cos(theta)          ;
nt   = 101 ; Tobs=zeros(nt,1);AC=Tobs;Xtip=Tobs;Ytip=Tobs;Ztip=Tobs;
Tobs = [100*P : P/(nt-1) : 101*P];
%
% RETARDED-TIME EQUATION ======================================================
Tsrc=zeros(NC,N,nt);
for nn=1:nt
   for k=1:N
      for j=1:NC
         Xsrc=XG(j,k);Ysrc=YG(j,k);
         [Tsrc_bis] = RTE(C0,P,W,Xsrc,Ysrc,Tobs(nn),X1,X2,X3,1e-5);
         Tsrc(j,k,nn) = Tsrc_bis;
      end
   end
end
%
% FIGURE ======================================================================
%
xmin = -1.1*Rdum;ymin=xmin;zmin=xmin;
xmax =  1.1*Rdum;ymax=xmax;zmax=xmax;
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
   for k=1:N
      for j=1:NC
         Rsrc      = sqrt( XG(j,k)^2 + YG(j,k)^2 );
         Phi       = atan( YG(j,k) / XG(j,k) );
         Xsrc(j,k) = Rsrc*cos( W*Tsrc(j,k,nn) +Phi );
         Ysrc(j,k) = Rsrc*sin( W*Tsrc(j,k,nn) +Phi );
         Zsrc(j,k) = 0.0 ;
      end     
   end
   %Acoustic Chord at blade tip
   AC(nn)=abs( Ysrc(NC,N) - Ysrc(1,N) );
   %Blade tip position for the acoustic circle
   Xtip(nn)=Xsrc((NC+1)/2,N);Ytip(nn)=Ysrc((NC+1)/2,N);Ztip(nn)=Zsrc((NC+1)/2,N);
   %
   %plot3([0 1],[0 0],[0 0],'k-')
   %plot3([0 0],[0 1],[0 0],'k-')
   %plot3([0 0],[0 0],[0 1],'k-')
   %%
   plot3(Xsrc,Ysrc,Zsrc,'k-',Xsrc',Ysrc',Zsrc','k-')
   hold on
   plot3(XG,YG,ZG,'-','color',0.6*[1 1 1])
   plot3(XG',YG',ZG','-','color',0.6*[1 1 1])
   %
   plot3(Xsrc(1:NC,N),Ysrc(1:NC,N),Zsrc(1:NC,N),'b-','linewidth',2);
   %
   plot3(Xtip(1:nn),Ytip(1:nn),Ztip(1:nn),'ko','markerfacecolor','k','markersize',3);
   %
   plot3(X1dummy,X2dummy,X3dummy,'o','color','r','markerfacecolor','r');
   title(['Blade ' BLADE '    Tip Mach : ' strMach '    Observer distance : ' strRobs ' m ' ],'FontWeight','bold')
   set(gca,'Xlim',[xmin xmax],'Ylim',[ymin ymax],'Zlim',[zmin zmax]);
   axis off;
   az = 00;
   al = 90;
   view(az,al);
   pause(0.1);
   hold off;
   %
end
text(-Rdum,Rdum,0,[' Maximum Acoustic Chord (relative to Geometric Chord) : ' num2str(max(AC(:))/Cmax*100)])
%
%fprintf(' Maximum Acoustic Chord (relative to Geometric Chord) : %3.1f \n',max(AC(:))/Cmax*100)
%
%plot(AC/C(end)*100,'b-','linewidth',2);
%set(gca,'Xlim',[1 nt])
