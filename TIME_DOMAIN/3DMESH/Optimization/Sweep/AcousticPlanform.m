clear all ; close all ; clc ; 
%
% BLADE GEOMETRIE =============================================================
NCP       = 4; CP_R=zeros(NCP,1);
CP_R(NCP) = 0.0875;
CP_R(1)   = 0.00875;
CP_R      = [CP_R(1):(CP_R(NCP)-CP_R(1))/(NCP-1):CP_R(NCP)];
CP_C      = [0.025 0.025 0.015 0.015];
NC        = 5;
%
% RANGE OF CHORD VARIATION
C_VAR_R1 = [-0.040:0.010:0.040];
C_VAR_R2 = [-0.040:0.010:0.040];
C_VAR_R3 = [-0.040:0.010:0.040];
C_VAR_R4 = [-0.040:0.010:0.040];
% Algorithme combinatoire
% Combinaisons sur la chord
sets = {C_VAR_R1, C_VAR_R2, C_VAR_R3, C_VAR_R4};
[x y z w] = ndgrid(sets{:});
COMB_C = [x(:) y(:) z(:) w(:)];
temp=size(COMB_C);NNC=temp(1);
%
dlmwrite('./RESULTATS/MATRICE.COMBINATOIRE.SWEEP.dat',COMB_C');
%
for ycp=1:NNC
   %  
   fprintf(' SWEEP # %i / %i \n',ycp,NNC)
   %  
   for cp=1:NCP
      CP_Fx(cp) = COMB_C(ycp,cp) ;
   end
   %
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
   RPM  = 8000 ;
   %
   W    = RPM * 2.0 * pi / 60; % Rotational speed
   P    = 60 / RPM ;           % Revolution time 
   Mach = W*R(end) / C0 ;      % Mach number at blade tip
   strMach=sprintf('%2.2f',Mach);
   % OBSERVATEUR =================================================================
   Robs       = 2.00 ;   % Observer Distance (m)
   theta      = 90.0 ;   % Elevation angle   (degres)
   phi        = 00.0 ;   % Azimuthal angle   (degres)
   %
   Rdum       = 0.10 ;              % Fake observer Distance for plotting
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
   Tobs = [2*P : P/(nt-1) : 3*P];
   %
   % RETARDED-TIME EQUATION ======================================================
   Tsrc=zeros(NC,N,nt);
   for nn=1:nt
      for k=1:N
         for j=1:NC
            Rsrc      = sqrt( XG(j,k)^2 + YG(j,k)^2 );
            Phi       = atan( YG(j,k) / XG(j,k) );
            [Tsrc] = RTE(C0,P,W,Rsrc,Phi,Tobs(nn),X1,X2,X3,1e-5);
            Ysrc(j,k) = Rsrc*sin( W*Tsrc + Phi );
         end
      end
      %Acoustic Chord at blade tip
      AC(nn)=abs( Ysrc(NC,N) - Ysrc(1,N) );
   end
   %
   ACmax(ycp)=max(AC(:))/Cmax*100;
   %
end
%
dlmwrite('./RESULTATS/ACmax.dat',ACmax');
