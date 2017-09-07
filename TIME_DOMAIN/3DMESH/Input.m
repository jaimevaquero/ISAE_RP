%
%BLADE = 'MAVion NB2';
%BLADE = 'MAVion NB3';
%BLADE = 'MAVion NB4';
%BLADE = 'DCAS NB2';
%BLADE = 'DCAS NB3';
BLADE = 'DCAS NB4';
%BLADE = 'Fukano NB2';
%BLADE = 'optim';
%BLADE = 'arbitrary';
%
switch BLADE
   case 'MAVion NB2'
      NCP       = 4; CP_R=zeros(NCP,1);
      CP_R(NCP) = 0.0875;
      CP_R(1)   = 0.00875;
      CP_R      = [CP_R(1):(CP_R(NCP)-CP_R(1))/(NCP-1):CP_R(NCP)];
      CP_C      = [0.025 0.025 0.015 0.015];
      CP_Fx     = [0.000 0.000 0.000 0.000];
      NC        = 5;
      RPM  = 7250 ;
   case 'MAVion NB3'
      NCP       = 4; CP_R=zeros(NCP,1);
      CP_R(NCP) = 0.0875;
      CP_R(1)   = 0.00875;
      CP_R      = [CP_R(1):(CP_R(NCP)-CP_R(1))/(NCP-1):CP_R(NCP)];
      CP_C      = [0.045 0.015 0.035 0.025];
      CP_Fx     = [0.000 0.000 0.000 0.000];
      NC        = 5;
      RPM  = 5000 ;
   case 'MAVion NB4'
      NCP       = 4; CP_R=zeros(NCP,1);
      CP_R(NCP) = 0.0875;
      CP_R(1)   = 0.00875;
      CP_R      = [CP_R(1):(CP_R(NCP)-CP_R(1))/(NCP-1):CP_R(NCP)];
      CP_C      = [0.025 0.015 0.025 0.005];
      CP_Fx     = [0.000 0.000 0.000 0.000];
      NC        = 5;
      RPM  = 6000 ;
   case 'DCAS NB2'
      NCP       = 4; CP_R=zeros(NCP,1);
      CP_R(NCP) = 0.127;
      CP_R(1)   = 0.00875;
      CP_R      = [CP_R(1):(CP_R(NCP)-CP_R(1))/(NCP-1):CP_R(NCP)];
      CP_C      = [0.025 0.03  0.02  0.01];
      CP_Fx     = [0.000 0.00  0.00  0.00];
      NC        = 5;
      RPM  = 4200 ;
   case 'DCAS NB3'
      NCP       = 4; CP_R=zeros(NCP,1);
      CP_R(NCP) = 0.127;
      CP_R(1)   = 0.00875;
      CP_R      = [CP_R(1):(CP_R(NCP)-CP_R(1))/(NCP-1):CP_R(NCP)];
      CP_C      = [0.04  0.01  0.01   0.01];
      CP_Fx     = [0.000 0.00  0.00  0.00];
      NC        = 5;
      RPM  = 4400 ;
   case 'DCAS NB4'
      NCP       = 4; CP_R=zeros(NCP,1);
      CP_R(NCP) = 0.127;
      CP_R(1)   = 0.00875;
      CP_R      = [CP_R(1):(CP_R(NCP)-CP_R(1))/(NCP-1):CP_R(NCP)];
      CP_C      = [0.01  0.01  0.04   0.01];
      CP_Fx     = [0.000 0.00  0.00  0.00];
      NC        = 5;
      %RPM  =62500 ;
      RPM  =22500 ;
      %RPM  =16500 ;
      %RPM  =12500 ;
      %RPM  = 2900 ;
   case 'Fukano NB2'
      NCP       = 4; CP_R=zeros(NCP,1);
      CP_R(NCP) = 0.0875;
      CP_R(1)   = 0.00875;
      CP_R      = [CP_R(1):(CP_R(NCP)-CP_R(1))/(NCP-1):CP_R(NCP)];
      CP_C      = [0.05  0.05  0.05   0.05];
      CP_Fx     = [0.000 0.000 0.000 0.000];
      NC        = 5;
      RPM  = 6300 ;
   case 'optim'
      NCP       = 4; CP_R=zeros(NCP,1);
      CP_R(NCP) = 0.0875;
      CP_R(1)   = 0.00875;
      CP_R      = [CP_R(1):(CP_R(NCP)-CP_R(1))/(NCP-1):CP_R(NCP)];
      CP_C      = [0.045 0.040 0.040 0.035];
      CP_Fx     = [0.040 0.040 0.040 -0.040];
      NC        = 5;
      RPM  = 5000 ;
   case 'arbitrary'
      NCP       = 4; CP_R=zeros(NCP,1);
      CP_R(NCP) = 0.0875;
      CP_R(1)   = 0.00875;
      CP_R      = [CP_R(1):(CP_R(NCP)-CP_R(1))/(NCP-1):CP_R(NCP)];
      CP_C      = [0.045 0.040 0.040 0.035];
      CP_Fx     = [0.040 0.040 0.040 -0.040];
      NC        = 5;
      RPM  = 5000 ;
end
%
% ROTATION SPEED ==========================================
%
% OBSERVER ================================================
Robs       = 2.00 ;   % Observer Distance (m)
theta      = 90.0 ;   % Elevation angle   (degres)
phi        = 20.0 ;   % Azimuthal angle   (degres)
%theta      = 90.0 ;   % Elevation angle   (degres)
%phi        = 00.0 ;   % Azimuthal angle   (degres)
