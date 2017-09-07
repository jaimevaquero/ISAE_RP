function [robs,tobs,pac]=Formulation1A(LengthTime,LengthSpace,direction,tsrc,SpaceX,SpaceY,SpaceZ,FonctionV1,FonctionV2,FonctionV3,FonctionP,x,dt,c0,p0,rho0);
%
%fprintf(' xobs = %g yobs = %g zobs = %g\n',x(1),x(2),x(3))
%
rgaz    = 287;
T0      = 293;
%
ddt =2*dt;%/ dn ; 
qdt =4*dt;%/ dn ; 
sdt =6*dt;%/ dn ; 
dt2 =dt*dt;
ddt2=ddt*ddt;
   %
   % INITIALISATION % JVG (Les Xnmi sont utilisés pour les dérivées, ainsi que les variables de dt montrées juste avant)
   % pression
   pnm4=zeros(LengthSpace,1);
   pnm3=pnm4;
   pnm2=pnm3;
   pnm1=pnm2;
   p   =pnm1;
   % densite
   rhonm4=zeros(LengthSpace,1);
   rhonm3=rhonm4;
   rhonm2=rhonm3;
   rhonm1=rhonm2;
   rho   =rhonm1;
   % source coordinates
   ynm4=zeros(LengthSpace,3);
   ynm3=ynm4;
   ynm2=ynm3;
   ynm1=ynm2;
   y   =ynm1;
   % flow velocity at surface
   unm4=zeros(LengthSpace,3);
   unm3=unm4;
   unm2=unm3;
   unm1=unm2;
   u   =unm1;
   % surface velocity
   vnm4=zeros(LengthSpace,3);
   vnm3=vnm4;
   vnm2=vnm3;
   vnm1=vnm2;
   v   =vnm1;
   % normal flow velocity at surface
   un_nm4=zeros(LengthSpace,1);
   un_nm3=un_nm4;
   un_nm2=un_nm3;
   un_nm1=un_nm2;
   un    =un_nm1;
   % normal surface velocity
   vn_nm4=zeros(LengthSpace,1);
   vn_nm3=vn_nm4;
   vn_nm2=vn_nm3;
   vn_nm1=vn_nm2;
   vn    =vn_nm1;
   % normal outward pointing normal
   n_nm4=zeros(3,1); % JVG La direction ne change pas, sinon il faut ajouter la dŽpendance avec ynm
   n_nm3=n_nm4;
   n_nm2=n_nm3;
   n_nm1=n_nm2;
   n    =n_nm1;
   % Flow rate
   Unm4=zeros(3,1);
   Unm3=Unm4; % JVG La direction ne change pas, sinon il faut ajouter la dŽpendance avec ynm
   Unm2=Unm3;
   Unm1=Unm2;
   U   =Unm1;
   % Force 
   Lnm4=zeros(3,1);
   Lnm3=Lnm4; % JVG La direction ne change pas, sinon il faut ajouter la dŽpendance avec ynm
   Lnm2=Lnm3;
   Lnm1=Lnm2;
   L   =Lnm1;
   %
   % INITIALISATION
tobs=zeros(LengthSpace,LengthTime);pac=tobs;
   %
   for k=1:LengthTime
      %
      %fprintf(' Reading  %s / %s \n',sprintf('%03d',k),sprintf('%03d',nt))
      %
      % STOCKAGE HISTORIQUE POUR d/dt
      pnm4=pnm3;
      pnm3=pnm2;
      pnm2=pnm1;
      pnm1=p;
      rhonm4=rhonm3;
      rhonm3=rhonm2;
      rhonm2=rhonm1;
      rhonm1=rho;
      for i=1:3
         ynm4(:,i)=ynm3(:,i);
         ynm3(:,i)=ynm2(:,i);
         ynm2(:,i)=ynm1(:,i);
         ynm1(:,i)=y(:,i);       % JVG Quel sens ça a? Si la surface se déplace?
      end 
      for i=1:3
         unm4(:,i)=unm3(:,i);
         unm3(:,i)=unm2(:,i);
         unm2(:,i)=unm1(:,i);
         unm1(:,i)=u(:,i);
      end 
      for i=1:3
         vnm4(:,i)=vnm3(:,i);
         vnm3(:,i)=vnm2(:,i);
         vnm2(:,i)=vnm1(:,i);
         vnm1(:,i)=v(:,i);
      end 
      un_nm4=un_nm3;
      un_nm3=un_nm2;
      un_nm2=un_nm1;
      un_nm1=un;
      vn_nm4=vn_nm3;
      vn_nm3=vn_nm2;
      vn_nm2=vn_nm1;
      vn_nm1=vn;
      n_nm4=n_nm3;
      n_nm3=n_nm2;
      n_nm2=n_nm1;
      n_nm1=n;
      Unm4=Unm3;
      Unm3=Unm2;
      Unm2=Unm1;
      Unm1=U;
      Lnm4=Lnm3;
      Lnm3=Lnm2;
      Lnm2=Lnm1;
      Lnm1=L;
      % 
      if (abs(direction)==1)
         for j=1:LengthSpace
            y(j,1)  = SpaceX(j);
            if (direction == 1)
                y(j,2)  = SpaceY(end);    % JVG c'est un vecteur? Peut-être pas si abs dir == 1. Seulement valable pour symmétrie de la surface par rapport aux champs des variables. Lorsque après on défini les normales, on tient compte de si on est sur SpaceY min ou max.
            elseif (direction == -1)
                y(j,2)  = SpaceY(1);
            end
            y(j,3)  = SpaceZ;
%             u(j,1) = 0;%FonctionV1(j,k) ;
%             u(j,2) = 0;%FonctionV2(j,k) ;
%             u(j,3) = 0;%FonctionV3(j,k) ;  % JVG pourquoi 0??
            u(j,1) = FonctionV1(j,k);
            u(j,2) = FonctionV2(j,k);
            u(j,3) = FonctionV3(j,k);
            un(j)  = u(j,2);
            v(j,1) = 0;
            v(j,2) = 0;
            v(j,3) = 0;
            vn(j)  = v(j,2);
            p(j)  = FonctionP(j,k) ;% + p0;
            rho(j) = p(j) / (rgaz*T0);
         end
      elseif (abs(direction)==2)
         for j=1:LengthSpace
             if (direction == 2)
                y(j,1)  = SpaceX(end);    % JVG c'est un vecteur? Peut-être pas si abs dir == 1. Seulement valable pour symmétrie de la surface par rapport aux champs des variables. Lorsque après on défini les normales, on tient compte de si on est sur SpaceY min ou max.
            elseif (direction == -2)
                y(j,1)  = SpaceX(1);
             end
            y(j,2)  = SpaceY(j);
            y(j,3)  = SpaceZ;
            u(j,1) = FonctionV1(j,k) ;
            u(j,2) = FonctionV2(j,k) ;
            u(j,3) = FonctionV3(j,k) ;
            un(j)  = u(j,1);
            v(j,1) = 0;
            v(j,2) = 0;
            v(j,3) = 0;
            vn(j)  = v(j,1);
            p(j)  = FonctionP(j,k) ;
            rho(j) = p(j) / (rgaz*T0);
         end
      end 
      %
      for j=1:LengthSpace
         %
         % INITIALISATIONS 
         dpdt=0;        % JVG on utilise dpdt??
         dudt=dpdt;
         r=dpdt;cost=r; % JVG on utilise cost??
         M=r;
         Mi=zeros(3,1);n=Mi;ri=Mi;rhi=Mi;
         dLdt=zeros(3,1);dUdt=dLdt;dndt=dLdt;dudt=dLdt;
         Un   = M ;Undt = M ;Udtn = M ;Mr   = M ;Mrdt = M ;Lr   = M ;Lrdt = M ;LM   = M ;
         %
         % Gradient de pression ======================================================
         %
         % Ordre 4 : (valeurs stockées a i-2)    % JVG Ici, dpdt = 0 pour les deux premiers instants temporels 
%          if (k==3)
%             dpdt=(p(j)-pnm1(j))/dt;
%          elseif ( k==4 )
%             dpdt=(p(j)-pnm2(j))/ddt;
%          elseif ( k>4 )
%             dpdt=(4/3)*(pnm1(j)-pnm3(j))/ddt + (-1/3)*(p(j)-pnm4(j))/qdt;
%          end
            
           if k == 2
               dpdt = (p(j) - pnm1(j))/dt;
           elseif k == 3
               dpdt = (1.5*p(j) -2.*pnm1(j) +0.5*pnm2(j))/dt;
           elseif k == 4
               dpdt = (11./6.*p(j) -3.*pnm1(j) + 1.5*pnm2(j) -1./3.*pnm3(j))/dt;
           else
               dpdt = (25./12.*p(j) - 4.*pnm1(j) + 3.*pnm2(j) - 4./3.*pnm3(j) + 0.25*pnm4(j))/dt;
           end
        
         %
         % VECTEUR NORMAL ============================================================
         if (direction==1)
            n(1) = 0;
            n(2) = 1;
            n(3) = 0;
         elseif (direction==2)
            n(1) = 1;
            n(2) = 0;
            n(3) = 0;
         elseif (direction==-1)
            n(1) =  0;
            n(2) = -1;
            n(3) =  0;
         elseif (direction==-2)
            n(1) = -1;
            n(2) = 0;
            n(3) = 0;
         end
         %
         %
         % VECTEUR RAYON =============================================================
         for i=1:3
            %ri(i) = x(i) - ynm2(j,i);   % JVG pourquoi ynm2? ynm2 prend
            %une valeur pour la premiere fois lorsque k = 3
            ri(i) = x(i) - y(j,i);
         end 
         %r  = norm(ri);
         r  = sqrt( ri(1)^2 + ri(2)^2 + ri(3)^3 );
         %
         for i=1:3
            %rhi(i) = (x(i) - ynm2(j,i)) / r;   % JVG pourquoi ynm2? ynm2 prend une valeur pour la premiere fois lorsque k = 3
            rhi(i) = (x(i) - y(j,i)) / r;
         end 
         %
         cost = dot(n,rhi);
         %
         % y(i) : source coordinates
         % u(i) : flow velocity at the surface
         % v(i) : velocity of the surface
         %
         for i=1:3
            U(i) = rho(j)*u(j,i) / rho0 + v(j,i)*(1.0 - rho(j)/rho0) ;                  % JVG Les variables stockées ˆ i
            L(i) = (p(j) - p0)*n(i) + rho(j)*u(j,i)*(un(j) - vn(j)) ;
            %U(i) = rhonm2(j)*unm2(j,i) / rho0 + vnm2(j,i)*(1.0 - rhonm2(j)/rho0) ;
            %L(i) = (pnm2(j) - p0)*n(i) + rhonm2(j)*unm2(j,i)*(un_nm2(j) - vn_nm2(j)) ; % JVG Les variables stockées ˆ i-2
         end
         %
         % DERIVEES TEMPORELLES ======================================================
         % 
         % Ordre 4 : (valeurs stockées a i-2)
%          if (k==3)
%             for i=1:3
%                dUdt(i)  = (U(i)   - Unm1(i)   )/dt;
%                dLdt(i)  = (L(i)   - Lnm1(i)   )/dt;
%                dndt(i)  = (n(i)   - n_nm1(i)  )/dt;
%                dudt(i)  = (u(j,i) - unm1(j,i) )/dt;
%             end
%          elseif ( k==4 )
%             for i=1:3
%                dUdt(i)  = (U(i)   - Unm2(i)   )/ddt;
%                dLdt(i)  = (L(i)   - Lnm2(i)   )/ddt;
%                dndt(i)  = (n(i)   - n_nm2(i)  )/ddt;
%                dudt(i)  = (u(j,i) - unm2(j,i) )/ddt;
%             end
%          elseif ( k>4 )
%             for i=1:3
%                dUdt(i) = (4/3)*(Unm1(i)   - Unm3(i)   )/ddt + (-1/3)*(U(i)   - Unm4(i)   )/qdt;
%                dLdt(i) = (4/3)*(Lnm1(i)   - Lnm3(i)   )/ddt + (-1/3)*(L(i)   - Lnm4(i)   )/qdt;
%                dndt(i) = (4/3)*(n_nm1(i)  - n_nm3(i)  )/ddt + (-1/3)*(n(i)   - n_nm4(i)  )/qdt;
%                dudt(i) = (4/3)*(unm1(j,i) - unm3(j,i) )/ddt + (-1/3)*(u(j,i) - unm4(j,i) )/qdt;
%             end
%          end

           for i = 1:3
                if k == 2
                    dUdt(i) = (U(i) - Unm1(i))/dt;
                    dLdt(i) = (L(i) - Lnm1(i))/dt;
                    dndt(i) = (n(i) - n_nm1(i))/dt;
                    dudt(i) = (u(j,i) - unm1(j,i))/dt;
                elseif k == 3
                    dUdt(i) = (1.5*U(i) -2.*Unm1(i) +0.5*Unm2(i))/dt;
                    dLdt(i) = (1.5*L(i) -2.*Lnm1(i) +0.5*Lnm2(i))/dt;
                    dndt(i) = (1.5*n(i) -2.*n_nm1(i) +0.5*n_nm2(i))/dt;
                    dudt(i) = (1.5*u(j,i) -2.*unm1(j,i) +0.5*unm2(j,i))/dt;
                elseif k == 4
                    dUdt(i) = (11./6.*U(i) -3.*Unm1(i) + 1.5*Unm2(i) -1./3.*Unm3(i))/dt;
                    dLdt(i) = (11./6.*L(i) -3.*Lnm1(i) + 1.5*Lnm2(i) -1./3.*Lnm3(i))/dt;
                    dndt(i) = (11./6.*n(i) -3.*n_nm1(i) + 1.5*n_nm2(i) -1./3.*n_nm3(i))/dt;
                    dudt(i) = (11./6.*u(j,i) -3.*unm1(j,i) + 1.5*unm2(j,i) -1./3.*unm3(j,i))/dt;
                else
                    dUdt(i) = (25./12.*U(i) - 4.*Unm1(i) + 3.*Unm2(i) - 4./3.*Unm3(i) + 0.25*Unm4(i))/dt;
                    dLdt(i) = (25./12.*L(i) - 4.*Lnm1(i) + 3.*Lnm2(i) - 4./3.*Lnm3(i) + 0.25*Lnm4(i))/dt;
                    dndt(i) = (25./12.*n(i) - 4.*n_nm1(i) + 3.*n_nm2(i) - 4./3.*n_nm3(i) + 0.25*n_nm4(i))/dt;
                    dudt(i) = (25./12.*u(j,i) - 4.*unm1(j,i) + 3.*unm2(j,i) - 4./3.*unm3(j,i) + 0.25*unm4(j,i))/dt;
                end
           end
            
         %
         % NOMBRES DE MACH ===========================================================
         for i=1:3
            %Mi(i)    = unm2(j,i) / c0;  % JVG pour quoi unm2 et dudt?
            Mi(i)    = u(j,i) / c0;
            dMidt(i) = dudt(i)   / c0;
         end 
         M    = norm(Mi);
         %
         Un   = dot(U,n);
         Undt = dot(U,dndt);            % JVG est-il le produit scalaire la correcte opération?
         Udtn = dot(dUdt,n);
         Mr   = dot(Mi,rhi);
         Mrdt = dot(dMidt,rhi);
         Lr   = dot(L,rhi);
         Lrdt = dot(dLdt,rhi);
         LM   = dot(L,Mi);
         %
         % integration sur s ou x,y,z ? 
         if (abs(direction)==1)
            if (j==1)
               %dSN = 0.5*(ynm2(2,1) - ynm2(1,1));           % JVG pour quoi ynm2?
               dSN = 0.5*(y(2,1) - y(1,1));
            elseif (j==LengthSpace)
               %dSN = 0.5*(ynm2(j,1) - ynm2(j-1,1));
               dSN = 0.5*(y(j,1) - y(j-1,1));
            else
               %dSN = 0.5*(ynm2(j+1,1) - ynm2(j-1,1));
               dSN = 0.5*(y(j+1,1) - y(j-1,1));
            end
         elseif (abs(direction)==2)
            if (j==1)
               %dSN = 0.5*(ynm2(2,2) - ynm2(1,2));
               dSN = 0.5*(y(2,2) - y(1,2));
            elseif (j==LengthSpace)
               %dSN = 0.5*(ynm2(2,2) - ynm2(1,2));
               dSN = 0.5*(y(j,2) - y(j-1,2));
            else
               %dSN = 0.5*(ynm2(2,2) - ynm2(1,2));
               dSN = 0.5*(y(j+1,2) - y(j-1,2));
            end
         end
         dS = dSN;%*dST;
         %
         % THICKNESS NOISE ===========================================================
         % 
         pQ1 = rho0*( Udtn + Undt ) / ( r*(1.0 - Mr)^2 );
         pQ1 = (1.0/(4.0*pi))*pQ1*dS;
         % 
         pQ2 = rho0*Un*( r*Mrdt + c0*Mr - c0*M^2 ) / ( r^2*(1.0 - Mr)^3 );
         pQ2 = (1.0/(4.0*pi))*pQ2*dS;
         %
         pQ = pQ1 + pQ2 ;
         %
         % LOADING NOISE =============================================================
         %
         pL1 = Lrdt / ( r*(1.0 - Mr)^2 );
         pL1 = (1.0/(4.0*pi*c0))*pL1*dS;
         % 
         pL2 = Lr*( r*Mrdt + c0*Mr - c0*M^2 ) / ( r^2*(1.0 - Mr)^3 );
         pL2 = (1.0/(4.0*pi*c0))*pL2*dS;
         %
         pL3 = (Lr - LM) / ( r^2*(1.0 - Mr)^2 );
         pL3 = (1.0/(4.0*pi))*pL3*dS;
         %
         pL = pL1 + pL2 + pL3 ; 
         %
         % TOTAL NOISE ===============================================================
         %
         robs(j,k)    = r;
         %tobs(j,k)    = tsrc(k) + r/c0;
         %pac(j,k)     = pL + pQ ;
         %if (k>2)
         %tobs(j,k)    = tsrc(k-2) + r/c0;
         tobs(j,k)    = tsrc(k)  + r/c0;
         %tobs(j,k)    = tsrc(k) - 2 *dt  + r/c0;
         pac(j,k)     = pL + pQ ;
         %else
         %tobs(j,k)    = 0;
         %pac(j,k)     = 0;
         %end
      end    % Boucle Profil
   end       % Boucle Temps 
