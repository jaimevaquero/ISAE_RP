clear all ; close all ; clc ; 

% ns  : curvilinear coordinates 
% nbe : blade element
% nt  : observer time

% ENVIRONNEMENT ===================================================================
Amplitude = 5;
lac       = 10;
%
Mach = 0.0;
rho0 = 1.22;
p0   = 1.0e5;
c0   = 340;
%
dt = lac/(100.0 * c0);
%
% FACTEUR 10 ? 
%
% DOMAINE SOURCE ==================================================================
xs = dlmread('./../../BDD_py/xs.dat');nxs=length(xs);
ys = dlmread('./../../BDD_py/ys.dat');nys=length(ys);
zs = 0;
% Surfaces de controle FW-H
xsurfp=max(xs);
xsurfm=min(xs);
%ysurfp=max(ys);
%ysurfm=min(ys);
% Elements de surface
elesurfX(1)     = 0.5*( xs(2)   - xs(1)     );
elesurfX(nxs)   = 0.5*( xs(nxs) - xs(nxs-1) );
for is=2:nxs-1
   elesurfX(is) = 0.5*abs(xs(is+1) - xs(is-1));
end
%elesurfY(1)     = 0.5*( ys(2)   - ys(1)     );
%elesurfY(nys)   = 0.5*( ys(nys) - ys(nys-1) );
%for is=2:nys-1
%   elesurfY(is) = 0.5*abs(ys(is+1) - ys(is-1));
%end
%
% DOMAINE OBSERVATEUR =============================================================
xo = 50;
yo = 100;
%yo = 95;
zo = 0 ;
%dxo = Lxo/(nxo-1);
%xo(1) = 0.25;
%for io=2:nxo
%   xo(io) = xo(io-1) + dxo;
%end
%%
%yo = ys((nys+1)/2);
%%
%dzo = Lzo/(nzo-1);
%zo(1) = 0.25;
%for jo=2:nzo
%   zo(jo) = zo(jo-1) + dzo;
%end
%
%% TESTS POUR INDICE OBS DE POSITION SURFACE 
%% TOP
%jo=1;
%while ( zo(jo) < zsurfm )
%   jo=jo+1;
%end
%jo_low=jo;
%% BOTTOM
%jo=nzo;
%while ( zo(jo) > zsurfp )
%   jo=jo-1;
%end
%jo_high=jo;
%% REAR 
%io=1;
%while ( xo(io) < xsurfm )
%   io=io+1;
%end
%io_low=io;
%% FRONT
%io=nxo;
%while ( xo(io) > xsurfp )
%   io=io-1;
%end
%io_high=io;
%
% DISCRETISATION TEMPORELLE =======================================================
%
% nt=120;
nt = 400 + 728 - 300;
nzExt = 11;
for k=1:nt
   tsrc(k)      = real(k)*dt;
end
% %
fprintf('\n SURFACES TOP \n')
%
A_u_bis = dlmread(['./../../BDD_py/ux.xtp.dat']);
A_v_bis = dlmread(['./../../BDD_py/uy.xtp.dat']);
A_p_bis = dlmread(['./../../BDD_py/pp.xtp.dat']);
%
for k=1:nt
   for is=1:nxs
      k_bis=mod((k-1),nt)+1;     % JVG: Reprend de 1 à 100 puis de 1 à 20!!!
      %A_uyp(is,k) = 0;
      %A_vyp(is,k) = 0;
      A_u(is,k) = A_u_bis(is+( k_bis-1)*nxs);
      A_v(is,k) = A_v_bis(is+( k_bis-1)*nxs);
      A_w(is,k) = 0;
      A_p(is,k) = A_p_bis(is+( k_bis-1)*nxs);% + p0 ;
   end
end
%
% [T_robs1p,T_tobs1p,T_pac1p] = Formulation1A(nt,nxs,+1,tsrc,xs,ys,zs,A_u,A_v,A_w,A_p,[xo yo zo],dt,c0,p0,rho0);
[T_robs1p,T_tobs1p,T_pac1p] = Formulation1A_extrusion_JVG(nt,nxs,nzExt,+1,tsrc,xs,ys,zs,A_u,A_v,A_w,A_p,[xo yo zo],dt,c0,p0,rho0);
% %------------------------------------------------
fprintf('\n SURFACES BOTTOM \n')
%
A_u_bis = dlmread(['./../../BDD_py/ux.xtm.dat']);
A_v_bis = dlmread(['./../../BDD_py/uy.xtm.dat']);
A_p_bis = dlmread(['./../../BDD_py/pp.xtm.dat']);
%
for k=1:nt
   for is=1:nxs
      k_bis=mod((k-1),nt)+1;     % JVG: Reprend de 1 à 100 puis de 1 à 20!!!
      %A_uyp(is,k) = 0;
      %A_vyp(is,k) = 0;
      A_u(is,k) = A_u_bis(is+( k_bis-1)*nxs);
      A_v(is,k) = A_v_bis(is+( k_bis-1)*nxs);
      A_w(is,k) = 0;
      A_p(is,k) = A_p_bis(is+( k_bis-1)*nxs);% + p0 ;
   end
end
%
% [T_robs1m,T_tobs1m,T_pac1m] = Formulation1A(nt,nxs,-1,tsrc,xs,ys,zs,A_u,A_v,A_w,A_p,[xo yo zo],dt,c0,p0,rho0);
[T_robs1m,T_tobs1m,T_pac1m] = Formulation1A_extrusion_JVG(nt,nxs,nzExt,-1,tsrc,xs,ys,zs,A_u,A_v,A_w,A_p,[xo yo zo],dt,c0,p0,rho0);
% %------------------------------------------------
fprintf('\n SURFACES RIGHT \n')
%
A_u_bis = dlmread(['./../../BDD_py/ux.ytp.dat']);
A_v_bis = dlmread(['./../../BDD_py/uy.ytp.dat']);
A_p_bis = dlmread(['./../../BDD_py/pp.ytp.dat']);
%
for k=1:nt
   for is=1:nxs
      k_bis=mod((k-1),nt)+1;     % JVG: Reprend de 1 à 100 puis de 1 à 20!!!
      %A_uyp(is,k) = 0;
      %A_vyp(is,k) = 0;
      A_u(is,k) = A_u_bis(is+( k_bis-1)*nxs);
      A_v(is,k) = A_v_bis(is+( k_bis-1)*nxs);
      A_w(is,k) = 0;
      A_p(is,k) = A_p_bis(is+( k_bis-1)*nxs);% + p0 ;
   end
end
%
% [T_robs2p,T_tobs2p,T_pac2p] = Formulation1A(nt,nxs,+2,tsrc,xs,ys,zs,A_u,A_v,A_w,A_p,[xo yo zo],dt,c0,p0,rho0);
[T_robs2p,T_tobs2p,T_pac2p] = Formulation1A_extrusion_JVG(nt,nxs,nzExt,+2,tsrc,xs,ys,zs,A_u,A_v,A_w,A_p,[xo yo zo],dt,c0,p0,rho0);
% %------------------------------------------------
fprintf('\n SURFACES LEFT \n')
%
A_u_bis = dlmread(['./../../BDD_py/ux.ytm.dat']);
A_v_bis = dlmread(['./../../BDD_py/uy.ytm.dat']);
A_p_bis = dlmread(['./../../BDD_py/pp.ytm.dat']);
%
for k=1:nt
   for is=1:nxs
      k_bis=mod((k-1),nt)+1;     % JVG: Reprend de 1 à 100 puis de 1 à 20!!!
      %A_uyp(is,k) = 0;
      %A_vyp(is,k) = 0;
      A_u(is,k) = A_u_bis(is+( k_bis-1)*nxs);
      A_v(is,k) = A_v_bis(is+( k_bis-1)*nxs);
      A_w(is,k) = 0;
      A_p(is,k) = A_p_bis(is+( k_bis-1)*nxs);% + p0 ;
   end
end
%
% [T_robs2m,T_tobs2m,T_pac2m] = Formulation1A(nt,nxs,-2,tsrc,xs,ys,zs,A_u,A_v,A_w,A_p,[xo yo zo],dt,c0,p0,rho0);
[T_robs2m,T_tobs2m,T_pac2m] = Formulation1A_extrusion_JVG(nt,nxs,nzExt,-2,tsrc,xs,ys,zs,A_u,A_v,A_w,A_p,[xo yo zo],dt,c0,p0,rho0);
% %

T_robs = [T_robs1p;T_robs1m;T_robs2p;T_robs2m];
T_tobs = [T_tobs1p;T_tobs1m;T_tobs2p;T_tobs2m];
T_pac  = [T_pac1p; T_pac1m; T_pac2p; T_pac2m ];

% dlmwrite('robs.dat',T_robs2m(:,10:109));
% dlmwrite('tobs.dat',T_tobs2m(:,10:109));
% dlmwrite('pac.dat',T_pac2m(:,10:109));    % JVG on enregistre la pression sur [xo yo zo] générée par tous les points source
                                        % aux instants temporels allant de 10 a 109
% dlmwrite('robs.dat',T_robs(:,:));
% dlmwrite('tobs.dat',T_tobs(:,:));
% dlmwrite('pac.dat',T_pac(:,:));

for l=1:nzExt
    for k =1:nt
        for i=1:2*nxs+2*nys
            T_robsWr(i+(l-1)*(2*nxs+2*nys),k)=T_robs(i,l,k);
            T_tobsWr(i+(l-1)*(2*nxs+2*nys),k)=T_tobs(i,l,k);
            T_pacWr (i+(l-1)*(2*nxs+2*nys),k)=T_pac (i,l,k);
        end
   end
end

dlmwrite('robs.dat',T_robsWr(:,:));
dlmwrite('tobs.dat',T_tobsWr(:,:));
dlmwrite('pac.dat',T_pacWr(:,:));