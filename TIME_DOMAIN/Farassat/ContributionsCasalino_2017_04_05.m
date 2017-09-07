clear all ; close all ; clc ; 
%
lac  = 10;
c0   = 340;
dt   = lac/(100.0 * c0);
%
N1=1;N2=4;nt=100*N2;
robs_real_bis = dlmread('robs.dat');robs_real = repmat(robs_real_bis,N1,N2); 
pobs_real_bis = dlmread('pac.dat'); %pobs_real = repmat(pobs_real_bis,N1,N2); 
tobs_real_bis = dlmread('tobs.dat');%tobs_real = repmat(tobs_real_bis,N1,N2); 
xs=dlmread('./../../BDD_py/xs.dat');nxs=length(xs);
ys=dlmread('./../../BDD_py/ys.dat');nys=length(ys);
for j=1:nt
% for i=1:nxs
for i=1:nxs
xs_bis(i,j) = xs(i);
end
for i = nxs+1:2*nxs
    i_bis = mod(i-1,nxs) + 1;
    xs_bis(i,j) = xs(i_bis) + 70 - 30;
end
for i = 2*nxs+1:3*nxs
    i_bis = mod(i-1,nxs) + 1;
    xs_bis(i,j) = xs(i_bis) - 70 + 30;
end
for i = 3*nxs+1:4*nxs
    i_bis = mod(i-1,nxs) + 1;
    xs_bis(i,j) = xs(i_bis) + 2*(70 - 30) ;
end
end
% for j=1:100
for j = 1:nt
% for i=1:nxs
for i=1:(2*nxs + 2*nys)
pobs_real(i,j) = pobs_real_bis(i,j);
tobs_real(i,j) = tobs_real_bis(i, j);
end
end
% REPRODUCTION DES RESULTATS 4 FOIS DANS LE TEMPS
% for j=101:200
% j_bis=mod((j-1),100) + 1;
% % for i=1:nxs
% for i=1:(2*nxs + 2*nys)
% pobs_real(i,j) = pobs_real_bis(i,j_bis);
% tobs_real(i,j) = tobs_real_bis(i,j_bis) + (tobs_real_bis(i,end)-tobs_real_bis(i,1)) + (tobs_real_bis(i,2)-tobs_real_bis(i,1)); % JVG on étend le temps de 4fois. Par contre sur les boucles 201:300 et 301:nt je pense que il ne faut pas multiplier dt(2-1) par 2 ou 3. On fait périodique de période 100
% %tobs_real(i,j) = tobs_real(i,j-1) + (tobs_real_bis(i,100) -tobs_real_bis(i, 1));
% end
% end
% %figure;plot(tobs_real(101,:))
% for j=201:300
% j_bis=mod((j-1),100) + 1;
% % for i=1:nxs
% for i=1:(2*nxs + 2*nys)
% pobs_real(i,j) = pobs_real_bis(i,j_bis);
% tobs_real(i,j) = tobs_real_bis(i,j_bis) + 2*((tobs_real_bis(i,end)-tobs_real_bis(i,1))+ (tobs_real_bis(i,2)-tobs_real_bis(i,1)));
% %tobs_real(i,j) = tobs_real(i,j-1) + (tobs_real_bis(i,100) -tobs_real_bis(i, 1));
% end
% end
% for j=301:nt
% j_bis=mod((j-1),100) + 1;
% % for i=1:nxs
% for i=1:(2*nxs + 2*nys)
% pobs_real(i,j) = pobs_real_bis(i,j_bis);
% tobs_real(i,j) = tobs_real_bis(i,j_bis) + 3*((tobs_real_bis(i,end)-tobs_real_bis(i,1))+ (tobs_real_bis(i,2)-tobs_real_bis(i,1)));
% %tobs_real(i,j) = tobs_real(i,j-1) + (tobs_real_bis(i,100) -tobs_real_bis(i, 1));
% end
% end
%for j=101:nt
%for i=1:nxs
%%block=floor(j/100)
%j_bis=mod(j,100)+1;
%tobs_real(i,j) = tobs_real_bis(i, j_bis) + tobs_real(i,j-1);
%end
%end
figure
%subplot(2,1,1)
%pcolor(xs_bis(:,1:100),tobs_real(:,1:100),pobs_real(:,1:100));shading interp;colorbar;
%xlabel(' xs (m) ')
%ylabel(' Tobs real (s) ')
%subplot(2,1,2)
%pcolor(xs_bis(:,101:nt),tobs_real(:,101:nt),pobs_real(:,101:nt));shading interp;colorbar;
pcolor(xs_bis,tobs_real,pobs_real);shading interp;colorbar;
xlabel(' xs (m) ')
ylabel(' Tobs real (s) ')
%break

Tmin = min(min(tobs_real));kmin=floor(Tmin/dt);
Tmax = max(max(tobs_real));kmax=floor(Tmax/dt);
% rmax =  sqrt( (50-30)^2 + (100-30)^2 );             % JVG 95 - 19.5? 100 au lieu de 95? 19.5 doit être la valeur de ys, si on ne varie que xs. Par contre si on ne change pas xs on na pas linfluence dautres cotés
% rmin =  sqrt( (50-50)^2 + (100-70)^2 );
rmax =  sqrt( (50-30)^2 + (100-30)^2 );             % JVG 95 - 19.5? 100 au lieu de 95? 19.5 doit être la valeur de ys, si on ne varie que xs. Par contre si on ne change pas xs on na pas linfluence dautres cotés
rmin =  sqrt( (50-50)^2 + (100-70)^2 );
Tc   = (kmax-kmin)*dt - (rmax-rmin)/c0;
fprintf(' Tmin        = %g s \n',Tmin)
fprintf(' Tmax        = %g s \n',Tmax)
fprintf(' Kmin        = %g   \n',kmin)
fprintf(' Kmax        = %g   \n',kmax)
fprintf(' Borne inf. (index)  = %g \n', floor(kmin + rmax/(c0*dt) -kmin))
fprintf(' Borne sup. (index)  = %g \n', floor(kmax + rmin/(c0*dt) -kmin))
fprintf(' Borne inf. (time )  = %g \n', kmin*dt + rmax/(c0) -kmin*dt)  % JVG d'ou ces valeures???Check
fprintf(' Borne sup. (time )  = %g \n', kmax*dt + rmin/(c0) -kmin*dt)
fprintf(' Tc          = %g s \n',Tc)
%
N=size(tobs_real);ns=N(1);nt=N(2);
%
%pobs_discret=zeros(ns,20000);wij=pobs_discret;
pobs_discret=zeros(ns,ns*nt);wij=pobs_discret;
%
for j = 1:nt
   for i = 1:ns
      jadv = floor( (tobs_real(i,j)) / (dt) );
      w    =        (tobs_real(i,j)) / (dt)    - jadv;
      %fprintf(' jadv = %i \n', jadv);
      %fprintf(' w = %g \n', w);
      %pause
      %toto(i,j) = jadv;%*dt;
      %tobs_discret(jadv) = jadv;%*dt;
      if (pobs_discret(i,jadv) == 0)
          pobs_discret(i,jadv) = pobs_real(i,j) ;
          wij(i,jadv)          = w;
      else
         if ((wij(i,jadv) - w) == 0)
         fprintf('jadv = %i \n',jadv);
         end
         pw                   = (pobs_discret(i,jadv) - pobs_real(i,j))/( wij(i,jadv) - w );
         pobs_discret(i,jadv) = pobs_real(i,j) - pw*w;
         wij(i,jadv)          = 0;
      end
   end
end
%figure
%plot(toto)
%figure
%pcolor(xs_bis,tobs_real,toto);shading interp;
%xlabel(' xs (m) ')
%ylabel(' Tobs real (s) ')
%
%
Nadv=size(pobs_discret);ntadv=Nadv(2);nsadv=Nadv(1);
%%
pac_bis=zeros(ntadv,1);
%for k=(1014-100)+1:1014
%for jadv = 1:ntadv
% for jadv = floor(kmin):floor(kmax)
for jadv = floor(rmax/(c0*dt)):floor(kmax )
for i=1:nsadv
   pac_bis(jadv)=pac_bis(jadv) + pobs_discret(i,jadv);
end
end
%pac=pac_bis;
pac = pac_bis(pac_bis~=0); % JVG Retourne les valeurs non nulles de pac_bis
%figure
%plot(pac)
%
p0   = 1.0e5;
xo=dlmread('./../../POSTRAITEMENT_py/xobs.dat');nxo=length(xo);
yo=dlmread('./../../POSTRAITEMENT_py/yobs.dat');nyo=length(yo);
for k=29:428
fid=fopen(['./../../ANALYTIQUE_py/pac.T',sprintf('%03d',k),'.bin'],'r');arg=fread(fid,1,'float');pana_bis=fread(fid,[nxo,nyo],'double');fclose(fid);
pana(k-28) = pana_bis(101,201)-p0;
end

figure
plot([0:(4/(400-1)):4],pana);%,'k-','linewidth',2);
hold on
%for j=1:nsadv
%plot([0:(1/(length(pij(j,:))-1)):N2],pij(j,:));
%end
plot([0:(N2/(length(pac)-1)):N2],pac,'r-');
%plot([0:(N2/(length(pac)-1)):N2],-lac*pac,'r-');

legend('Analytique','Farassat','Location','NorthEast')

