function [Aext_u,Aext_v,Aext_w,Aext_p,nzs,xsExt,ysExt,zsExt,elesurfExt] = ...
    Extrusion2Dto3Dz(nxs,nt,FonctionU,FonctionV,FonctionW,FonctionP,xs,ys,zs,ExtrusionDirection,nzExtr,dzExt);

fprintf('\n Executing Extrusion in z \n')
fprintf('  Direction k:   %3g \n',ExtrusionDirection)
fprintf('  Points:        %3g \n',nzExtr)
fprintf('  Spacing:       %3.1g \n',dzExt)
fprintf('  Length:        %3.1g \n',(nzExtr-1)*dzExt)

nzs      = nzExtr;
% xsExt    = zeros(nxs,nzExtr);
xsExt    = zeros(1,nzExtr);
ysExt    = zeros(nxs,nzExtr);
zsExt    = zeros(nzExtr,1);
Aext_u   = zeros(nxs,nzExtr,nt); 
Aext_v   = zeros(nxs,nzExtr,nt);
Aext_w   = zeros(nxs,nzExtr,nt);
Aext_p   = zeros(nxs,nzExtr,nt);

xsExt(1,:) = xs;%

for i = 1:nxs
%     xsExt(i,:)    = xs(i);
    ysExt(i,:)    = ys(i);
end

for k = 1:nt
    for i = 1:nxs
        Aext_u(i,:,k)   = FonctionU(i,k);
        Aext_v(i,:,k)   = FonctionV(i,k);
        Aext_w(i,:,k)   = FonctionW(i,k);
        Aext_p(i,:,k)   = FonctionP(i,k);
    end
end

zsExt(1) = zs;
for i = 2:nzExtr
    zsExt(i) = zsExt(i-1) + ExtrusionDirection * dzExt;
end

elesurfZ      = zeros(nzExtr,1);
elesurfZ(1)   = 0.5*(zsExt(2) - zsExt(1));
elesurfZ(end) = 0.5*(zsExt(end) - zsExt(end-1));
for i = 2:nzExtr-1
    elesurfZ(i) = 0.5*(zsExt(i+1)-zsExt(i-1));
end

% It should be changed when considering a closed surface!
% elesurfXY      = zeros(nxs,1);
% elesurfXY(1)   = 0.5*((xs(2)-xs(1))^2. + (ys(2)-ys(1))^2.)^0.5;
% elesurfXY(end) = 0.5*((xs(end)-xs(end-1))^2. + (ys(end)-ys(end-1))^2.)^0.5;
% for i = 2:nxs-1
%     elesurfXY(i) = 0.5*((xs(i+1)-xs(i-1))^2. + (ys(i+1)-ys(i-1))^2.)^0.5;
% end

elesurfXY      = zeros(nxs,1);
elesurfXY(1)   = 0.5*(ys(2)-ys(1));
elesurfXY(end) = 0.5*(ys(end)-ys(end-1));
for i = 2:nxs-1
    elesurfXY(i) = 0.5*(ys(i+1)-ys(i-1));
end


elesurfExt = zeros(nxs,nzExtr);
for j = 1:nzExtr
    for i = 1:nxs
        elesurfExt(i,j) = elesurfXY(i)*elesurfZ(j);
    end
end

fprintf('\n Extrusion Completed \n\n')
end