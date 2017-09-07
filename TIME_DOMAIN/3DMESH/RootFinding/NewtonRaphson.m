function [] = NewtonRaphson(xmin,xmax,eps)
%
% INITIALISATION 
fonction = 1000;
TAU_P    = 0.5*abs(xmax - xmin);
nn       = 0;
%
while (abs(fonction)>eps)
   %
   nn=nn+1;
   %
   fonction = a*TAU_P^2 + b*TAU_P + c;
   %
   fonction_p = 2*a*TAU_P + b;
   %
   % MODIF
   TAU_PP = TAU_P - fonction/fonction_p;
   %
   TAU_P = TAU_PP;
   % 
end
%
TAU = TAU_PP;
%
disp('Newton Raphson OK')
%
figure
plot(x,f)
hold on;
%plot([c c],[-200 1000],'r-')
%plot([x(1) x(nx)],[fc fc],'r-')
plot([TAU TAU],[-200 1000],'r-')
plot([x(1) x(nx)],[0 0],'r-')
