clear all ; close all ; clc ; 

p1 = 10;
p2 = 2;
p3 = -150;

x = [0 : 0.01 : 10];nx=length(x);f=zeros(nx,1);

f = p1.*x.^2 + p2.*x + p3 ; 
%
eps      = 1e-3;
%% DICHOTOMIE =================================
ecart=1000;
a = x(1); 
b = x(nx);
%
while (ecart>1000*eps);
   %
   i = find(abs(x(:)-a)==min(abs(x(:)-a)));
   fa = f(i);
   i = find(abs(x(:)-b)==min(abs(x(:)-b)));
   fb = f(i);
   %
   c = (a + b)/2; 
   %
   i = find(abs(x(:)-c)==min(abs(x(:)-c)));
   fc = f(i);
   %
   if (fa<0 && fb<0)
      disp('Change interval bounds\n')
   elseif(fa<0 && fb>0)
      if (fc<0)
         a = c;
      else
         b = c;
      end
   elseif(fa>0 && fb<0)
      if (fc<0)
         b = c;
      else
         a = c;
      end
   end
   ecart=abs(fc);
end
disp('Dichotomie OK')
% NEWTON - RAPHSON ===========================
fonction = 1000;
% 
TAU_P = 0.5*abs(x(nx) - x(1));
%
nn=0;nn_bounds=nn;
%
while (abs(fonction)>eps)
   %
   nn=nn+1;
   %
   fonction = p1*TAU_P^2 + p2*TAU_P + p3;
   %
   fonction_p = 2*p1*TAU_P + p2;
   %
   % MODIF
   TAU_PP = TAU_P - fonction/fonction_p;
   %
   TAU_P = TAU_PP;
   % 
   if (TAU_PP > TAU_P)
      nn_bounds=nn_bouds+1;
      if (nn_bounds > 20)
         break
      end
   end
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
plot([c c],[-200 1000],'g-')
plot([x(1) x(nx)],[fc fc],'g-')
plot([TAU TAU],[-200 1000],'r-')
plot([x(1) x(nx)],[0 0],'r-')
