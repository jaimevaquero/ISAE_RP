clear all ; close all ; clc ; 

a = 10;
b = 2;
c = -150;

x = [0 : 0.01 : 10];nx=length(x);f=zeros(nx,1);

f = a.*x.^2 + b.*x + c ; 
%
eps      = 1e-3;
%%% DICHOTOMIE =================================
%ecart=1000;
%a = x(1); 
%b = x(nx);
%%
%while (ecart>1000*eps);
%   %
%   i = find(abs(x(:)-a)==min(abs(x(:)-a)));
%   fa = f(i);
%   i = find(abs(x(:)-b)==min(abs(x(:)-b)));
%   fb = f(i);
%   %
%   c = (a + b)/2; 
%   %
%   i = find(abs(x(:)-c)==min(abs(x(:)-c)));
%   fc = f(i);
%   %
%   if (fa<0 && fb<0)
%      disp('Change interval bounds\n')
%   elseif(fa<0 && fb>0)
%      if (fc<0)
%         a = c;
%      else
%         b = c;
%      end
%   elseif(fa>0 && fb<0)
%      if (fc<0)
%         b = c;
%      else
%         a = c;
%      end
%   end
%   ecart=abs(fc)
%end
disp('Dichotomie OK')
% NEWTON - RAPHSON ===========================
fonction = 1000;
% 
TAU_P = 0.5*abs(x(nx) - x(1));
%
nn=0;
%
while (abs(fonction)>eps)
   %
   nn=nn+1
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
