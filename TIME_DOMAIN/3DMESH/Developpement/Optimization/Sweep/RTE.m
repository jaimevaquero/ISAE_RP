% ENTREES
% SORTIES 
%
function [TAU] = RTE(C0,P,W,R,Psi,Tobs,X1,X2,X3,eps);
   %
   % INITIALISATION :
   TAU_P    = 1000*P;
   fonction = 100.0;
   %
   while (abs(fonction)>eps)
      %
      %fprintf(' ecart = %g \n',abs(fonction))
      %
      Y1 = R*cos(W*TAU_P + Psi);
      Y2 = R*sin(W*TAU_P + Psi);
      Y3 = 0.0 ;
      %
      alpha = (X1 - Y1)^2;
      beta  = (X2 - Y2)^2;
      gam   = (X3 - Y3)^2;
      %
      h        = alpha + beta + gam;
      g        = sqrt(h);
      %
      fonction = TAU_P - Tobs + (g/C0);
      %
      Y1_p = - R*W*sin(W*TAU_P + Psi);
      Y2_p =   R*W*cos(W*TAU_P + Psi);
      Y3_p = 0.0 ;
      %
      alpha_p = -2.0*(X1-Y1)*Y1_p;
      beta_p  = -2.0*(X2-Y2)*Y2_p;
      gam_p   = -2.0*(X3-Y3)*Y3_p;
      %
      h_p        = alpha_p + beta_p + gam_p;
      g_p        = h_p/(2.0*g);
      %
      fonction_p = 1.0 + (g_p/C0);
      %
      % MODIF
      TAU_PP = TAU_P - fonction/fonction_p;
      %
      TAU_P = TAU_PP;
      %
   end
   %
   TAU = TAU_PP;
   %print *,' NEWTON : OK'
   %
end
