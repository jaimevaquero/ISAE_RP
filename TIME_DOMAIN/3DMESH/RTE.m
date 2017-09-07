% ENTREES
% SORTIES 
%
function [TAU] = RTE(Tini,C0,P,W,Xsrc,Ysrc,Tobs,X1,X2,X3,eps);
   %
   R   = sqrt(Xsrc^2 + Ysrc^2);
   Psi = atan( Ysrc/Xsrc );
   %
   % INITIALISATION :
   % POUR NEWTON RAPHSON, NE PAS AVOIR TEMPS INITIAL TROP ELOIGNE
   % DOIT ETRE LE PLUS PRES POSSIBLE DE LA RACINE
   % AUTRE SOLUTION (cf. Numerical Recipes) : 
   % COMBINER NEWTON-RAPHSON avec methode Bisection (dichotomie)
   %TAU_P    = Tini;
   %TAU_P    = 100*P;
   %TAU_P    = 101*P;
   fonction = 100.0;
   % 
   D=sqrt( (X1-Xsrc)^2 + (X2-Ysrc)^2 + X3^2 );
   TAU_P = Tini - D / C0;
   %
   nn=0;
   %
   while (abs(fonction)>eps)
      %
      nn=nn+1;
      %fprintf(' Tobs : %g \n',Tobs)
      %fprintf(' ecart = %g (TAU = %g) \n',abs(fonction),TAU_P)
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
   %
   %if (nn>1000)
   %break
   %end
      %
   end
   %
   TAU = TAU_PP;
   %print *,' NEWTON : OK'
   %
end
