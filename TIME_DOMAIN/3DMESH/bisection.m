function [bmid] = bisection(P,W,Tobs,R,Psi,X1,X2,X3,C0,eps);
%
ecart=1000;fbmid=ecart;
%
binf = P;
bsup = 200*P;
%
while (ecart>eps);
   %
   % BORNE INF.
   Y1binf = R*cos(W*binf + Psi);
   Y2binf = R*sin(W*binf + Psi);
   Y3binf = 0.0 ;
   r1binf = (X1 - Y1binf)^2;
   r2binf = (X2 - Y2binf)^2;
   r3binf = (X3 - Y3binf)^2;
   gbinf  = sqrt(r1binf + r2binf + r3binf);
   fbinf  = binf - Tobs + (gbinf/C0);
   % BORNE SUP.
   Y1bsup = R*cos(W*bsup + Psi);
   Y2bsup = R*sin(W*bsup + Psi);
   Y3bsup = 0.0 ;
   r1bsup = (X1 - Y1bsup)^2;
   r2bsup = (X2 - Y2bsup)^2;
   r3bsup = (X3 - Y3bsup)^2;
   gbsup  = sqrt(r1bsup + r2bsup + r3bsup);
   fbsup  = bsup - Tobs + (gbsup/C0);
   %
   % BORNE INTERMEDIAIRE
   bmid   = (binf + bsup)/2;
   Y1bmid = R*cos(W*bmid + Psi);
   Y2bmid = R*sin(W*bmid + Psi);
   Y3bmid = 0.0 ;
   r1bmid = (X1 - Y1bmid)^2;
   r2bmid = (X2 - Y2bmid)^2;
   r3bmid = (X3 - Y3bmid)^2;
   gbmid  = sqrt(r1bmid + r2bmid + r3bmid);
   fbmidnm1 = fbmid;
   fbmid  = bmid - Tobs + (gbmid/C0);
   %
   if (fbinf<0 && fbsup<0)
      disp('Change interval bounds\n')
   elseif(fbinf<0 && fbsup>0)
      if (fbmid<0)
         binf = bmid;
      else
         bsup = bmid;
      end
   elseif(fbinf>0 && fbsup<0)
      if (fbmid<0)
         bsup = bmid;
      else
         binf = bmid;
      end
   end
   ecart=abs(fbmid - fbmidnm1);
end
