clear all ; close all ; clc ; 
%
C_ACmax = dlmread('./Chord/RESULTATS/ACmax.dat');
C_COMB  = dlmread('./Chord/RESULTATS/MATRICE.COMBINATOIRE.CHORD.dat');C_NN=size(C_COMB);C_NCB=C_NN(2);C_NCP=C_NN(1);
F_ACmax = dlmread('./Sweep/RESULTATS/ACmax.dat');
F_COMB  = dlmread('./Sweep/RESULTATS/MATRICE.COMBINATOIRE.SWEEP.dat');F_NN=size(F_COMB);F_NCB=F_NN(2);F_NCP=F_NN(1);
%
C_imax=find(C_ACmax(:)==min(C_ACmax(:)));
%
F_imax=find(F_ACmax(:)==min(F_ACmax(:)));
%
fprintf('CP C for minimum Acoustic Chord : %2.4f %2.4f %2.4f %2.4f\n',C_COMB(1,C_imax),C_COMB(2,C_imax),C_COMB(3,C_imax),C_COMB(4,C_imax))
fprintf('CP F for minimum Acoustic Chord : %2.4f %2.4f %2.4f %2.4f\n',F_COMB(1,F_imax),F_COMB(2,F_imax),F_COMB(3,F_imax),F_COMB(4,F_imax))
%for ycp=1:C_NCB
%   
%end
