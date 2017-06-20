function [K_N, Ku_Nvar1,P]=makeNoisyKspacefromdynamics_noav(K,P)
% load k-space data (2D-dyn)
%outputs K_N: noisy full k-space 
%Ku_Nvar1: undersampled and averaged k-space

K=squeeze(K);
%make masks
[pdf,~] = genPDF(size(K(:,:,1)),4,1/P.acc,1,0,0);
Mfull=genEllipse(size(K,1),size(K,2));
M=genSampling(pdf,10,100).*Mfull;

% make weighting scheme 
if P.jjj==1
    MNSA=ceil(1./pdf);
elseif P.jjj==2
    MNSA=ceil(pdf*P.acc);
elseif P.jjj==3
    MNSA=P.acc*ones(size(pdf));
elseif P.jjj==4
    MNSA=ones(size(pdf));
elseif P.jjj==5
    MNSA=P.acc*ones(size(pdf))*P.usedyns;
if max(MNSA)>size(K,3);
error('not enough dynamics!')  
end
end;

P.MNSA=MNSA; %export MNSA
P.M=M;
P.pdf=pdf;

%% add noise to kspace
for iii=1:max(MNSA(:)) %Matrix of NSA values
K_N=K(:,:,iii);
Ku_N1=squeeze(M(:,:).*(MNSA(:,:)>=iii)).*K_N;
Ku_N2(:,:,1,iii)=Ku_N1;
end
Ku_Nvar1=permute(Ku_N2,[1 2 3 5 6 7 8 9 10 4]); 
disp('check if this makes snense still;')



end 