function [K_N, Ku_Nvar1,P]=makeNoisyKspacefromdynamics(K,P)
% load k-space data (2D-dyn)
%outputs K_N: noisy full k-space 
%Ku_Nvar1: undersampled and averaged k-space

K=squeeze(K);
%make masks
rng('default');
rng(1)

[pdf,~] = genPDF(size(K(:,:,1)),4,1/P.acc,2,0,0);
Mfull=genEllipse(size(K,1),size(K,2));
M=genSampling(pdf,10,100).*Mfull;
P.pdf=pdf;
P.M=M;

MNSA=makeMNSA(P);
if max(MNSA(:))>size(K,3)
error('not enough dynamics')

end
sum(M(:).*MNSA(:))
P.MNSA=MNSA; %export MNSA
P.M=M;

%% add noise to kspace
for iii=1:max(MNSA(:)) %Matrix of NSA values
K_N=K(:,:,iii);
Ku_N1=squeeze(M(:,:).*(MNSA(:,:)>=iii)).*K_N;
Ku_N2(:,:,1,iii)=Ku_N1;
end
Ku_Nvar1=sum(Ku_N2,4)./MNSA; 


end 