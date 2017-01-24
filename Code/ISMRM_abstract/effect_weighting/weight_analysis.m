% GOAL: SYSTEMATICALLY DEFINE EFFECT OF ADDING A COVARIANCE MATRIX IN
% L2_NORM IN CS EQN. AND RECONSTRUCTION
%DO THIS WITH A SYNTHETIC KSPACE

clear all close all clc;

%synthetic kspace folder
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/MRIPhantomv0-8'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Code'));
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/CS Excercise'))
addpath(genpath('/opt/amc/bart-0.3.00')); vars;
%params; 

N_Coils=1
res_traj=64

%make kspace
% [Ksamp,sensmaps]=genPhantomKspace(res_traj,N_Coils);

Ksamp=fftshift(fftshift(fft2(phantom(256)),1),2)

%%
% add noise 
acc=6
nave=10;
nnoise=30

for i=1:nnoise
    for j=1:2
        for q=1:nave
        if j==1;
            P.noNSAcorr=true ;
        else
            P.noNSAcorr=false ;
        end
P.NoiseLevel=1e-3*(i-1);

P.acc=acc;
P.jjj=2
[K_N, Ku_Nvar1,Ku_N2,P]=makeNoisyKspace(Ksamp,P); %UNABLE TO DEAL WITH COILS 
W=P.MNSA;

Ku=squeeze(Ku_N2);
Ku=permute(Ku,[4 1 2 5 3]);

figure(99); imshow(abs(W),[]);
figure(100);
% reconstruct with and without W 
P=struct;
P.outeriter=6;
P.Itnlim=8;

P.TVWeight=(3e-4);
P.TGVfactor=0;
P.xfmWeight=0;
P.squareksp=true;
P.resultsfolder='/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Code/ISMRM_abstract/effect_weighting'

P=reconVarNSA((Ku),P) %UNABLE TO DEAL WITH 2D DATA 

R{i,j,q}=P.recon;
        end
    end
end

%% 
%analysis MSE: 
GoldenRecon=phantom;
clear Error
for i=1:nnoise
    for j=1:2
        for q=1:nave
            % Error(i,j)=sqrt(sum(abs(R{i,j}(:)-GoldenRecon(:)).^2))
            r=abs((R{i,j,q}));
%             r=r./max(r(:));
%             Error(i,j,q)=ssim(abs(r),abs(GoldenRecon))
          Error(i,j,q)=immse(GoldenRecon,r)
        end
    end
end

%%
NoiseLevel=1e-3*[0:1:29]
figure;
plot(NoiseLevel,max(Error,[],3)); legend('no W','W')
title('effect of  weighting matrix on reconstruction quality  (5x undersampled, center averaging)')
xlabel('added noise level')
ylabel('SSIM')