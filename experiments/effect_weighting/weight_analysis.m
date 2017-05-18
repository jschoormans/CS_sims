% GOAL: SYSTEMATICALLY DEFINE EFFECT OF ADDING A COVARIANCE MATRIX IN
% L2_NORM IN CS EQN. AND RECONSTRUCTION
%DO THIS WITH A SYNTHETIC KSPACE

clear all close all clc;

%synthetic kspace folder
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/MRIPhantomv0-8'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Code'));
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/CS Excercise'))
addpath(genpath('/opt/amc/bart-0.3.00')); vars;
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/radialavg'))
%params; 

N_Coils=1;
res_traj=64;

%make kspace
% [Ksamp,sensmaps]=genPhantomKspace(res_traj,N_Coils);

Ksamp=fftshift(fftshift(fft2(phantom(256)),1),2);

%%
% add noise 
acc=3
nave=1;
nnoise=2
sigma=1e1
for i=1:nnoise
    for j=1:2
        for q=1:nave
        if j==1;
            P.noNSAcorr=true ;
        else
            P.noNSAcorr=false ;
        end
P.NoiseLevel=sigma*(i-1);

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
P.resultsfolder='/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/experiments/effect_weighting'

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
        for q=1:nave
            %% 
            for j=1:2
                
                % Error(i,j)=sqrt(sum(abs(R{i,j}(:)-GoldenRecon(:)).^2))
                r=abs((R{i,j,q}));
                %             r=r./max(r(:));
                Error(i,j,q)=ssim(abs(r),abs(GoldenRecon))
                %           Error(i,j,q)=immse(GoldenRecon,r)
                %% POWER SPECTRAL DENSITY FUNCTION
                
                PSD=fft2c(abs(r-GoldenRecon));
                
                figure(10);
                ax1= subplot(3,2,j*1);
                imshow(abs(r-GoldenRecon),[]); axis off; 
                colormap(ax1,bone); 
                if j==1; title('no W'); else; title('with W'); end
                
                ax2=subplot(3,2,2+j);
                imshow(abs(PSD)); axis off;
                colormap(ax2,jet); drawnow; 
                
                %plot radial average 
                [Zr(j,:), Rx] = radialavg(PSD,25)
                ax3=subplot(3,2,4+j) 
                plot(Rx,abs(Zr(j,:))); 
                
            end
    end
end

%%
close all
NoiseLevel=1e-3*[0:nnoise-1];
figure(1);
plot(NoiseLevel,max(Error,[],3)); legend('no W','W')
title('effect of  weighting matrix on reconstruction quality  (5x undersampled, center averaging)')
xlabel('added noise level')
ylabel('SSIM')

