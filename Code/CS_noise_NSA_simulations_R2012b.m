clear all; close all; clc; 
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations')
addpath(genpath('/opt/amc/bart')); vars;

%VARS
ny=288
nz=84

Im=imread('image.bmp');   %image in image domain
Im=double(Im(:,:,1))./255

%MAKE SENS MAPS AND K_SPACE (SIMULATION DATA)
sens=ones(1,ny,nz,1);
K=bart('fft 3',Im);

load('mask')
M=mask'
acceleration=sum(M(:)==1)/((ny/2)*(nz/2)*(pi))
% MAKE ELLIPSE MASK 
%
Mfull=genEllipse(ny,nz);

K=squeeze(Mfull').*K
Ku=squeeze(M).*K;      %undersampled k-space%
Ku=permute(Ku,[3 2 1]);
%CS
ImCS=bart('pics -RW:6:0:0.1 -i100',Ku,sens);
Imlin=ifftshift(ifft2((squeeze(Ku))));

figure(1);
subplot(221); imshow(Im'); subplot(222); imshow(squeeze(M))
subplot(223); imshow(abs(squeeze(Imlin)),[]); subplot(224); imshow(abs(squeeze(ImCS)),[])
%%
clear K_N Ku_N

% from experiment: 
SNR_estimator =0.03 %(one coil tho)
E_k_mean=mean(mean(mean(abs(Ku).^2)))
E_I=sum(sum(sum(abs(Im(:,:,:)).^2,1),2),3) % in experiment: Ei =Ek (verschil door ellipse?)
sigma=E_k_mean*SNR_estimator

for i=1:4  %MAKE  four noisy undersampled meas (and one fully sampled noisy measurement)
K_N=addNoise(K,sqrt(sigma));
Ku_N(:,:,:,:,i)=squeeze(M).*K_N;
end
K_N=permute(K_N,[3 2 1]);
Ku_N=permute(Ku_N,[3 2 1,4,5]);
Ku_N4=mean(Ku_N,5);

ImCS_Nus=bart('pics -RW:6:0:0.01 -e -i100',Ku_N(:,:,:,1),sens);
ImCS_Nus4=bart('pics -RW:6:0:0.01 -e -i100',Ku_N4,sens);
ImCS_Nfs=bart('pics -RW:6:0:0.01 -e -i100',K_N,sens);

D1=abs(squeeze(ImCS_Nus)-Im')
D2=abs(squeeze(ImCS_Nus4)-Im');
D3=abs(squeeze(ImCS_Nfs)-Im');

figure(2);
subplot(231); imshow(abs(squeeze(ImCS_Nus)),[0 1]); 
title('1 dyn')
subplot(232); imshow(abs(squeeze(ImCS_Nus4)),[]); 
title('4 dyn')
subplot(233); imshow(abs(squeeze(ImCS_Nfs)),[0 1]);
title('full s')

subplot(234); hold on;imshow(D1,[0 1]); text(20,20,num2str(sum(D1(:))),'Color','yellow'); hold off
subplot(235); hold on;imshow(D2,[0 1]); text(10,20,num2str(sum(D2(:))),'Color','yellow'); hold off
subplot(236); hold on;imshow(D1,[0 1]); text(10,20,num2str(sum(D3(:))),'Color','yellow'); hold off

%% DO SAME EXPERIMENT FOR DIFFERENT NOISE LEVELS 
NoiseLevels=[[1:1:9],[10:2:30],[35:5:100]]
[Recon, D]=analyse_NS(Im,K,M,sens,NoiseLevels);

%%
figure(3);
plot(NoiseLevels,D','+-');
legend('one NSA','4 NSA','fully sampled')
xlabel('noise std added in k-space (a.u.)')
ylabel('abs sum of pixel differences')

figure(4);
cc=1
for ii=[1,5,15,20,25,31]
subplot(2,3,cc); cc=cc+1;
hold on 
    imshow(abs(squeeze((cat(4,Recon(ii,1,:,:),Recon(ii,2,:,:),Recon(ii,3,:,:))))),[])
    text(10,10,num2str(NoiseLevels(ii)),'Color','yellow')
end

%% DO SAME EXPERIMENT FOR DIFFERENT NOISE LEVELS  & REGS 
NoiseLevels=[[1:1:9],[10:2:30],[35:5:100]]
regs=[1e-2]
[Recon, D,D2,D3]=analyse_NSreg(Im,K,M,sens,NoiseLevels,regs);

%%
figure(5);
subplot(131)
plot(NoiseLevels,squeeze(D(:,:,1)'),'-*');
title('MSE')
legend('CS: R=4 NSA1','CS: R=4 NSA4','full')
subplot(132)
plot(NoiseLevels,squeeze(D2(:,:,1)'),'-*'); 
title('SSIM')
legend('CS: R=4 NSA1','CS: R=4 NSA4','full')
subplot(133)
plot(NoiseLevels,squeeze(D3(:,:)'),'-*'); 
title('noise estimator')
legend('CS: R=4 NSA1','CS: R=4 NSA4','full')


%%  
%WHAT IS THE EFFECT OF UNDERSAMPLING?!?!?

addpath(genpath('/home/jschoormans/lood_storage/divi/projects/cosart/Matlab Collection/utils'))

NoiseLevels=[1,10,40];
unders=[0.9:-0.05:0.05]
[Recon3, Du,Du2,Du3,acc]=analyse_NSregus(Im,K,sens,NoiseLevels,1e-2,unders); %maybe with vardens undersamplen????

%%

figure(6)
subplot(331)
plot(acc,squeeze(Du3(:,1,:))','-*');
title('MSE')
subplot(332)
plot(acc,squeeze(Du3(:,2,:))','-*');
title('MSE')
subplot(333)
plot(acc,squeeze(Du3(:,3,:))','-*');
title('MSE')

subplot(334)
plot(acc,squeeze(Du2(:,1,:))','-*');
title('SSIm')
subplot(335)
plot(acc,squeeze(Du2(:,2,:))','-*');
title('SSIM')
subplot(336)
plot(acc,squeeze(Du2(:,3,:))','-*');
title('SSIM')

subplot(337)
plot(acc,squeeze(Du3(:,1,:))','-*');
title('std')
subplot(338)
plot(acc,squeeze(Du3(:,2,:))','-*');
title('std')
subplot(339)
plot(acc,squeeze(Du3(:,3,:))','-*');
title('sdtd')

%% NOGMAALS VOOR MEER NOISELEVELS EN ZONDER NSA=4 
NoiseLevels=[1,5,10,15,20,25,30];
unders=[0.9:-0.05:0.05]
[Recon4,Dv,Dv2,Dv3,acc]=analyse_NSregus2(Im,K,sens,NoiseLevels,1e-2,unders); %maybe with vardens undersamplen????

%%

figure(7)
hold on
for ii=1:7
plot(acc,squeeze(Dv2(1,ii,:)),'r')
plot(acc,squeeze(Dv2(2,ii,:)),'b--')

end
hold off
title('SSIM for different noiselevels (CS and full sampling)')
hold off

figure(8);
subplot(511)
us=1;
imshow([squeeze(Recon4(1,2,:,:,us)) squeeze(Recon4(1,1,:,:,us)) squeeze(Recon4(2,1,:,:,us)) squeeze(Recon4(3,1,:,:,us)) squeeze(Recon4(4,1,:,:,us)) squeeze(Recon4(5,1,:,:,us)) squeeze(Recon4(6,1,:,:,us)) squeeze(Recon4(7,1,:,:,us))])
axis off
subplot(512)
us=5;
imshow([squeeze(Recon4(1,2,:,:,us)) squeeze(Recon4(1,1,:,:,us)) squeeze(Recon4(2,1,:,:,us)) squeeze(Recon4(3,1,:,:,us)) squeeze(Recon4(4,1,:,:,us)) squeeze(Recon4(5,1,:,:,us)) squeeze(Recon4(6,1,:,:,us)) squeeze(Recon4(7,1,:,:,us))])
axis off
subplot(513)
us=10;
imshow([squeeze(Recon4(1,2,:,:,us)) squeeze(Recon4(1,1,:,:,us)) squeeze(Recon4(2,1,:,:,us)) squeeze(Recon4(3,1,:,:,us)) squeeze(Recon4(4,1,:,:,us)) squeeze(Recon4(5,1,:,:,us)) squeeze(Recon4(6,1,:,:,us)) squeeze(Recon4(7,1,:,:,us))])
axis off
subplot(514)
us=14;
imshow([squeeze(Recon4(1,2,:,:,us)) squeeze(Recon4(1,1,:,:,us)) squeeze(Recon4(2,1,:,:,us)) squeeze(Recon4(3,1,:,:,us)) squeeze(Recon4(4,1,:,:,us)) squeeze(Recon4(5,1,:,:,us)) squeeze(Recon4(6,1,:,:,us)) squeeze(Recon4(7,1,:,:,us))])
axis off
subplot(515)
us=18;
imshow([squeeze(Recon4(1,2,:,:,us)) squeeze(Recon4(1,1,:,:,us)) squeeze(Recon4(2,1,:,:,us)) squeeze(Recon4(3,1,:,:,us)) squeeze(Recon4(4,1,:,:,us)) squeeze(Recon4(5,1,:,:,us)) squeeze(Recon4(6,1,:,:,us)) squeeze(Recon4(7,1,:,:,us))])
axis off

figure(9);
imshow([squeeze(Recon4(1,2,:,:,1)) squeeze(Recon4(2,2,:,:,1)) squeeze(Recon4(3,2,:,:,1)) squeeze(Recon4(4,2,:,:,1)) squeeze(Recon4(5,2,:,:,1)) squeeze(Recon4(6,2,:,:,1)) squeeze(Recon4(7,2,:,:,1))],[])
axis off

%%
us=5
figure(10); %plot one SNR with varying us

SNR=1
PLOT=[]
for iii=1:9
    us=2*iii-1
PLOT=[PLOT,squeeze(Recon4(SNR,2,:,:,us))];
end
    imshow(PLOT)
    text(10,10,strcat('SNR: ',num2str(NoiseLevels(SNR))),'Color','red');
    text(20,20,strcat('undersampling',num2str(unders(1:2:17))),'Color','red')
    axis off
    export_fig 10_plot_us_SNR1.eps -native

