% NEW SIMULATIONS FOR CS 

clear all; close all; clc; 
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/exportfig'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/tightfig/'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/CS Exercise'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/MRIPhantomv0-8'))
CC=clock; CCC=[num2str(CC(2)),'-',num2str(CC(3)),'-',num2str(CC(4)),'-',num2str(CC(5))];
figfolder=['/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Figures/',CCC]
dir=mkdir(figfolder);cd(figfolder)

%% simulate k-space
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Code')
if false
Im=imread('head256.jpg');   %image in image domain
Im=double(Im(:,:,1))./255;
K=bart('fft 3',Im);
sens=ones(size(K));

else
[K,sens]=genPhantomKspace(128,1);
sens=permute(sens,[1 2 4 3]);
end
K=K./max(K(:)); %normalize kspace
%%
accvector=[10,9,8,7,6,5,4,3,2,1];
noisevector=[0 1 5 7 10 15 20 40 80 100 120 150]*1e-4; %(0;1/100;1/10 1/5 % of cksp)
for ii=1:length(accvector);
for jj=1:length(noisevector)
P.jjj=3; %weighting scheme (3=normal)
P.acc=accvector(ii);
P.NoiseLevel=noisevector(jj);
[KD{1,ii,jj}, KD{2,ii,jj}, KD{3,ii,jj}]=makeNoisyKspace(K,P);
end
end

%% Recons
reg=0.01;

ImRef=squeeze(bart(['pics -RW:7:0:',num2str(reg),' -S -e -i100'],K,sens)); %to compare image
for ii=1:length(accvector);
    for jj=1:length(noisevector)
        if KD{3,ii,jj}.jjj==3 || KD{3,ii,jj}.jjj==4
            R{2,ii,jj}=bart(['pics -RW:7:0:',num2str(reg),' -S -e -i100'],KD{2,ii,jj},sens)
        else 
        %other recon
        end
    end
end
%% Estimate errors
clear MSE SSIM hfen
for ii=1:length(accvector);
    for jj=1:length(noisevector)
        MSE(2,ii,jj)=mean(mean(abs(squeeze(R{2,ii,jj})-ImRef).^2,1),2);
        SSIM(2,ii,jj)=ssim(double(abs(squeeze(R{2,ii,jj}))),double(abs(ImRef)));
        PSNR(2,ii,jj)=psnr(double(abs(squeeze(R{2,ii,jj}))),double(abs(ImRef)));

    end
end

%% calculate SNR

sigcoordsx=[7]
sigcoordsy=[58:70]
noisecoordsx=[1:20]
noisecoordsy=[1:20]

for ii=1:length(accvector);
    for jj=1:length(noisevector)
ROISignal= abs(squeeze(R{2,ii,jj}(sigcoordsx,sigcoordsy)))
ROINoise =abs(squeeze(R{2,ii,jj}(noisecoordsx,noisecoordsy)))
SNR{2,ii,jj}=mean(ROISignal(:))/std(ROINoise(:))
    end
end

ii=length(accvector) ;for jj=1:length(noisevector)
SNRvector(jj)=SNR{2,ii,jj}
end
%% Visualize
cd(figfolder)
figure(101)
plot(sqrt(SNRvector(2:end)),squeeze(PSNR(2,[2,5,8,10],1:end)).','.-')
set(gcf,'color','w');

xlabel('sqrt SNR')
ylabel('PSNR')
export_fig -native '1_SSIM.eps'
export_fig -native '1_SSIM.png'
legend('9x','6x','3x','full','Location','NorthWest')
%%
figure(102)
plot(1./accvector,squeeze(PSNR(2,:,[1,4,6:8])),'*-');
set(gcf,'color','w');

title('PSNR')
xlabel('NSA/undersampling')
ylabel('PSNR')
export_fig -native '2_PSNR.eps'
export_fig -native '2_PSNR.png'
legend('SNR=inf','SNR=150','SNR=70','SNR=50','SNR=25')


figure(1022)
plot(1./accvector,squeeze(SSIM(2,:,[1,4,6:8])),'*-');
set(gcf,'color','w');

title('SSIM')
xlabel('NSA/undersampling')
ylabel('SSIM')
export_fig -native '2_SSIM.eps'
export_fig -native '2_SSIM.png'
legend('SNR=inf','SNR=150','SNR=70','SNR=50','SNR=25')


figure(1023)
plot(1./accvector,squeeze(MSE(2,:,[1,4,6:8])),'*-');
set(gcf,'color','w');

title('MSE')
xlabel('NSA/undersampling')
ylabel('MSE')
export_fig -native '2_MSE.eps'
export_fig -native '2_MSE.png'
legend('SNR=inf','SNR=150','SNR=70','SNR=50','SNR=25')

%%
figure(501)
plot(1./accvector,squeeze(PSNR(2,:,[1])),'*-');
ylim([45 65])
set(gcf,'color','w');
title('PSNR')
xlabel('NSA/undersampling')
ylabel('PSNR')
export_fig -native '2_PSNR_1.eps'
export_fig -native '2_PSNR_1.png'

%legend('SNR=inf','SNR=150','SNR=70','SNR=50','SNR=25')




%%
figure(104); count=1;
for ii=1:length(accvector); for jj=1:length(noisevector); subplot(length(accvector),length(noisevector),count); count=count+1;
imshow(abs(R{2,ii,jj}),[]); axis off; end; end;tightfig;
export_fig -native '4_RECONS.eps'
export_fig -native '4_RECONS.png'

figure(105); count=1;
for ii=1:length(accvector); for jj=1:length(noisevector); subplot(length(accvector),length(noisevector),count); count=count+1;
imshow(abs(R{2,ii,jj}(1:20,1:20)),[]); axis off;
end; end;tightfig;
export_fig -native '5_RECONS.eps'
export_fig -native '5_RECONS.png'



