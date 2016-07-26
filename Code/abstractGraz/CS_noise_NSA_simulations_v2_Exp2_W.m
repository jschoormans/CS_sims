% NEW SIMULATIONS FOR CS 

clear all; close all; clc; 
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/exportfig'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/tightfig/'))
CC=clock; CCC=[num2str(CC(2)),'-',num2str(CC(3)),'-',num2str(CC(4)),'-',num2str(CC(5))];
figfolder=['/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Figures/',CCC]
dir=mkdir(figfolder);cd(figfolder)

%% simulate k-space
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Code')
if true
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
accvector=[10,5,4,3,1];
noisevector=[0 1 10 20 50]*1e-4; %(0;1/100;1/10 1/5 % of cksp)
for ii=1:length(accvector);
for jj=1:length(noisevector)
P.jjj=3; %weighting scheme (3=normal)
P.acc=accvector(ii);
P.NoiseLevel=noisevector(jj);
[KD{1,ii,jj}, KD{2,ii,jj}, KD{3,ii,jj}]=makeNoisyKspace(K,P);
end
end

%% Recons
reg=0.005;

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
        SSIM(2,ii,jj)=ssim(abs(squeeze(R{2,ii,jj})),abs(ImRef));
        hfen(2,ii,jj)=HFEN(abs(squeeze(R{2,ii,jj})),abs(ImRef))
    end
end

%% Visualize
cd(figfolder)
figure(101)
plot(1./accvector,squeeze(MSE(2,:,:)),'.-')
title('MSE')
xlabel('sampling fraction')
ylabel('MSE')
export_fig -native '1_MSE.eps'
export_fig -native '1_MSE.png'

figure(102)
plot(1./accvector,squeeze(SSIM(2,:,:)),'.-')
title('SSIM')
xlabel('acceleration')
ylabel('SSIM')
export_fig -native '2_SSIM.eps'
export_fig -native '2_SSIM.png'

figure(103)
plot(1./accvector,squeeze(hfen(2,:,:)),'.-')
title('HFEN')
xlabel('acceleration')
ylabel('HFEN')
export_fig -native '3_HFEN.eps'
export_fig -native '3_HFEN.png'

figure(104); count=1;
for ii=1:length(accvector); for jj=1:length(noisevector); subplot(length(accvector),length(noisevector),count); count=count+1;
imshow(abs(R{2,ii,jj}),[]); axis off; end; end;tightfig;
export_fig -native '4_RECONS.eps'
export_fig -native '4_RECONS.png'

figure(105); count=1;
for ii=1:length(accvector); for jj=1:length(noisevector); subplot(length(accvector),length(noisevector),count); count=count+1;
imshow(abs(R{2,ii,jj}([150:200],[80:130])),[]); axis off; end; end;tightfig;
export_fig -native '5_RECONS.eps'
export_fig -native '5_RECONS.png'

figure(106)
plot(noisevector,squeeze(SSIM(2,:,:)).','.-')
title('SSIM')
xlabel('noise')
ylabel('SSIM')
export_fig -native '6_SSIM.eps'
export_fig -native '6_SSIM.png'

