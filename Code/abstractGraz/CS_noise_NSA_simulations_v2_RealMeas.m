% NEW SIMULATIONS FOR CS 
clear all; close all; clc; 
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/exportfig'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/tightfig/'))
addpath(genpath('/opt/amc/matlab/toolbox/MRecon-3.0.519'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/'))
addpath(genpath('/opt/amc/bart')); vars
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab'))
CC=clock; CCC=[num2str(CC(2)),'-',num2str(CC(3)),'-',num2str(CC(4)),'-',num2str(CC(5))];
figfolder=['/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Figures/',CCC]
dir=mkdir(figfolder);cd(figfolder)

%% load k-space
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Code')
folder='/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/Reconstructions/kiwi/' 
% file='ph_25072016_2030563_13_2_wipt1wffekiwi12810dynclearV4.raw' 
% file='ph_25072016_2034581_15_2_wipt1wffekiwi128100dynclearV4.raw' %actually 256 10
% file='ph_25072016_2036025_16_2_wipt1wffekiwi256100dynclearV4.raw'
% file='ph_25072016_2133129_28_2_wipt1wffekiwi51210dynclearV4.raw' %512 10 dyns
file='ph_25072016_2125045_27_2_wipt1wffekiwi256100dynclearV4.raw'  
MR=MRecon([folder,file]); 
MR.Parameter.Parameter2Read.typ=1;
MR.Parameter.Parameter2Read.chan=2;
MR.ReadData;
MR.RandomPhaseCorrection;
MR.PDACorrection;
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;
MR.SortData;
MR.Average;
MR.RemoveOversampling;
MR.RingingFilter;
MR.ZeroFill;
% MR.K2I;
% MR.ShowData
K=MR.Data;

K=K./max(K(:)); %normalize kspace
sens=ones(size(K));
%%
accvector=[10,9,8,7,5,4,3,2,1];
noisevector=[0 1 5 10 20 30 50]*1e-4; %(0;1/100;1/10 1/5 % of cksp)
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

ImRef=fftshift(squeeze(bart(['pics -RW:7:0:',num2str(reg),' -S -e -i100'],K,sens)),2); %to compare image
for ii=1:length(accvector);
    for jj=1:length(noisevector)
        if KD{3,ii,jj}.jjj==3 || KD{3,ii,jj}.jjj==4
            R{2,ii,jj}=fftshift(bart(['pics -RW:7:0:',num2str(reg),' -S -e -i100'],KD{2,ii,jj},sens),2)
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
imshow(abs(R{2,ii,jj}([200:300],[200:300])),[]); axis off; end; end; tightfig; 
export_fig -native '5_RECONS.eps'
export_fig -native '5_RECONS.png'

figure(106)
plot(noisevector,squeeze(SSIM(2,:,:)).','.-')
title('SSIM')
xlabel('noise')
ylabel('SSIM')
export_fig -native '6_SSIM.eps'
export_fig -native '6_SSIM.png'

