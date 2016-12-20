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
folder='/home/jschoormans/lood_storage/divi/Projects/cosart/scans/kiwi/'
% folder='/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/Reconstructions/kiwi/' 
% % file='ph_25072016_2030563_13_2_wipt1wffekiwi12810dynclearV4.raw' 
% % file='ph_25072016_2034581_15_2_wipt1wffekiwi128100dynclearV4.raw' %actually 256 10
% % file='ph_25072016_2036025_16_2_wipt1wffekiwi256100dynclearV4.raw'
% % file='ph_25072016_2133129_28_2_wipt1wffekiwi51210dynclearV4.raw' %512 10 dyns
% file='ph_25072016_2125045_27_2_wipt1wffekiwi256100dynclearV4.raw'  
% file='ph_25072016_2137350_29_2_wipt1wffekiwi51230dynclearV4.raw'
% file='ph_26072016_1752392_7_2_wipt1wffe256thinsliceclearV4.raw' 

%%GRAPEFRUIT
% folder='/home/jschoormans/lood_storage/divi/Projects/cosart/scans/fruit287/'
% file='fr_28072016_1937549_9_2_wipffegrapefruit256100clearV4.raw' 

%KIWI 
% folder='/home/jschoormans/lood_storage/divi/Projects/cosart/Scans/fruit287/'
% file='fr_28072016_1952148_15_2_wipffekiwi256100clearV4.raw'
% file='fr_28072016_1949258_14_2_wipffekiwi128100clearV4.raw'

folder='/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/Reconstructions/kiwi/'
file='ph_25072016_2150347_31_2_wipt1wffekiwi200100dynnewsliceclearV4.raw'
file='ph_25072016_2157057_32_2_wipt1wffekiwi240100dynnewsliceclearV4.raw'

MR=MRecon([folder,file]); 
MR.Perform
MR.I2K
Kn=MR.Data;
Kn=Kn./max(Kn(:)); %normalize kspace
sigma=0.002;
for i=1:100
   K(:,:,1,1,1,1,1,1,1,1,1,i)=Kn+randn(size(Kn)).*sigma;
end
K=padarray(K,[8 8],0,'both');
sens=ones(size(K,1),size(K,2));
%% make ref image
reg=0.02;
Kref=padarray(Kn,[8 8],0,'both'); %still use coils though
ImRef=(squeeze(bart(['pics -RW:7:0:',num2str(reg),' -S -e -i500'],Kref,sens))); %to compare image
ImRef=abs(ImRef);
ImRef=ImRef./(max(ImRef(:)));
sens=bart('ecalib -m1 -r20',Kref); %should be [nx,ny,1,nc] ?

%% make k-spaces
% accvector=[1,1.25,1.5,2,2.5,3,4,5,6];
accvector=[3,5,7]
noisevector=[5,6,7] % in this case the different W 
for jj=1:length(noisevector)
for ii=1:length(accvector);
P.jjj=noisevector(jj); %weighting scheme (3=normal)
P.usedyns=1; %for example
P.acc=accvector(ii);
[KD{1,ii,jj}, KD{2,ii,jj}, KD{3,ii,jj}]=makeNoisyKspacefromdynamics(K,P);
end
end

%% PARAMETER USED FOR CONJUGATE GRADIENT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TVWeight = 0; 	% Weight for TV penalty
xfmWeight = 0.0002;	% Weight for Transform L1 penalty
Itnlim = 8;		% Number of iterations



%generate transform operator
XFM = Wavelet('Daubechies',4,4);	% Wavelet

% initialize Parameters for reconstruction
param = init;
param.XFM = XFM;
param.TV = TVOP;
param.TVWeight =0;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;
%%
N=[size(K,1) size(K,2)]

for ii=1:length(accvector);
    for jj=1:length(noisevector)
%         if KD{3,ii,jj}.jjj==3 || KD{3,ii,jj}.jjj==4 ||KD{3,ii,jj}.jjj==5
%             R{2,ii,jj}=fftshift(bart(['pics -RW:7:0:',num2str(reg),' -d5 -S -e -i100'],KD{2,ii,jj},sens),2)
%             R{2,ii,jj}=abs(R{2,ii,jj});
%             R{2,ii,jj}=R{2,ii,jj}./max(R{2,ii,jj}(:)); %scaling!
%         else
            %generate Fourier sampling operator
            param.data=KD{2,ii,jj};
            FT = p2DFT(KD{3,ii,jj}.M, N, 1, 2);
            param.FT = FT;
            param.xfmWeight=xfmWeight
            im_dc2 = FT'*(param.data.*KD{3,ii,jj}.M./KD{3,ii,jj}.pdf); %linear recon; scale data to prevent low-pass filtering
            res = XFM*im_dc2;
            param.V=(KD{3,ii,jj}.MNSA.*KD{3,ii,jj}.M);
            param.Debug=0;
            param.xfmWeight=xfmWeight*(mean(param.V(KD{3,ii,jj}.M~=0)))
            
            for n=1:3
                res = fnlCg_test(res,param);
                R{2,ii,jj} = XFM'*res;
            end
            R{2,ii,jj}=(abs(R{2,ii,jj}));
            R{2,ii,jj}=R{2,ii,jj}./max(R{2,ii,jj}(:));
%         end
    end
end
%% Estimate errors
clear MSE SSIM hfen
for ii=1:length(accvector);
    for jj=1:length(noisevector)

        MSE(2,ii,jj)=mean(mean(abs(squeeze(R{2,ii,jj})-ImRef).^2,1),2);
        SSIM(2,ii,jj)=ssim(abs(squeeze(R{2,ii,jj})),abs(ImRef));
        hfen(2,ii,jj)=HFEN(abs(squeeze(R{2,ii,jj})),abs(ImRef));
        PSNR(2,ii,jj)=psnr(abs(squeeze(R{2,ii,jj})),abs(ImRef))
    end
end
SSIM
%% Visualize
cd(figfolder)
figure(101)
subplot(221)
plot(1./accvector,squeeze(MSE(2,:,:)),'.-')
title('MSE')
xlabel('sampling fraction')
ylabel('MSE')

subplot(222)
plot(1./accvector,squeeze(SSIM(2,:,:)),'.-')
title('SSIM')
xlabel('acceleration')
ylabel('SSIM')

subplot(223)
plot(1./accvector,squeeze(hfen(2,:,:)),'.-')
title('HFEN')
xlabel('acceleration')
ylabel('HFEN')

subplot(224)
plot(1./accvector,squeeze(PSNR(2,:,:)),'.-')
title('PSNR')
xlabel('acceleration')
ylabel('PSNR')

legend('sampling 1','sampling 2','sampling 3')
export_fig -native '1_errors.eps'
export_fig -native '1_errors.png'

figure(102)
subplot(221)
plot(noisevector,squeeze(MSE(2,:,:)).','.')
title('MSE')
xlabel('sampling pattern')
ylabel('MSE')

subplot(222)
plot(noisevector,squeeze(SSIM(2,:,:)).','.')
title('SSIM')
xlabel('sampling pattern')
ylabel('SSIM')

subplot(223)
plot(noisevector,squeeze(hfen(2,:,:)).','.')
title('HFEN')
xlabel('sampling pattern')
ylabel('HFEN')

subplot(224)
plot(noisevector,squeeze(PSNR(2,:,:)).','.')
title('PSNR')
xlabel('sampling pattern')
ylabel('PSNR')

export_fig -native '2_errors.eps'
export_fig -native '2_errors.png'

figure(103)
% QQ PLOT
for ii=1:length(accvector); 
for jj=1:length(noisevector);
    subplot(length(accvector),length(noisevector),(ii-1)*length(noisevector)+jj)
temp=real(squeeze(R{2,ii,jj})-ImRef);
qqplot(temp(:))
end
end
export_fig -native '3_QQ.eps'
export_fig -native '3_QQ.png'



I4=[]; I4temp=[];
for ii=1:length(accvector); 
for jj=1:length(noisevector);
I4temp=[I4temp,abs(R{2,ii,jj})./max(abs(R{2,ii,jj}(:)))];
end;
I4=[I4;I4temp];
clear I4temp; I4temp=[];
end;
figure(104);
imshow(I4,[])
export_fig -native '4_RECONS.eps'
export_fig -native '4_RECONS.png'

xv=120:200; yv=xv;
I5=[]; I5temp=[];
for ii=1:length(accvector); 
for jj=1:length(noisevector);
I5temp=[I5temp,abs(R{2,ii,jj}(xv,yv))./max(abs(R{2,ii,jj}(:)))];
end;
I5=[I5;I5temp];
clear I5temp; I5temp=[];
end;
figure(105);
imshow(I5,[])
export_fig -native '5_RECONS.eps'
export_fig -native '5_RECONS.png'


figure(106)
plot(noisevector,squeeze(SSIM(2,:,:)).','.-')
title('SSIM')
xlabel('1/noise')
ylabel('SSIM')
export_fig -native '6_SSIM.eps'
export_fig -native '6_SSIM.png'

figure(107)
imshow(abs(ImRef),[]);
export_fig -native '7_imRef.eps'
export_fig -native '7_imRef.png'

I5=[]; I5temp=[];
for ii=1:length(accvector); 
for jj=1:length(noisevector);
I5temp=[I5temp,KD{3,ii,jj}.M.*KD{3,ii,jj}.MNSA];
end;
I5=[I5;I5temp];
clear I5temp; I5temp=[]
end;
figure(108)
cmap=hot(200); cmap=cmap([1,20:200],:);
imshow(I5,[]);colormap(cmap)
export_fig -native '8_MASKS.eps'
export_fig -native '8_MASKS.png'

figure(109);
imshow(abs(ImRef(xv,yv)),[]); axis off
export_fig -native '9_ImRefZ.eps'
export_fig -native '9_ImRefZ.png'


