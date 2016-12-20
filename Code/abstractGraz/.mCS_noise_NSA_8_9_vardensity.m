% NEW SIMULATIONS FOR CS 
%%
% SCANS: 256/512/1024 2D SCANS OF PHANTOM WITH HOLES
% ALSO: 256x256x?? 3D SCAN
%GOAL: UNDERSAMPLE WITH VAR DENSITY (1 CHANNEL)
% EVALUATE PERFORMANCE

%GOAL II: UNDERSAMPLE AND EVALUATE 3D SCANS

% GOAL III: INFLUENCE OF CHANNELS ON EVALUATION

% NEXT GOAL: DO OVER WITH PHANTOM FRUIT SCANS/ANATOMICAL SCANS OF STATIC
% TISSUE (BRAIN --> PROSPECTIVE THO IN SIMILAR TIME FRAMES)


%%
clear all; close all; clc; 
disp('adding paths and setting up stuff...')
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/exportfig'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/tightfig/'))
addpath(genpath('/opt/amc/matlab/toolbox/MRecon-3.0.519'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/'))
addpath(genpath('/opt/amc/bart')); vars
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab'))
CC=clock; CCC=[num2str(CC(2)),'-',num2str(CC(3)),'-',num2str(CC(4)),'-',num2str(CC(5))];
figfolder=['/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Figures/',CCC]
dir=mkdir(figfolder);cd(figfolder)

%% defining parameters 
clear dir;
disp('setting parameters')
P=struct;
P.chan=1;
P.folder='/home/jschoormans/lood_storage/divi/Projects/cosart/scans/Graz1-phantoms_NSA_resolution/'
cd(P.folder)
files=dir('*.raw');
P.file=files(7).name;

cd(figfolder) %CD figfolder to save figures in good folder

%% PHANTOM RESOLUTION

MR=MRecon([P.folder,P.file]); 
MR.Parameter.Parameter2Read.typ=1;
disp(['picking channel',num2str(P.chan)])
allchans=MR.Parameter.Parameter2Read.chan;
MR.Parameter.Parameter2Read.chan=allchans(P.chan);
MR.Parameter.Recon.ImmediateAveraging='no'
MR.ReadData;
disp('random phase corr')
MR.RandomPhaseCorrection;
MR.PDACorrection;
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;
disp('sorting')
MR.SortData;

% MR.Average;
MR.RemoveOversampling;
% MR.RingingFilter;
MR.ZeroFill;
% MR.K2I;
% MR.ShowData
 disp('data to k-space')
K=MR.Data;
K=K./max(K(:)); %normalize kspace
sens=ones(size(K,1),size(K,2));
%% make ref image
disp('making ref image')
reg=0.02;
Kref=mean(K,12); %still use coils though
ImRef=fftshift(squeeze(bart(['pics -RW:7:0:',num2str(reg),' -S -e -i500'],Kref,sens)),2); %to compare image
ImRef=abs(ImRef);
ImRef=ImRef./(max(ImRef(:)));
disp('calculating sense maps')
sens=bart('ecalib -m1 -r20',Kref); %should be [nx,ny,1,nc] ?
disp('sense maps finsihed')

%% make undersampled k-spaces

disp('making the undersampled k-spaces')
accvector=[5]
% accvector=[1,4,8]
noisevector=[5,6,7] % in this case the different W 
for jj=1:length(noisevector)
for ii=1:length(accvector);
P.jjj=noisevector(jj); %weighting scheme (3=normal)
P.usedyns=1; %for example
P.acc=accvector(ii);
[KD{1,ii,jj}, KD{2,ii,jj}, KD{3,ii,jj}]=makeNoisyKspacefromdynamics(K,P);
end
end

%% CONJUGATE GRADIENT RECONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TVWeight = 0; 	% Weight for TV penalty
xfmWeight = 0.001;	% Weight for Transform L1 penalty
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
            param.xfmWeight=xfmWeight;
            im_dc2 = FT'*(param.data.*KD{3,ii,jj}.M./KD{3,ii,jj}.pdf); %linear recon; scale data to prevent low-pass filtering
            res = XFM*im_dc2;
            param.V=(KD{3,ii,jj}.MNSA.*KD{3,ii,jj}.M);
            param.Debug=0;
            param.xfmWeight=xfmWeight*(mean(param.V(KD{3,ii,jj}.M~=0)));
            
            for n=1:3
                res = fnlCg_test(res,param);
                R{2,ii,jj} = XFM'*res;
            end
            R{2,ii,jj}=fftshift(abs(R{2,ii,jj}),2);
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
        
        coordsx=665:690; coordsy=440
        Edge=abs(squeeze(R{2,ii,jj}(coordsx,coordsy)))
        Center=14 %whatever
        S(2,ii,jj)=SigmoidFitting(Edge,Center,1)
        
    end
end



%% Visualization code 
cd(figfolder)
figure(101)
set(gcf,'color','w');
subplot(221)
plot(1./accvector,squeeze(MSE(2,:,:)),'.-')
title('MSE')
xlabel('sampling fraction')
ylabel('MSE')

subplot(222)
plot(1./accvector,squeeze(SSIM(2,:,:)),'.-')
title('SSIM')
xlabel('sampling fraction')
ylabel('SSIM')

subplot(223)
plot(1./accvector,squeeze(hfen(2,:,:)),'.-')
title('HFEN')
xlabel('sampling fraction')
ylabel('HFEN')

subplot(224)
plot(1./accvector,squeeze(PSNR(2,:,:)),'.-')
title('PSNR')
xlabel('sampling fraction')
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
axis off
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

