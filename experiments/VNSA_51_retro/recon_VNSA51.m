%% KLADBLOK EXPERIMENT VNSA 51
clear all; close all; clc; 
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/exportfig'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/tightfig/'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/'))
addpath(genpath('/opt/amc/bart')); vars
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab'))

CC=clock; CCC=[num2str(CC(2)),'-',num2str(CC(3)),'-',num2str(CC(4)),'-',num2str(CC(5))];
figfolder=['/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/experiments/VNSA_51_retro/Figures/',CCC]
DD=mkdir(figfolder);

%% load k-space
cd('/scratch/jschoormans/VNSA_51/')
files=dir('*.mat')
load(files(2).name)

%% make sense maps
Kref=mean(K,4); %still use coils though
Kref=permute(Kref,[4 1 2 3]);

sens=bart('ecalib -m1 -r20',Kref); %should be [nx,ny,1,nc] ?
%% Imref
reg=0;
ImRef=fftshift(squeeze(bart(['pics -RW:7:0:',num2str(reg),' -S -e -i100'],Kref,sens)),2); %to compare image
ImRef=abs(ImRef);
ImRef=ImRef./(max(ImRef(:)));
figure(1); imshow(squeeze(abs(ImRef)),[])
%% make k-spaces
accvector=[3];
P.jjj=5; %weighting scheme (3=normal)
P.usedyns=10; %for example
for jj=1:3
for ii=1:length(accvector);
P.usedyns=3; %for example
P.acc=accvector(ii);
P.jj=jj
P.noiselevel=0
[KD{1,ii,jj}, KD{2,ii,jj}, KD{3,ii,jj}]=makeNoisyKspacefromdynamics(K,P);
end
end

%% RECON HERE
%K should be [kx ky kz ncoils nNSA]

PR=struct;
PR.outeriter=4;
PR.Itnlim=4;
PR.noNSAcorr=false;
PR.TVWeight=(6e-3);
PR.TGVfactor=0;
PR.xfmWeight=0;
PR.reconslices=1;
PR.squareksp=false;
PR.resultsfolder='/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/experiments/VNSA_50_retrospective_new/Results/'
for jj=1:3
for ii=1:length(accvector);
Ko=KD{3,ii,jj}.Ku_N2;
Ko=permute(K,[5 1 2 3 4]);
R{2,ii,jj}=reconVarNSA(Ko,PR)
end
end
%% Estimate errors
clear MSE SSIM hfen
for ii=1:length(accvector);
    for jj=1:3
        recon=ifftshift(R{2,ii,jj}.recon,1);
        
        MSE(2,ii,jj)=mean(mean(abs(squeeze(recon)-ImRef).^2,1),2)
        SSIM(2,ii,jj)=ssim(abs(squeeze(recon)),abs(ImRef))
        hfen(2,ii,jj)=HFEN(abs(squeeze(recon)),abs(ImRef))
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

xv=170:250; yv=xv;
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