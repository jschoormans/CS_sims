%% KLADBLOK EXPERIMENT VNSA 50
clear all; close all; clc; 

addpath(genpath('L:\basic\divi\Projects\cosart\CS_simulations'))
addpath(genpath('C:\Users\jschoormans\Dropbox\phD\bart-0.3.01')); vars 

cd('L:\basic\divi\Ima\parrec\Jasper\VNSA\VNSA_50\VNSA_50')
files=dir('*.mat')
load(files(5).name)

%% make sense maps
Kref=mean(K,4); %still use coils though
Kref=permute(Kref,[4 1 2 3]);

sens=bart('ecalib -m1 -r20',Kref); %should be [nx,ny,1,nc] ?
%% Imref
reg=0.001;
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

%{
%% PARAMETER USED FOR CONJUGATE GRADIENT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TVWeight = 0; 	% Weight for TV penalty
xfmWeight = 0.00001;	% Weight for Transform L1 penalty
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
xfmWeight = 0.001*25;	% Weight for Transform L1 penalty
N=[size(K,1) size(K,2)]

for ii=1:length(accvector);
    for jj=1:length(noisevector)
        if KD{3,ii,jj}.jjj==3 || KD{3,ii,jj}.jjj==4 ||KD{3,ii,jj}.jjj==5
            reg=noisevector(jj)
            R{2,ii,jj}=fftshift(bart(['pics -RW:7:0:',num2str(reg),' -d5 -S -e -i100'],KD{2,ii,jj},sens),2)
            R{2,ii,jj}=abs(R{2,ii,jj});
            R{2,ii,jj}=R{2,ii,jj}./max(R{2,ii,jj}(:)); %scaling!
        else
            %generate Fourier sampling operator
            param.data=KD{2,ii,jj};
            FT = p2DFT(KD{3,ii,jj}.M, N, 1, 2);
            param.FT = FT;
            param.xfmWeight=xfmWeight
            im_dc2 = FT'*(param.data.*KD{3,ii,jj}.M./KD{3,ii,jj}.pdf); %linear recon; scale data to prevent low-pass filtering
            res = XFM*im_dc2;
            param.V=sqrt(KD{3,ii,jj}.MNSA);
            param.Debug=0;
            param.xfmWeight=xfmWeight*(mean(param.V(KD{3,ii,jj}.M~=0)))
            
            for n=1:6
                res = fnlCg_test(res,param);
                R{2,ii,jj} = XFM'*res;
            end
            R{2,ii,jj}=fftshift(abs(R{2,ii,jj}),2)
            R{2,ii,jj}=R{2,ii,jj}./max(R{2,ii,jj}(:));
        end
    end
end
%}
%% RECON HERE
%K should be [kx ky kz ncoils nNSA]

PR=struct;
PR.outeriter=4;
PR.Itnlim=15;
PR.noNSAcorr=false;
PR.TVWeight=(0);
PR.TGVfactor=0;
PR.xfmWeight=0%.3e-2;
PR.reconslices=1;
PR.squareksp=true;
PR.resultsfolder='/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/experiments/VNSA_50_retrospective_new/Results/'
PR.filename=files(5).name; disp('fix!!!')
PR.sensemapsprovided=0;
for jj=1:3
for ii=1:length(accvector);
K=KD{3,ii,jj}.Ku_N2;
K=permute(K,[5 1 2 3 4]);
PR=reconVarNSA(K,PR)
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