%% KLADBLOK EXPERIMENT VNSA 51
clear all; close all; clc;

%% load k-space
if ispc()
    addpath(genpath('L:\basic\divi\Projects\cosart\CS_simulations'))
    addpath(genpath('L:\basic\divi\Projects\cosart\CS_simulations\experiments\VNSA_51_retro')) % for this: only the local code!!
    addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\Wavelab850'))
    addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\exportfig'))
    experimentfolder=['L:\basic\divi\Projects\cosart\CS_simulations\experiments\VNSA_51_retro\',date]
    mkdir(experimentfolder)
    addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\spot-master'))
    addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\Wavelab850\'))
    rmpath(genpath('L:\basic\divi\Projects\cosart\CS_simulations\experiments\VNSA_51_retro\ReconCode'))
    vars
    cd('L:\basic\divi\Ima\parrec\Jasper\VNSA\VNSA_51\VNSA_51')
else
    addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/exportfig'))
        addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/imagine'))
    addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/tightfig/'))
    addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/'))
    addpath(genpath('/opt/amc/bart')); vars
    addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/experiments/VNSA_51_retro' ))
    experimentfolder=['/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/experiments/VNSA_51_retro/',date,'2']
    mkdir(experimentfolder)
    
    cd('/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/VNSA/VNSA_51/VNSA_51')
end
files=dir('*.mat')
load(files(3).name)

 %% FIX CHECKERBOARD ISSUE
 %%%%%%%%%%%%%%%%%TEMP TEMP TEMP
 rr=mod([1:720],2); rr=double(rr)-double(rr==0);
 checkerboard=repmat(rr,[720 1]);
 checkerboard=repmat(checkerboard,[1 1 8 100]);
 K_ch=checkerboard.*K;
 %%%%%%%%%%%%%%%%%%

%% make sense maps
disp('to do: use this in recon- fix checkerboard here as well')
K_ch_crop=K_ch(183:538,183:538,:,:); 

Kref=mean(K_ch_crop,4); %still use coils though
Kref=permute(Kref,[4 1 2 3]);
sens=bart('ecalib -m1 -r20',Kref); %should be [nx,ny,1,nc] ?
imagine(sens)
%% Imref
reg=0.01;
ImRef=(squeeze(bart(['pics -RW:7:0:',num2str(reg),' -S -e -i100'],Kref,sens))); %to compare image
ImRef=abs(ImRef);
ImRef=ImRef./(max(ImRef(:)));
figure(1); imshow(squeeze(abs(ImRef)),[])

%% RECON HERE
%K should be [kx ky kz ncoils nNSA]
% FOR CERTAIN ACC FACTORS - LOOP OVERNIGHT
xfmWeight=2e-4;
PR=struct;
PR.outeriter=1;
PR.Itnlim=20;
PR.noNSAcorr=false;
PR.TVWeight=0;
PR.TGVfactor=0;
PR.xfmWeight=xfmWeight;
PR.reconslices=1;
PR.squareksp=false;
PR.resultsfolder=''
PR.sensemaps=squeeze(sens); 
PR.sensemapsprovided=1

PR.visualize_nlcg=0;
PR.debug_nlcg=0;
PR.VNSAlambdaCorrection=1;
PR.WeightedL2=1;
PR.Scaling=true;
PR.VNorm=1; %power of weighting mat  rix (number of NSA)^p; 

%%

NNSA=[6];
Nacc=6 ;
Nave=3;
% SSSIM(ii,jj,kk,nave)=zeros(Nacc,3,NNSA,Nave);
rng('default');
rng(1)

for kk=NNSA%:nNSA
    for jj=1:3  %:3
        for ii=1:Nacc  %:length(accvector);
            
            P.usedyns=kk; %for example
            P.acc=ii;
            P.jjj=jj+4 % 5 6 7
            P.noiselevel=0
            for nave=1:Nave
                [~,~, KD]=makeNoisyKspacefromdynamics(K_ch_crop,P);
                
                Ko=KD.Ku_N2;
                Ko=permute(Ko,[5 1 2 3 4]);
                parfor rr=[1:6] %regularizaiton multiplication factors...
                    R=reconVarNSA(Ko,PR,rr);
                    Recon(ii,jj,kk,nave,rr,:,:)=abs(R.recon);
%                     SSSIM(ii,jj,kk,nave,rr)=ssim(ImRef,abs(R.recon)) %
%                     can be done afterwards
                end
            end
        end
    end
end
%%
% figure(5);subplot(221);plot(mean(SSSIM(:,:,6,:),4));  title('NSA =1'); legend('equal','periphery','center')
kk=6
figure(5);
% errorbar(mean(SSSIM(:,:,kk,:),4),10*var(SSSIM(:,:,kk,:),[],4));
plot(max(SSSIM(:,:,kk,:),[],4))
title('NSA =1'); legend('equal','periphery','center')

%%
clear QQ Q
figure(6); 
for jj=1:3
    QQ=[]
    for ii=1:Nacc
        QQ=[QQ,squeeze(Recon(ii,jj,6,1,:,:))]; 
    end
    Q{jj}=QQ;
end
imshow(cat(1,Q{1},Q{2},Q{3}))


%%
regionx=100:250;
nave=1
for kk=NNSA%:nNSA
        for jj=1:3  %:3
            for ii=1:Nacc  %:length(accvector);
                
                MSE(ii,jj,kk)=sum(sum((squeeze(Recon(ii,jj,kk,nave,regionx,regionx))-ImRef(regionx,regionx)).^2))
                SSIM2(ii,jj,kk)=ssim(ImRef(regionx,regionx),abs(squeeze(Recon(ii,jj,kk,nave,regionx,regionx))))
            end
        end
end

figure(7); plot(MSE(:,:,6));  legend('equal','periphery','center')
figure(8); plot(SSIM2(:,:,6));legend('equal','periphery','center')




%%
    save([experimentfolder,'\R_',num2str(ii),'_',num2str(jj),'_',num2str(kk),'.mat'],'R')
