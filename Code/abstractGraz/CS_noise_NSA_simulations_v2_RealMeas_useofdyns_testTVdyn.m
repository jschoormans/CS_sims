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
% file='ph_25072016_2030563_13_2_wipt1wffekiwi12810dynclearV4.raw' 
% file='ph_25072016_2034581_15_2_wipt1wffekiwi128100dynclearV4.raw' %actually 256 10
% file='ph_25072016_2036025_16_2_wipt1wffekiwi256100dynclearV4.raw'
% file='ph_25072016_2133129_28_2_wipt1wffekiwi51210dynclearV4.raw' %512 10 dyns
% file='ph_25072016_2125045_27_2_wipt1wffekiwi256100dynclearV4.raw'  
file='ph_25072016_2137350_29_2_wipt1wffekiwi51230dynclearV4.raw'
MR=MRecon([folder,file]); 
MR.Parameter.Parameter2Read.typ=1;
MR.Parameter.Parameter2Read.chan=2;
MR.Parameter.Recon.ImmediateAveraging='no'
MR.ReadData;
MR.RandomPhaseCorrection;
MR.PDACorrection;
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;
MR.SortData;
% MR.Average;
MR.RemoveOversampling;
MR.RingingFilter;
MR.ZeroFill;
% MR.K2I;
% MR.ShowData
K=MR.Data;

K=K./max(K(:)); %normalize kspace
sens=ones(size(K,1),size(K,2));
%% make ref image
reg=0.02;
Kref=mean(K,12); %still use coils though
ImRef=fftshift(squeeze(bart(['pics -RW:7:0:',num2str(reg),' -S -e -i100'],Kref,sens)),2); %to compare image
ImRef=abs(ImRef);
ImRef=ImRef./(max(ImRef(:)));
sens=bart('ecalib -m1 -r20',Kref); %should be [nx,ny,1,nc] ?

%% make k-spaces
accvector=[5,4,3,2,1];
noisevector=[1]
P.jjj=5;
P.usedyns=3; %for example
for jj=1:length(noisevector)
for ii=1:length(accvector);
P.acc=accvector(ii);
[KD{1,ii,jj}, KD{2,ii,jj}, KD{3,ii,jj}]=makeNoisyKspacefromdynamics_noav(K,P);
end
end

%%
reg=0.02;
for ii=1:length(accvector);
    for jj=1:length(noisevector)
        if KD{3,ii,jj}.jjj==3 || KD{3,ii,jj}.jjj==4 ||KD{3,ii,jj}.jjj==5
            R{2,ii,jj}=fftshift(bart(['pics -RW:7:0:',num2str(reg),',RT:1024:0:0.01 -d5 -S -e -i100'],KD{2,ii,jj},sens),2)
            R{2,ii,jj}=abs(R{2,ii,jj});
            R{2,ii,jj}=R{2,ii,jj}./max(R{2,ii,jj}(:)); %scaling!
        else
            
        end
    end
end

%% average dynamics

for ii=1:length(accvector);
    for jj=1:length(noisevector)
        R{2,ii,jj}=mean(R{2,ii,jj},10)
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
