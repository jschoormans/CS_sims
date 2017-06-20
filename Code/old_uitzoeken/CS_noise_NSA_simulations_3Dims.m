clear all; close all; clc; 
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations')
addpath(genpath('/opt/amc/bart')); vars;
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab Collection/arrShow-develop')); %array show
addpath(genpath('/home/jschoormans/lood_storage/divi/projects/cosart/Matlab Collection/nifti'))
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations')

%VARS
folder='/home/jschoormans/lood_storage/divi/Temp/jasper/20151104_3DimSDE_ISMRM/v2/'
% file='20_22102015_1512550_9_2_wip3d05mm35mmacc4nsa4senseV4.raw'
file='20_22102015_1504222_6_2_wip3d07mm40mmfullsenseV4.raw'
MR=CS_Recon([folder,file])
MR.Parameter.Recon.ArrayCompression='yes'
MR.Parameter.Recon.ACNrVirtualChannels='1'
MR.Perform1; 
MR.SortData;


%MAKE SENS MAPS AND K_SPACE (SIMULATION DATA)

K=MR.Data;
[nx,ny,nz,nc]=size(K)
sens=ones(nx,ny,nz,nc); %temp for now;

mask=bart('poisson -Y206 -Z68 -y2 -z2 -C20 -v');
M=repmat(mask,[nx 1 1 nc]);
acceleration=sum(mask(:)==1)/((ny/2)*(nz/2)*(pi))

Ku=M.*double(K);      %undersampled k-space%
%CS
ImCS=bart('pics -RW:6:0:0.1 -i20 -d5',Ku,sens);
Imlin=ifftn(squeeze(Ku));

figure(1);
subplot(221);; subplot(222); imshow(squeeze(mask))
subplot(223); imshow(abs(squeeze(Imlin(:,:,50,1))),[]); subplot(224); imshow(abs(squeeze(ImCS(:,:,50,1))),[])
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
NoiseLevels=[1,10,50]
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
NoiseLevels=[1:20:101]
regs=[1e-3 5e-3 1e-2 5e-2 1e-1 5e-1]
[Recon, D]=analyse_NSreg(Im,K,M,sens,NoiseLevels,regs);

%%
figure(5);
subplot(121)
plot(squeeze(D(:,:,6)'));
subplot(122)
plot(squeeze(D(:,:,6)'));







