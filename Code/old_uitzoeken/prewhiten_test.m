% prewhiten data test; 
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations')
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/MRIPhantomv0-8'))
clear all; close all; clc;



%simulate k-space
[K]=genPhantomKspace(512,1);
ny=size(K,2);
nz=size(K,1);
sens=ones(ny,nz); %change to sense based on phantoms...

%%

acc=5
%sampling pattern
[pdf,val] = genPDF(size(K),5,1/acc,2,0,0);
Mfull=genEllipse(size(K,1),size(K,2));
M=genSampling(pdf,10,100).*Mfull;



MNSA=ceil(1./pdf);
% MNSA=ceil(pdf*acc);

MNSA=acc*ones(size(MNSA));
%% add noise to kspace
clear Ku_N2
NoiseLevel=4e-4;
for iii=1:max(MNSA(:)) %Matrix of NSA values
K_N=addNoise(K,NoiseLevel);
Ku_N1=squeeze(M(:,:).*(MNSA(:,:)>=iii)).*K_N;
Ku_N2(1,:,:,1,1,iii)=permute(Ku_N1,[3,2,1,4,5]);
end
Ku_Nvar1=sum(Ku_N2,6)./permute(MNSA,[3 2 1]);
disp('check dimensionssss')

sum(MNSA(:).*M(:))
sum(Mfull(:))
%% RECONS 
% ORDINARY RECON
reg=0.05
R1=bart(['pics -RW:7:0:',num2str(reg),' -S -e -i100'],K_N.*Mfull,sens);

%% RECON R2: WITHOUT PREAVERAGING

kspacefull=[]; trajfull=[];
for ii=1:size(Ku_N2,6)
    
[kspace, traj]=calctrajBART(squeeze(Ku_N2(1,:,:,1,1,ii))); 

kspacefull=[kspacefull;kspace];
trajfull=[trajfull,traj];
end

%
kspace=kspacefull.';
traj=trajfull.*(ny/2);
R2=bart(['pics -RW:7:0:',num2str(reg),' -S  -e -i100 -t'],traj,kspace,sens);

%% RECON R3: same traj, with preaveraging

 
[kspace, traj]=calctrajBART(squeeze(Ku_Nvar1)); 
traj=traj.*(ny/2);
R3=bart(['pics -RW:7:0:',num2str(reg),' -S  -e -i100 -t'],traj,kspace.',sens);



%%
% VISUALIZATION

figure(1)
subplot(131)
imshow(abs(squeeze(R1)),[])
subplot(132)
imshow(abs(squeeze(R2)),[])
subplot(133)
imshow(abs(squeeze(R3)),[])

figure(3)
hold on
% plot(abs(squeeze(R1(:,100))),'r')
plot(abs(squeeze(R2(250,:))),'k')
plot(abs(squeeze(R3(250,:))),'y')
hold off
