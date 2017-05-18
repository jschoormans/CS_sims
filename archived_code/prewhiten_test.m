% prewhiten data test; 
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations')
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/MRIPhantomv0-8'))

%simulate k-space
[K]=genPhantomKspace(256,1)
ny=size(K,2)
nz=size(K,1)
sens=ones(ny,nz); %change to sense based on phantoms...

%sampling pattern
[pdf,val] = genPDF(size(K),5,1/3,2,0,0);
Mfull=genEllipse(size(K,1),size(K,2));
M=genSampling(pdf,10,100).*Mfull;
MNSA=1./pdf;

%% add noise to kspace
NoiseLevel=2e-4;
for iii=1:max(MNSA(:,:)) %Matrix of NSA values
K_N=addNoise(K,NoiseLevel);
Ku_N1=squeeze(M(:,:).*(MNSA(:,:)>=iii)).*K_N;
Ku_N2(1,:,:,1,1,iii)=permute(Ku_N1,[3,2,1,4,5]);
end
Ku_Nvar1=sum(Ku_N2,6)./permute(MNSA(:,:),[3 2 1]); %mean of measures....

%% RECONS 
% ORDINARY RECON
reg=0.05
R1=bart(['pics -RW:7:0:',num2str(reg),' -S -e -i100'],Ku_Nvar1,sens);

%% RECON R2: WITHOUT PREAVERAGING
kspacefull=[]; trajfull=[];
for ii=1:3
    
[kspace, traj]=calctrajBART(squeeze(Ku_N2(1,:,:,1,1,ii))); 

kspacefull=[kspacefull;kspace];
trajfull=[trajfull,traj];
end

%
kspace=kspacefull.';
traj=trajfull.*128;
R2=bart(['pics -RW:7:0:',num2str(reg),' -S  -e -i100 -t'],traj,kspace,sens);

%% RECON R3: same traj, with preaveraging

 
[kspace, traj]=calctrajBART(squeeze(Ku_Nvar1)); 
traj=traj.*128;
R3=bart(['pics -RW:7:0:',num2str(reg),' -S  -e -i100 -t'],traj,kspace.',sens);



%%
% VISUALIZATION

figure(1)
subplot(131)
imshow(abs(squeeze(R1)).',[])
subplot(132)
imshow(abs(squeeze(R2)),[])
subplot(133)
imshow(abs(squeeze(R3)),[])



%% RECON WITH MRICS

% C=cov(K_N(:))
% V=sqrt(MNSA(:)).^(-1).*C; % variances;
% COV=diag(V)
% [S,D,notused] = svd(COV);
% inverse(V)

% Kuwr=Ku_Nvar1
% R2=bart(['pics -RW:7:0:',num2str(reg),' -S -e -i100'],Kuwr,sens);
%% MODIFY MRICS




