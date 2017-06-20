% VARIABLE NSA PHANTOM RECON FILE
%COMBINES CODE FROM PREWHITEN_TEST_VAR_NSA AND
%DEMO_BRAIN_TEST_VARIANCE_MATRIX

%GOAL: COMPARE DIFFERENT APPROACHES OF RECONSTRUCTING A CS -KSPACE WITH A
%NON-UNIFORM NOISE PROFILE; OBTAINED BY VARING THE NUMBER OF SIGNAL
%AVERAGES OVER SPACE

cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Code')
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/MRIPhantomv0-8'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/tightfig'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/WaveLab850'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/sparseMRI_v0.2'))

clear all; close all; clc;


%% simulate k-space
phantom=1; %toggle simulation/loading of data
if phantom==1
    ny=512;nz=512;
    N=[ny nz]
    [K,sens]=genPhantomKspace(N(1),1);
    sens=permute(sens,[1 2 4 3]);
        nc=size(K,4);

else
    run exp_4_7_pomegranate.m
    K=K(:,:,:,1,:,:,:,:,:,:,:,:);  
    K=padarray(K,[56 56],'both');
    ny=size(K,1);
    nz=size(K,2); N=[ny,nz]
    nc=size(K,4);
    sens=ones(ny,nz,nc); %temptemptemp
    
end


acc=5

%sampling pattern
[pdf,val] = genPDF(N,5,1/acc,2,0,0);
Mfull=genEllipse(ny,nz);
Mfull=repmat(Mfull,[1 1 nc]); %add coils
M=repmat(genSampling(pdf,10,100),[1 1 nc]).*Mfull;
%% PARAMETER USED FOR CONJUGATE GRADIENT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TVWeight = 0; 	% Weight for TV penalty
xfmWeight = 0.00001;	% Weight for Transform L1 penalty
Itnlim = 8;		% Number of iterations

%generate Fourier sampling operator
FT = p2DFT(M, N, 1, 2);

%generate transform operator
XFM = Wavelet('Daubechies',4,4);	% Wavelet

% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.TVWeight =0;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;


%% SIGNAL AVERAGING APPROACHES
jjj=1 %option to try different averaging approaches
if jjj==1
    MNSA=ceil(1./pdf);
elseif jjj==2
    MNSA=ceil(pdf*acc);
elseif jjj==3
    MNSA=acc*ones(size(pdf));
elseif jjj==4 %extreme center-heavy
    MNSA=ones(size(pdf));
    MNSA(1+ny/4:3*ny/4,1+nz/4:3*nz/4)=100*ones(ny/2,nz/2);
elseif jjj==5 %extreme outside-heavy
    MNSA=100*ones(size(pdf));
    MNSA(1+ny/4:3*ny/4,1+nz/4:3*nz/4)=ones(ny/2,nz/2);
elseif jjj==6
        MNSA=min(ceil(1./pdf).^2,50);
elseif jjj==7
        MNSA=ceil(pdf*acc).^2;
end


%% add noise to kspace
clear Ku_N2

if phantom==1
NoiseLevel=8e-4;
for iii=1:max(MNSA(:)) %Matrix of NSA values
K_N=addNoise(K,NoiseLevel);
Ku_N1=repmat(squeeze(M(:,:).*(MNSA(:,:)>=iii)),[1 1 size(K,3)]).*K_N;
Ku_N2(1,:,:,:,1,iii)=permute(Ku_N1,[1,2,3,4]);
end
else 
    for iii=1:max(MNSA(:)) %Matrix of NSA values
    K_N=K(:,:,:,1,:,:,:,:,:,:,iii);
    Ku_N1=repmat(squeeze(M(:,:).*(MNSA(:,:)>=iii)),[1 1 size(K,3)]).*K_N;
    Ku_N2(1,:,:,:,1,iii)=permute(Ku_N1,[1,2,3,4]);
    end
    end

Ku_Nvar1=sum(Ku_N2,6)./permute(repmat(MNSA,[1 1 nc 1]),[4 1 2 3]);

% sum(MNSA(:).*M(:)) %number of sampling points
% sum(Mfull(:)) %number of sampling points
data=squeeze(Ku_Nvar1); %noisy 2D kspace data

param.data = data./max(abs(data(:)));  %give data to parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  RECONS 
%% ORDINARY RECON
reg=0.05
R1{jjj}=bart(['pics -RW:7:0:',num2str(reg),' -S -e -i20 -d5'],K_N.*Mfull,sens(end:-1:1,end:-1:1,:,:));
R1{jjj}=R1{jjj}./max(R1{jjj}(:));
%% RECON R2: WITHOUT PREAVERAGING
clear traj2;
ADMMreg=0.5
[kspace, traj]=calctrajBART(permute(Ku_N2,[1 2 3 4 6 5])); 
traj2(1,1,:)=traj(3,1,:); traj2(2,1,:)=traj(2,1,:); traj2(3,1,:)=traj(1,1,:); %FOR 2D signals; when we do not want ANY frequency encoding!
R2{jjj}=bart(['pics -RW:7:0:0.02 -S -u',num2str(ADMMreg),' -m -i30 -d5 -t'],traj2,kspace,sens);
R2{jjj}=R2{jjj}./max(R2{jjj}(:));

%% RECON R3: same traj, with preaveraging
clear traj2;
[kspace, traj]=calctrajBART((Ku_Nvar1)); 
traj2(1,1,:)=traj(3,1,:); traj2(2,1,:)=traj(2,1,:); traj2(3,1,:)=traj(1,1,:); %FOR 2D signals; when we do not want ANY frequency encoding!
R3{jjj}=bart(['pics -RW:7:0:0.02 -S -u',num2str(ADMMreg./mean(MNSA(M~=0))),' -m -i30 -d5 -t'],traj2,kspace,sens);
R3{jjj}=R3{jjj}./max(R3{jjj}(:));

%% RECON R4: conjugate gradient method - with pre-averaging but no variance matrix

xfmWeight = 0.001*25;	% Weight for Transform L1 penalty

param.xfmWeight=xfmWeight
im_dc2 = FT'*(param.data.*M./pdf); %linear recon; scale data to prevent low-pass filtering
res_orig = XFM*im_dc2;
for n=1:6
    res_orig = fnlCg(res_orig,param);
	R4{jjj} = XFM'*res_orig;
end
R4{jjj}=R4{jjj}./max(R4{jjj}(:));

%% RECON R5
res = XFM*im_dc2;
param.V=sqrt(MNSA);
param.Debug=0;
param.xfmWeight=xfmWeight*(mean(param.V(M~=0)))

for n=1:6
	res = fnlCg_test(res,param);
	R5{jjj} = XFM'*res;
end
R5{jjj}=R5{jjj}./max(R5{jjj}(:));



%% RECONS R6 and R7
pattern=param.V;
R6{jjj}=bart(['pics -RW:7:0:0.02 -S -i20 -d5 -p'],pattern,squeeze(Ku_Nvar1),sens); %all combinations of those do not seem to work...
R7{jjj}=bart(['pics -RW:7:0:0.02 -S -u',num2str(ADMMreg),'-i30 -d5 -m'],squeeze(Ku_Nvar1),sens);
R6{jjj}=R6{jjj}./max(R6{jjj}(:));R6{jjj}=bart(['pics -RW:7:0:0.02 -S -i20 -d5 -p'],squeeze(Ku_Nvar1),sens); %all combinations of those do not seem to work...

R7{jjj}=R7{jjj}./max(R7{jjj}(:));
figure(300)
imshow(abs(cat(2,R6{jjj},R7{jjj})),[])





%% VISUALIZATION 

figure(100);
imshow(abs(cat(2,R1{jjj}./max(R1{jjj}(:)),R2{jjj},R3{jjj},R4{jjj},R5{jjj})),[])

figure(101);
aax=200; bby=200;
imshow(abs(cat(2,R1{jjj}(aax:aax+50,bby:bby+50)./max(R1{jjj}(:)),R2{jjj}(aax:aax+50,bby:bby+50),R3{jjj}(aax:aax+50,bby:bby+50),R4{jjj}(aax:aax+50,bby:bby+50),R5{jjj}(aax:aax+50,bby:bby+50))),[])


%% SHOW NOISE PROFILE OF RECONS

No{1,jjj}=fft2c(R1{jjj}).*M-param.data;
No{2,jjj}=fft2c(R2{jjj}).*M-param.data;
No{3,jjj}=fft2c(R3{jjj}).*M-param.data;
No{4,jjj}=fft2c(R4{jjj}).*M-param.data;
No{5,jjj}=fft2c(R5{jjj}).*M-param.data;

%calculate noise stds for certain regions
for i=1:ny
    for j=1:nz
    Distance(i,j)=sqrt(abs(i-(ny/2 +1)).^2+abs(j-(nz/2+1)).^2);
    end
end

[distances,indexdist]=sort(Distance(:),'ascend');
indexdist=indexdist(M(:)>0);
distances=distances(M(:)>0);
figure(200)
hold on
plot(distances,mean((No{1,jjj}(indexdist)).^2)*(MNSA(indexdist)),'k--') %added noise profile
plot(distances,smooth((No{1,jjj}(indexdist)).^2,10000),'k')
plot(distances,smooth((No{2,jjj}(indexdist)).^2,10000),'b')
plot(distances,smooth((No{3,jjj}(indexdist)).^2,10000),'r')
plot(distances,smooth((No{4,jjj}(indexdist)).^2,10000),'y')
plot(distances,smooth((No{5,jjj}(indexdist)).^2,10000),'g')
hold off






