clear all; close all; clc;
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Code')
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/MRIPhantomv0-8'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/tightfig'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/WaveLab850'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/sparseMRI_v0.2'))


folder='/home/jschoormans/lood_storage/divi/Projects/cosart/scans/Reconstructions/20160710_CS/'
file='20_11072016_1724136_13_2_wipt1wffeclearV4.raw'  
MR=MRecon([folder,file])
%%
MR.Parameter.Recon.ImmediateAveraging='No' 
MR.Parameter.Parameter2Read.typ=1
MR.Parameter.Parameter2Read.chan=[3]'

MR.ReadData;
MR.DcOffsetCorrection;
MR.PDACorrection;
MR.RandomPhaseCorrection;
MR.MeasPhaseCorrection;
MR.SortData;
MR.RemoveOversampling;

%% sample from data
acc=6;
ksp=MR.Data;
[ksp] = squareksp(ksp);
ksp=padarray(ksp,[0,256],'both');

[pdf,val] = genPDF(size(ksp(:,:,1)),5,1/acc,2,0,0);
Mfull=genEllipse(size(ksp,1),size(ksp,2));
% Mfull=repmat(Mfull,[size(MR.Data,1) 1 1]); %add coils
M=genSampling(pdf,10,100).*Mfull;

jjj=2%option to try different averaging approaches
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

end

MNSA=10*MNSA;
for iii=1:max(MNSA(:)) %Matrix of NSA values
Ku_N1=squeeze(M.*(MNSA(:,:)>=iii)).*ksp(:,:,:,1,1,1,1,1,1,1,1,iii);
Ku_N2(:,:,iii)=Ku_N1;
end
Ku_Nvar1=sum(Ku_N2,3)./MNSA;

%%
% Ku_Nvar1_FT=fftshift(fft(Ku_Nvar1,[],1));
Ku_Nvar1_FT=Ku_Nvar1;
%% Recons
param = init;

param.data=squeeze(Ku_Nvar1_FT)
N=size(param.data)
% M1=squeeze(M(1,:,:)); 
M1=M;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TVWeight = 0; 	% Weight for TV penalty
xfmWeight = 0.01;	% Weight for Transform L1 penalty
Itnlim = 8;		% Number of iterations

%generate Fourier sampling operator
FT = p2DFT(M, N, 1, 2);

%generate transform operator
XFM = Wavelet('Daubechies',4,4);	% Wavelet

% initialize Parameters for reconstruction

param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.TVWeight =0;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;
xfmWeight = 2;	% Weight for Transform L1 penalty

param.xfmWeight=xfmWeight
im_dc2 = FT'*(param.data.*M1./pdf); %linear recon; scale data to prevent low-pass filtering
res_orig = XFM*im_dc2;
for n=1:4
    res_orig = fnlCg(res_orig,param);
	R4{jjj} = XFM'*res_orig;
end
R4{jjj}=R4{jjj}./max(R4{jjj}(:));

%% RECON R5
res = XFM*im_dc2;
param.V=(MNSA);
param.Debug=0;
param.xfmWeight=xfmWeight*(mean(param.V(M~=0)))

for n=1:4
	res = fnlCg_test(res,param);
	R5{jjj} = XFM'*res;
end
R5{jjj}=R5{jjj}./max(R5{jjj}(:));

%% VISUALIZATION 
figure(100);
imshow(abs(cat(2,R4{jjj},R5{jjj})),[])





