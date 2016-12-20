%%
clear all; close all; clc; 
disp('adding paths and setting up stuff...')

CC=clock; CCC=[num2str(CC(2)),'-',num2str(CC(3)),'-',num2str(CC(4)),'-',num2str(CC(5))];
figfolder=['L:\basic\divi\Projects\cosart\CS_simulations\Figures\',CCC]
figdir=mkdir(figfolder);

%% LOAD KSPACE

cd('L:\basic\divi\Projects\cosart\scans\CS-3D-prospective-grapefruit')

filenumber=6
filename=['K',num2str(filenumber),'.mat']
load(filename)

%% FFT IN MEASUREMENT DIRECTIONS
K=ifft(K,[],1);
K=K./max(K(:)); %normalize kspace
K=squeeze(K);
%% RECON
sl=50;
for ch=1:size(K,4)
    ch
%% make mask
Ks=squeeze(K(sl,:,:,ch,:));
data=mean(Ks,3);   %mean of data over NSA
data=squareksp(data); %make square and power of n 
mask=double(data~=0); %2d mask
fullmask=Ks~=0;
MNSA=sum(fullmask,3); %sum over NSA
MNSA=squareksp(MNSA);
pdf=ones(size(data)); %could be fixed for faster convergence!!!

%% CONJUGATE GRADIENT RECONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TVWeight = 0; 	% Weight for TV penalty
xfmWeight = 0.0001;	% Weight for Transform L1 penalty
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
N=[size(data,1) size(data,2)]


            %generate Fourier sampling operator
            param.data=data;
            FT = p2DFT(mask, N, 1, 2);
            param.FT = FT;
            param.xfmWeight=xfmWeight;
            im_dc2 = FT'*(param.data.*mask./pdf); %linear recon; scale data to prevent low-pass filtering
            res = XFM*im_dc2;
            param.V=(MNSA.*mask);
            param.Debug=0;
            param.xfmWeight=xfmWeight*(mean(param.V(mask~=0)));
            
            for n=1:6
                res = fnlCg_test(res,param);
            end
            recon = XFM'*res;
            
               recon=fftshift(fftshift(abs(recon),2),1);
               recon=recon./max(recon(:));
            
            R(:,:,ch)=recon;
end

%% COMBINE CHANNELS
Rfull=squeeze(sqrt(sum(R.^2,3)));

save(['R-sl50-',num2str(filenumber),'.mat'],'Rfull')
figure(1); imshow(Rfull,[])
%%

cd(figfolder)
