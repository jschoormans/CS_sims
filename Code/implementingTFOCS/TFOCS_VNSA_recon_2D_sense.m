
clear all; close all; clc; 
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/spot-master'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/Wavelab850'))
addpath /home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/TFOCS-1.4
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/imagine'))
addpath(genpath('/opt/amc/bart/')); vars; 

%%
cd('/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/VNSA/VNSA11')
disp('loading k-space...')
load('K_n_10022017_1149564_6_2_wipvnsaacc5v3senseV4_VC6')
K_orig=K; 

%% 1: simple 2D example; iFFT in z-direction, one NSA

sl=1

K=mean(K_orig,5); % only take first NSA for now (mean of 5 NSAs)
K=ifft(K,[],1); % FOR NOW: iFFT along read dim
K=squeeze(K(sl,:,:,:));

mask=K(:,:,1)~=0;


figure(1); imshow(abs([K(:,:,1)./max(K(:)),mask]),[0 1]); axis off; title('kspace and mask')

[n1,n2,ncoils]=size(K)
for i=1:ncoils
    Kc=K(:,:,i);
Ku(:,i)=Kc(mask);
end; clear Kc;
Ku=vec(Ku);

%% make operators 

mat = @(x) reshape(x,n1,n2,ncoils); %function that reshapes vector to matrix 
mat2d = @(x) reshape(x,n1,n2*ncoils); % functions that reshapes in 2d matrix for visualization purposes 
matcc =@(x) reshape(x,n1,n2); %function that reshapes vector to matrix 

AA=opDFT2(n1,n2,1) %DFT2 MATRIX

l=[1:n1*n2];
mask2=(l(~mask(:))); %indices of rows in FDT2 matrix that should be removed (were not sampled!)
R=opExcise(AA,mask2,'rows')
E1=opBlockDiag(R,R,R,R,R,R)


%% perform linear recon 

linear_recon=E1'*Ku;

sensmaps=bart('ecalib -r20 -m1 -S -k5',permute(K,[4 1 2 3]));  %% make sense maps 
sensmaps=fftshift(fftshift(sensmaps,3),2);
% note: SENSE MAPS LOOK VERY BAD (LOW SNR/ TOO SMALL KERNEL)?


figure(2)
imshow(abs([mat2d(linear_recon./max(linear_recon(:)));mat2d(sensmaps./max(sensmaps(:)))]),[]); % linear recon for all coils apart 


%% make sense maps operators

clear S
for i=1:ncoils; S{i}=opDiag(conj(sensmaps(1,:,:,i))); end; 
S2= horzcat(S{i},S{2},S{3},S{4},S{5},S{6})
E2=E1*S2'

%linear recon s
linear_recon_s=E2'*Ku;
figure(3); imshow(abs(matcc(linear_recon_s)),[]); axis off; 

%% Call the TFOCS solver

mu              = 5;
er              = @(x) norm(x(:)-linear_recon_s(:))/norm(linear_recon_s(:)); %NOT SURE ABOUT THIS: WE DONT KNOW SOLUTION YET...
opts = [];
opts.errFcn     = @(f,dual,primal) er(primal);
opts.maxIts     = 200;
opts.printEvery = 20;
opts.tol        = 1e-8;
opts.stopcrit   = 4;
opts.alg='GRA'  

x0 = linear_recon_s;  %first guess
z0  = [];           % we don't have a good guess for the dual (is this true)


a=0;
% build operators:
% A           = linop_handles([N,N], @(x)fft2(x), @(x)ifft2(x),'c2c');
A=linop_spot(E2,'c2c')
normA2      = 1;
W_wavelet   = linop_spot(opWavelet2(n1,n2,'daubechies'),'c2c');
normWavelet      = linop_normest( W_wavelet );

contOpts            = [];
contOpts.maxIts     = 4;


EPS=2e5
[x_wavelets,out_wave] = solver_sBPDN_W( A, W_wavelet, double(Ku), EPS, mu, ...
    x0(:), z0, opts, contOpts);


 %possible reweighting here
 
 figure(3); imshow(abs([matcc(x_wavelets./max(x_wavelets(:))),matcc(linear_recon_s./max(linear_recon_s(:)))]),[])
 
%% LASSO 
%{
EPS=1e-3

er          = @(x) norm(x-x0)/norm(x0);    % AND THIS?
obj_ref = EPS*norm(x0,1) + sum( (E2*x0 - double(Ku)).^2)/2; %WHAT IS THIS EXACTLY?
    
    
opts = struct('restart',-Inf,'tol',1e-13,'maxits',1000);
opts.errFcn     = { @(f,primal) er(primal), ...
                    @(f,primal) f - obj_ref   }; 
                
[ x_wavelets,out_wave ] = solver_L1RLS( A, double(Ku), EPS, x0, opts );
     
                
                
% VISUALIZE
[x_wavelets,out_wave]= solver_L1RLS( A, double(Ku),EPS,x0, opts ); % THIS IS THE ONE, BUT X SHOULD BE IN SPARSIFYING FORM (WAVELET)

%}
