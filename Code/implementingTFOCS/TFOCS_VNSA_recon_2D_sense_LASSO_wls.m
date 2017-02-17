
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/spot-master'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/Wavelab850'))
addpath /home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/TFOCS-1.4
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/imagine'))
addpath(genpath('/opt/amc/bart/')); vars; 
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/TFOCS-1.4/examples'))
%%
clear all; close all; clc; 
disp('loading k-space...')
% load('K_scan13_nocc')
cd('/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/VNSA/VNSA12')
disp('loading k-space...')
% load('K_n_12022017_1411245_2_2_wipvnsaprescanfullfa2senseV4_noVC')
load('K_n_12022017_1416122_3_2_wipvnsafa2acc5v1senseV4_noVC')
% load('K_n_12022017_1420335_4_2_wipvnsafa2acc5v3senseV4_noVC')
K_orig=K; 

%% 1: simple 2D example; iFFT in z-direction, one NSA

sl=1

K=mean(K_orig,5); % only take first NSA for now (mean of 5 NSAs)
K=ifft(K,[],1); % FOR NOW: iFFT along read dim
K=squeeze(K(sl,:,:,:));

mask=K(:,:,1)~=0;


figure(1); imshow(abs([K(:,:,1)./max(K(:)),mask]),[0 1]); axis off; title('kspace and mask')

clear Ku
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
E1=opBlockDiag(R,R,R,R,R,R,R,R,R,R,R,R,R) % can this be done easier?

%% perform linear recon 

linear_recon=E1'*Ku;

sensmaps=bart('ecalib -r20 -m1 -k5',permute(K,[4 1 2 3]));  %% make sense maps 
sensmaps=fftshift(fftshift(sensmaps,3),2);
% note: SENSE MAPS LOOK VERY BAD (LOW SNR/ TOO SMALL KERNEL)?


figure(2)
imshow(abs([mat2d(linear_recon./max(linear_recon(:)));mat2d(sensmaps./max(sensmaps(:)))]),[]); % linear recon for all coils apart 


% make sense maps operators  
clear S
for i=1:ncoils; S{i}=opDiag(conj(sensmaps(1,:,:,i))); end; 
S2= horzcat(S{i},S{2},S{3},S{4},S{5},S{6},S{7},S{8},S{9},S{10},S{11},S{12},S{13})

E2=E1*S2'


pdf=estPDF(double(mask)); %remove lowpass filter effect 
FiltOp=opDiag(1./pdf(mask==1));
FiltOp=opBlockDiag(FiltOp,FiltOp,FiltOp,FiltOp,FiltOp,FiltOp,FiltOp,FiltOp,FiltOp,FiltOp,FiltOp,FiltOp,FiltOp)

linear_recon_s=E2'*FiltOp*Ku;
figure(3); imshow(abs(matcc(linear_recon_s)),[]); axis off; 

%% Wavelet operator

W=opWavelet2(n1,n2,'daubechies')
E3=E2*W'
%%  scaling operator 
scaling =false
if scaling==true
MNSA=K_orig(100,:,:,1,:); MNSA=sum(squeeze(MNSA)~=0,3);
else
    MNSA=ones(n1,n2);
end

figure(99); imshow(MNSA,[]);
ScalOp=opDiag(vec(MNSA(MNSA>0)));
ScalOp=opBlockDiag(ScalOp,ScalOp,ScalOp,ScalOp,ScalOp,ScalOp,ScalOp,ScalOp,ScalOp,ScalOp,ScalOp,ScalOp,ScalOp);
E4=ScalOp*E3;

Ksc= ScalOp*double(Ku);% scale data; 

%% Call the TFOCS solver
%% LASSO 
EPS=9e0



x0 = W*linear_recon_s;  %first guess



A=linop_spot(E3,'c2c')

% normA2      = 1;

er          = @(x) norm(x-x0)/norm(x0);    % AND THIS?
obj_ref = EPS*norm(x0,1) + sum(abs(E4*x0 - double(Ku)).^2)/2; %WHAT IS THIS EXACTLY?
    
opts = struct('restart',-Inf,'tol',1e-13,'maxits',50);
% opts.errFcn     = { @(f,primal) er(primal), ...
%                     @(f,primal) f - obj_ref   };  %WHAT IS THIS?? IT IS
%                     OPTIONAL...
opts.printEvery=10;
opts.nonneg=false
opts.alg='TS' %TS CG GRA


[ x_wavelets,out_wave,op] = solver_L1RLS( A, double(Ku), EPS, x0, opts );
im_recon=W'*x_wavelets; %inverse wavelet operator to get image 
% VISUALIZE

figure(3); imshow(abs([matcc(im_recon./max(im_recon(:))),matcc(linear_recon_s./max(linear_recon_s(:)))]),[]); axis off;

%
if false
    weights = findWeights( x_wavelets, .98 );
    W_weights1=opDiag(weights);
    E4=E2*W'*W_weights1
    A=linop_spot(E4,'c2c')

    
    [ x_wavelets_w,out_wave_w,op] = solver_L1RLS( A, Ksc, EPS, x0, opts );
    im_recon_w=W'*W_weights1*x_wavelets; %inverse wavelet operator to get image 

    figure(4); imshow(abs([matcc(im_recon_w./max(im_recon_w(:))),matcc(linear_recon_s./max(linear_recon_s(:)))]),[]); axis off;

end


%% ADMM LASSO
b=Ku; 
A=E3;

lambda_max = norm( A'*b, 'inf' );
lambda = 0.1*lambda_max;


[x history] = lassoADMM(A, b, lambda, 1.0, 1.0);
im_recon_admm=W'*x;

figure(3); imshow(abs(matcc(im_recon_admm)),[]); axis off;









