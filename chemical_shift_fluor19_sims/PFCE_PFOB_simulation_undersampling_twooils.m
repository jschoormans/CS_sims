% phantom experiemnt Fluor imaging/unfolding 
clear; close all
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/chemical_shift_fluor19_sims')
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/spot-master'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/exportfig'));
%% parameters
Nx=64;  %number of voxels
Gx=0.0004 %T/m gradient strength 
Lx=0.1; %FOV size 

rr = @(I) reshape(I,[Nx,Nx])

%% create phantom 

P1 = phantom([100,0.4,0.2,0.5,0.25,60], [64]);
P2 = phantom([80,0.2,0.4,-0.5,0.5,10], [64]);
P4 = phantom([80,0.3,0.2,0.5,-0.25,10], [64]);
P5 = phantom([50,0.3,0.3,-0.5,-0.5,10], [64]);

Image1=P1+P2;
Image2=P4+P5;
I=[Image1,Image2]

figure(1); imshow(I,[]); title('both oils')

%%
[PFCE,PFCE_alpha,PFOB,PFOB_alpha]=calcspectra(Gx,Nx,Lx)

Spectrum=zeros(1,15);
Spectrum(round(PFOB)+abs(min(round(PFOB)))+1)=PFOB_alpha./sum(PFOB_alpha(:))

offset=-(min(round(PFOB)));
A=opConvolve(Nx,Nx,Spectrum,[0 offset],'cyclic') %cyclic/ truncated?
A2=opConvolve(Nx,Nx,Spectrum.',[0 offset],'cyclic') %cyclic/ truncated?

Ic=A2*Image1(:);
figure(2);  imshow(rr(Ic),[]); axis off;

Spectrum2=[1]; 
B=opConvolve(64,64,Spectrum2,[0 0],'cyclic') %cyclic/ truncated?
B2=opConvolve(64,64,Spectrum2.',[0 0],'cyclic')

I_corrupt1=A*Image1(:)+B*Image2(:);
I_corrupt2=A2*Image1(:)+B2*Image2(:);

figure(3); imshow(rr(I_corrupt1),[])


%% inverse; kspace
F=opDFT2(64,64,1); % FOURIER OP 

%UNDERSAMPLING OPERATOR
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/CS Excercise'))
pdf=genPDF([64,64],2,0.4,2,0.2,1);
mask=genSampling(pdf,10,10);
imshow(squeeze(mask));
maskf=fftshift(mask);
l=[1:64*64];
mask2=(l(~maskf(:))); %indices of rows in FDT2 matrix that should be removed (were not sampled!)
R=opExcise(F,mask2,'rows');




K1=R*I_corrupt1;
K2=R*I_corrupt2;

K=[K1;K2];

sigma2=4;
K=K+randn(size(K))*sigma2; %add noise to K

M1=[R*A,R*B];
M2=[R*A2,R*B2];
M=[M1;M2] %measurement operator 
%%
CGK=nl_conjgrad_fluor(M,K,zeros(size(I)),50,I(:),1e-1,64,128); 

%% recon in k-space (do spectrum shift in k-space, no need for FFT all the time??) 
rr = @(I) reshape(I,[Nx,Nx*2])
rr1 = @(I) reshape(I,[Nx,Nx])

cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/chemical_shift_fluor19_sims/figs')

figure(3); imshow(rr(I),[]); axis off;
export_fig phantom.pdf -native
figure(4); imshow(rr1(R'*K(1:1630)),[]); axis off;    
export_fig Icorrupt1.pdf -native
figure(4); imshow(rr1(R'*K(64^2+1:end)),[]); axis off;
export_fig Icorrupt2.pdf -native


figure(6); imshow(rr(abs(CGK)),[]); axis off;
export_fig CG.pdf -native
