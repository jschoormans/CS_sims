% phantom experiemnt Fluor imaging/unfolding 
% SPOT TOOLBOX 
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/spot-master'))


%

%% model problem

rr = @(I) reshape(I,[64,64])

I=phantom(64);
figure(1); subplot(311); imshow(I);axis off;

Spectrum=[1,0,0.8,0.6,0, 0, 0,0.5,0.1 0 0 0.44]; 
% Spectrum=[1,0.2]; 

figure(1); subplot(312); stem(Spectrum);axis off;

A=opConvolve(64,64,Spectrum,[0 0],'cyclic') %cyclic/ truncated?
Ic=A*I(:);
figure(1); subplot(313);  imshow(rr(Ic)); axis off;


%% inverse im space
Ir=opInverse(A)*Ic(:);
figure(2); subplot(221); imshow(rr(Ir),[]); axis off; title('pseudo-inverse')

%% inverse; kspace

F=opDFT2(64,64,1);
K=F*A*I(:);
A2=F*A; 
A2I=opInverse(A2)

Klin=opInverse(A2)*K;
figure(2); subplot(222); imshow(rr(Klin),[]); axis off; title('inverse (linear recon) -kspace')


%% recon
CGI=nl_conjgrad_fluor(A,Ic,zeros(size(Ir)),200,I(:),0,64,64); 
figure(2); subplot(223); imshow(rr(abs(CGI)),[]); axis off; title('CG-image space')

%% recon in k-space (do spectrum shift in k-space, no need for FFT all the time??) 

CGK=nl_conjgrad_fluor(A2,K,zeros(size(K)),200,I(:),0,64,64); 
figure(2); subplot(224); imshow(rr(abs(CGK)),[]); axis off; title('CG-k space')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TWO SPECTRA FOR TWO IMAGES IN ONE KSPACE
rr2 = @(I) reshape(I,[64,128])

P1 = phantom([100,0.5,0.2,0.5,0.25,60], [64]);
P2 = phantom([80,0.2,0.5,-0.5,0.5,10], [64]);
P3 = phantom([80,0.045,0.1,-0.25,0.35,10], [64]);


P4 = phantom([80,0.3,0.2,0.5,-0.25,10], [64]);
P5 = phantom([50,0.3,0.3,-0.5,-0.5,10], [64]);

Image1=P1+P2+P3;
Image2=P4+P5;

figure(4); subplot(221);hold on ; imshow([Image1,Image2],[]); title('both oils')

Spectrum2=[1 0.8 0.2 0.77 0 0 0.2 0.21 0.22 0.23 0.25 0.3 0.1 0.12]; 

B=opConvolve(64,64,Spectrum2,[0 0],'cyclic') %cyclic/ truncated?

I_corrupt=A*Image1(:)+B*Image2(:);
figure(4); subplot(223); hold on; stem(Spectrum); stem(Spectrum2); hold off; axis off;
figure(4); subplot(222); imshow(rr(I_corrupt),[])

%%  SOLVE THE TWO - OILS PROBLEM

addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Code/implementingTFOCS'))


S=[A,B]

% CG 
lambda=2
CGI2=nl_conjgrad_fluor(S,I_corrupt,(size([I_corrupt; I_corrupt])),200,[Image1(:); Image2(:)],lambda,64,128); 

figure(4); subplot(224); imshow(abs(rr2(CGI2)),[]); title('two reconstructed images (CG)')

%%
figure(5); 
imshow(rr(CGI2(1:4096,1))+rr(CGI2(4097:8192,1)),[])




















