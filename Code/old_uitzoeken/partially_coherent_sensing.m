% partially coherent sensing
clear all; close all;

N=128;

% normal variable density masks
PDF=genPDF([N,N],4,0.25)
M=genSampling(PDF,10,10)
PSF=fftshift(fftshift(abs(fft2(M)),1),2);

figure(1)
subplot(121);
imshow(M)
subplot(122)
imshow(log(PSF),[]); colormap('summer')

%% MAKE LOCALLY COHERENT SENSING MATRIX

MP=M;
tol=0.08
for i=1:5:128
   MP(i+1,:)=MP(i,:);
   MP(i+2,:)=MP(i+1,:);
   MP(i+3,:)=MP(i+2,:);
   MP(i+4,:)=MP(i+3,:);
end
MP=MP(1:N,1:N) ; %temp fix
randmask=rand(N,N)<tol;
MP=MP+randmask; 
PSFP=fftshift(fftshift(abs(fft2(MP)),1),2);

figure(1)
subplot(121);
imshow(MP)
subplot(122)
imshow(log(PSFP),[]); colormap('summer')
%% MAKE LOCALLY COHERENT SENSING MATRIX TRY III

MP=M;
tol=0.08
for i=1:5:128
   MP(i+1,:)=MP(i,:);
   MP(i+2,:)=MP(i+1,:);
   MP(i+3,:)=MP(i+2,:);
   MP(i+4,:)=MP(i+3,:);
end
MP=MP(1:N,1:N) ; %temp fix

MP=MP(randperm(N),:);
randmask=rand(N,N)<tol;
MP=MP+randmask; 
PSFP=fftshift(fftshift(abs(fft2(MP)),1),2);

figure(1)
subplot(121);
imshow(MP)
subplot(122)
imshow(log(PSFP),[]); colormap('summer')

%% MAKE LOCALLY COHERENT SENSING MATRIX TRY II
MP=M;
tol=0.08

for i=1:3:128
    MP(i+1,:)=MP(i,:);
    MP(i+2,:)=MP(i+1,:);
%     MP(i+3,:)=MP(i+2,:);
%     MP(i+4,:)=MP(i+3,:);
end

for j=1:3:128
    MP(:,j+1)=MP(:,j);
    MP(:,j+2)=MP(:,j);
%     MP(:,j+3)=MP(:,j);
%     MP(:,j+4)=MP(:,j);
end

MP=MP(1:N,1:N) ; %temp fix
MP=MP(1:N,1:N) ; %temp fix
randmask=rand(N,N)<tol;
MP=MP+randmask; 
PSFP=fftshift(fftshift(abs(fft2(MP)),1),2);

figure(1)
subplot(121);
imshow(MP)
subplot(122)
imshow(log(PSFP),[]); colormap('summer')



%% INCOHERENCE
for i=1:N
   Incoherence(i)=corr(M(i,:).',M(1,:).') 
   Incoherence_par(i)=corr(MP(i,:).',MP(1,:).') 
end
figure(4); hold on; plot(Incoherence,'r.-'); 
plot(Incoherence_par,'k.-'); hold off

%% DO RECON
clear K
K=fftshift(fftshift(fft2(phantom(N,N)),1),2);
Image=zeros(N,N);Image(10,12)=1; Image(22,25)=1; Image(33,12)=1; Image(60,60)=0.5; Image(80,30)=1;
K=fftshift(fftshift(fft2(phantom(N,N)),1),2);
K=fftshift(fftshift(fft2(Image),1),2);


% K=ones(N,N);
K=K+100e-1.*randn(N,N);

recon=bart('pics -l1 -d5 -r0.01 -i350',double(M).*K,ones(size(M)));
recon_par=bart('pics -l1 -r0.01 -d5 -i350',double(MP).*K,ones(size(MP)));
recon=fftshift(fftshift(recon,1),2);
recon_par=fftshift(fftshift(recon_par,1),2);


figure(3); subplot(221); imshow(abs(recon),[]);
subplot(222); imshow(abs(recon_par),[]);
subplot(223); plot(abs(recon(N/2+1,:)));
subplot(224); plot(abs(recon_par(N/2+1,:)));

