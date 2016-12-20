%KLADBLOK

% recondata %256 256 8
% sensemaps %1 256 256 8 

%% from mc kspace to image
%1 mc kspace to mc image
tic
image=bart('fftmod 7',bart('fft 7',bart('fftmod 7',recondata)));
figure;imshow(squeeze(abs(image(1,:,:,4))),[])
toc
%% mc image to combined image
tic
imagec=bart('fmac -C -s8',image, sensemaps);
size(imagec)
figure(85);imshow(squeeze(abs(imagec(1,:,:))),[])
toc
%% from image to mc kspace

fksp=bart('fakeksp',imagec,recondata,sensemaps);        %without replacing
% fksp=bart('fakeksp -r',imagec,recondata,sensemaps); %with replacing

%% WAVELET REPRESENTATION OF imagec
tic
WI=bart('cdf97 6',imagec);
figure(86);imshow(squeeze(abs(WI(1,:,:))),[])
toc
IWI=bart('cdf97 -i 6',WI);
figure(87);imshow(squeeze(abs(IWI(1,:,:))),[])
toc

%% 

