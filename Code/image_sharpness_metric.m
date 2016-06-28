% Image sharpness using eigenvalues; 
%based on: ´image sharpness meausre using eigenvalues, Wee et al.'

%
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations')
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/fmeasure/'))




Image=double(imread('brain.png'))./256;
Image=squeeze(double(Image(1:589,1:589,1)));

%% EXP 1: BLURRING
clear Me
for sigmab=1:10;
I=imgaussfilt(Image,sigmab);
% [Me(sigmab)]=Sharpness_EV(I)
Me(sigmab)=fmeasure(I,'LAPD')
end
figure(1)
plot([1:10],Me)
xlabel('gaussian filter blur sigma')
figure(2)
imshow(I)
%% EXP 2: NOISE AD BLURRING
clear Me
for sigmab=1:10;
for sigman=1:3;
I=imgaussfilt(Image,sigmab);
Noise=randn(589,589)*((sigman-1)/20);
I=I+Noise;
% [Me(sigmab,sigman)]=Sharpness_EV(I)
Me(sigmab,sigman)=fmeasure(I,'GLLV')

end
end

figure(3)
plot([1:10],Me)
xlabel('gaussian filter blur sigma')
legend('noise sigma=1','noise sigma=2','noise sigma=3','noise sigma=4','noise sigma=5')
figure(4); imshow(I)