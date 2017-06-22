
clear all; close all
varsVNSA
if ispc()
    folder=('L:\basic\divi\Ima\parrec\Jasper\VNSA\VNSA_71\2017_06_16\VN_5469\' )
else
    folder=('/home/jschoormans/lood_storage/divi/ima/parrec/jasper/VNSA/VNSA_71/2017_06_16/VN_5469/' )
end
files=dir([folder,'/*.raw'])
filenumber=5;
%%
MR=Recon_varNSA_CS(strcat(folder,files(filenumber).name));
MR.Perform1;

%% figure
MR.P.TVWeight=2e-5
MR.P.xfmWeight=4e-5

MR.P.reconslices=[20];
MR.P.WeightedL2=1;
MR.P.VNSAlambdaCorrection=1;
MR.ReconCS
R11=MR.P.Recon;

MR.P.WeightedL2=1;
MR.P.VNSAlambdaCorrection=0;
MR.ReconCS
R10=MR.P.Recon;

MR.P.WeightedL2=0;
MR.P.VNSAlambdaCorrection=0;
MR.ReconCS
R00=MR.P.Recon;

Nfactor=4;
MR.P.TVWeight=1e-5*Nfactor
MR.P.xfmWeight=2e-5*Nfactor
MR.P.WeightedL2=0;
MR.P.VNSAlambdaCorrection=0;

MR.ReconCS
RN=MR.P.Recon;

%
close all
figure(1);
imshow(abs(cat(2,squeeze(R00),squeeze(R10),squeeze(R11),squeeze(RN))),[])
title('W/lambda, 00/ 10 /  5$\lambda$')

% calculate noise spectrum II
ACF = @(x) conv(x,-x)
PSD = @(x) fftshift(fft(ifftshift(ACF(x))))
nn=18;
figure(3); imshow(squeeze(R00(1,:,1:nn)),[])

for ii=1:nn
PSD00(:,ii)=PSD(squeeze(R00(1,:,ii)));
PSD11(:,ii)=PSD(squeeze(R11(1,:,ii)));
PSD10(:,ii)=PSD(squeeze(R10(1,:,ii)));
PSDN(:,ii)=PSD(squeeze(RN(1,:,ii)));
end
figure(4); hold on;
plot(sum(abs(PSD00),2),'k');
plot(sum(abs(PSD11),2),'c')
plot(sum(abs(PSD10),2),'g')
plot(sum(abs(PSDN),2),'y')
hold off

%
freqs=linspace(-1,1,599)
x10=100*((sum(abs(PSD00),2))-(sum(abs(PSD10),2)))%./(sum(abs(PSD00),2))
x11=100*((sum(abs(PSD00),2))-(sum(abs(PSD11),2)))%./(sum(abs(PSD00),2))
xN=100*((sum(abs(PSD00),2))-(sum(abs(PSDN),2)))%./(sum(abs(PSD00),2))

figure(5); 
hold on 
plot(freqs,(x10),'g','LineWidth',1.5)
plot(freqs,(x11),'r','LineWidth',1.5)
plot(freqs,(xN),'k','LineWidth',1.5)

xlabel('spatial frequency')
ylabel('difference (%)')
title('ratio of  Weightedl2 norm and normal l2 norm')
legend('same \lambda','adapted \lambda','extra \lambda')
%% compare quality by line plots
y=40;
figure(6); clf
plot(abs(squeeze(R00(1,:,y))),'k');
hold on; 
plot(abs(squeeze(R10(1,:,y))),'g');
% plot(abs(squeeze(R11(1,:,y))),'r');
% plot(abs(squeeze(RN(1,:,y))),'b');
