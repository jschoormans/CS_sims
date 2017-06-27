clear all; close all
varsVNSA
if ispc()
    folder=('L:\basic\divi\Ima\parrec\Jasper\VNSA\VNSA_71\2017_06_16\VN_5469\' )
    cd('L:\basic\divi\Projects\cosart\CS_simulations\experiments\effect_weighting')
else
    folder=('/home/jschoormans/lood_storage/divi/ima/parrec/jasper/VNSA/VNSA_71/2017_06_16/VN_5469/' )
    cd('L:\basic\divi\Projects\cosart\CS_simulations\experiments\effect_weighting') % to change 
end
files=dir([folder,'/*.raw'])
filenumber=5;
%%
MR=Recon_varNSA_CS(strcat(folder,files(filenumber).name));
MR.Perform1;

%% figure
TVWeight=0
xfmWeight=.015;%2*TVWeight; %0.015 seems good?
MR.P.TVWeight=TVWeight
MR.P.xfmWeight=xfmWeight
MR.P.debug_nlcg=true;
MR.P.visualize_nlcg=false
MR.P.Itnlim=15
slices=[20,150]
for nslices=1:2
    
MR.P.reconslices=[slices(nslices)];
MR.P.WeightedL2=1;
MR.P.VNSAlambdaCorrection=1;
MR.ReconCS
R11{nslices}=MR.P.Recon;

MR.P.WeightedL2=1;
MR.P.VNSAlambdaCorrection=0;
MR.ReconCS
R10{nslices}=MR.P.Recon;

MR.P.WeightedL2=0;
MR.P.VNSAlambdaCorrection=0;
MR.ReconCS
R00{nslices}=MR.P.Recon;

Nfactor=7;
MR.P.TVWeight=TVWeight*Nfactor;
MR.P.xfmWeight=xfmWeight*Nfactor;
MR.P.WeightedL2=0;
MR.P.VNSAlambdaCorrection=0;

MR.ReconCS
RN{nslices}=MR.P.Recon;
end

% % PART 2: no lambda
MR.P.TVWeight=0
MR.P.xfmWeight=0
MR.P.Itnlim=7

slices=[20,180]
for nslices=1:2
    
MR.P.reconslices=[slices(nslices)];
MR.P.WeightedL2=1;
MR.P.VNSAlambdaCorrection=1;
MR.ReconCS
R11_noL{nslices}=MR.P.Recon;

MR.P.WeightedL2=0;
MR.P.VNSAlambdaCorrection=0;
MR.ReconCS
R00_noL{nslices}=MR.P.Recon;
end

% % calculate noise power spectra
ACF = @(x) conv(x,-x)
PSD = @(x) abs(fftshift(fft(ifftshift(ACF(x))))).^2
nn=30;

for ii=1:nn
PSD00(:,ii)=PSD(squeeze(R00{1}(1,:,ii)));
PSD11(:,ii)=PSD(squeeze(R11{1}(1,:,ii)));
PSD10(:,ii)=PSD(squeeze(R10{1}(1,:,ii)));
PSDN(:,ii)=PSD(squeeze(RN{1}(1,:,ii)));
end
for ii=1:nn
PSD00_noL(:,ii)=PSD(squeeze(R00_noL{1}(1,:,ii)));
PSD11_noL(:,ii)=PSD(squeeze(R11_noL{1}(1,:,ii)));
end

% figures 
resolution=0.7; % in mm; 
freqs=linspace(-1,1,599)./resolution; % spatial frequencies'

close all
figure(1);
imshow(abs(cat(2,squeeze(R00{2})./max(R00{2}(:)),squeeze(R10{2})./max(R10{2}(:)),squeeze(R11{2})./max(R11{2}(:)),squeeze(RN{2})./max(RN{2}(:)))),[])
% title('no W ,\lambda_0  |   W, \lambda_0  |  W, \lambda_c  |  no W,\lambda_N,')
export_fig '1.eps' -eps
export_fig '1.tiff' -eps

figure(3); imshow(squeeze(R00{1}(1,:,1:nn)),[]); title('image section used for calculating noise spectral density')
export_fig '3.eps' -eps

figure(4); hold on;
plot(freqs,log(sum(abs(PSD00),2)),'k');
plot(freqs,log(sum(abs(PSD10),2)),'c')
plot(freqs,log(sum(abs(PSD11),2)),'g')
plot(freqs,log(sum(abs(PSDN),2)),'y')
hold off
title('PSD'); 
xlabel('spatial frequency [mm^{-1}]')
ylabel('log(PSD)')
legend('no W ,\lambda_0','W, \lambda_0','W, \lambda_c','no W, \lambda_N')
export_fig '4.eps' -eps

%{
x10=100*((sum(abs(PSD00),2))-(sum(abs(PSD10),2)))./(sum(abs(PSD00),2));
x11=100*((sum(abs(PSD00),2))-(sum(abs(PSD11),2)))./(sum(abs(PSD00),2));
xN=100*((sum(abs(PSD00),2))-(sum(abs(PSDN),2)))./(sum(abs(PSD00),2));

figure(5); 
hold on 
plot(freqs,(x10),'g','LineWidth',1.5)
plot(freqs,(x11),'r','LineWidth',1.5)
plot(freqs,(xN),'k','LineWidth',1.5)

xlabel('spatial frequency [mm^{-1}]')
ylabel('difference (%)')
title('ratio of  Weightedl2 norm and normal l2 norm')
legend('same \lambda','adapted \lambda','extra \lambda')
export_fig '5.eps' -eps
%}

figure(11);
imshow(abs(cat(2,squeeze(R00_noL{2}),squeeze(R11_noL{2}))),[])
export_fig '11.eps' -eps
export_fig '11.tiff' -eps

figure(14); hold on;
plot(freqs,log(sum(abs(PSD00_noL),2)),'k');
plot(freqs,log(sum(abs(PSD11_noL),2)),'c')
hold off
title('PSD'); 
xlabel('spatial frequency [mm^{-1}]')
ylabel('log(PSD)')
legend('unweighted l_2 norm','weighted l_2 norm')
export_fig '14.eps' -eps


%{
x11=100*((sum(abs(PSD00_noL),2))-(sum(abs(PSD11_noL),2)))./(sum(abs(PSD00_noL),2))
figure(15); 
plot(freqs,(x11),'r','LineWidth',1.5)
xlabel('spatial frequency [mm^{-1}]')
ylabel('difference (%)')
title('ratio of  Weightedl2 norm and normal l2 norm')
legend('same \lambda','adapted \lambda','extra \lambda')
export_fig '15.eps' -eps

%}

figure(16) % mask 
imshow(abs(MR.P.MNSA),[])
export_fig '16.eps' -eps
export_fig '16.tiff' -tiff




