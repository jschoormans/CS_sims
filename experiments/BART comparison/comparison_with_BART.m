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

% figure
TVWeight=0;%15e-5
xfmWeight=0;%01e-4;
MR.P.TVWeight=TVWeight
MR.P.TGVfactor=0;
MR.P.xfmWeight=xfmWeight
MR.P.debug_nlcg=true;
MR.P.visualize_nlcg=false
MR.P.outeriter=1;

MR.P.reconslices=150

MR.P.WeightedL2=0;
MR.P.VNSAlambdaCorrection=0;
MR.ReconCS
R00=MR.P.Recon;

 
%% DO SCALING AS WELL 
MRC=MR.Copy;

%%
ksptemp=ifft2c(squeeze(MR.Data(150,:,:,:))); %dont really care about ffshift for this 
tmp = dimnorm(ksptemp, 3);
tmpnorm2 = sort(tmp(:), 'ascend');
% match convention used in BART
p100 = tmpnorm2(end);
p90 = tmpnorm2(round(.9 * length(tmpnorm2)));
p50 = tmpnorm2(round(.5 * length(tmpnorm2)));
if (p100 - p90) < 2 * (p90 - p50)
    scaling = p90;
else
    scaling = p100;
end
fprintf('\nScaling: %f\n\n', scaling);

ksp = MR.Data ./ scaling;
smaps=(MR.P.sensemaps);

% sensemapsorig=MR.P.sensemaps;
MRC.Data=ksp;
sensmag=reshape(vecnorm(reshape(smaps,[],8).'),[300 300 91]);

smaps=bsxfun(@rdivide,smaps,sensmag);
MRC.P.sensemaps=smaps;
MRC.ReconCS;
R00C=MRC.P.Recon;

%
RBart=bart('pics -S -d5 -i10',MR.Data(150,:,:,:),conj(MR.P.sensemaps(150,:,:,:)));
% RBart=bart('pics -S -d5 -i20 ',MRC.Data(150,:,:,:),conj(MRC.P.sensemaps(150,:,:,:)));

close all
figure(1);
imshow(abs(cat(2,squeeze(R00)./max(R00(:)),squeeze(R00C)./max(R00C(:)),squeeze(RBart./max(RBart(:))))),[])
title('CG sense comparison: non-scaled,  scaled, BART')
%% CG TV
MRC.P.outeriter=1
MRC.P.Itnlim=6;
MR.P.Itnlim=6;

lambda=0.0125; 
MR.P.TVWeight=0;
MR.P.xfmWeight=lambda.*scaling*.4
MRC.P.xfmWeight=0
MRC.P.TVWeight=lambda;
MRC.P.TGVfactor=0;

MR.ReconCS;
R00=MR.P.Recon;

MRC.P.WeightedL2=0
MRC.P.VNSAlambdaCorrection=true;

MRC.ReconCS;
R00C=MRC.P.Recon;

RBart=bart(['pics -S -d5 -i50 -RT:7:0:',num2str(lambda/5)],MR.Data(150,:,:,:),conj(MR.P.sensemaps(150,:,:,:)));
RBartscale=bart(['pics -S -d5 -i50 -RT:7:0:',num2str(lambda/5)],MRC.Data(150,:,:,:),conj(MRC.P.sensemaps(150,:,:,:)));

figure(2);
imshow(abs(cat(2,squeeze(R00)./max(R00(:)),squeeze(R00C)./max(R00C(:)),squeeze(RBart./max(RBart(:))),squeeze(RBartscale./max(RBartscale(:))))),[])
title('CG sense + TV comparison: non-scaled,  scaled, BART. BART scaled')





  