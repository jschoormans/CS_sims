MR=Recon_varNSA_CS(strcat(folder,files(iter).name));
MR.Parameter.Recon.ArrayCompression='Yes'; %coil compression
MR.Parameter.Recon.ACNrVirtualChannels=4;

MR.Perform1;
MR.P.WeightedL2=1;
MR.P.resize=false;
MR.P.fixsliceintensities=false;
MR.P.reconslices=[41:100];

MR.P.TVWeight=1e-5;
MR.P.xfmWeight=2e-5;
MR.P.Itnlim=25;
MR.P.outeriter=1;
MR.P.TGVfactor=0;
MR.ReconCS
MR.fixsliceintensities
MR.resizerecon 

imagine(abs(MR.P.Recon))

%% figure

MR.P.reconslices=[150];
MR.P.WeightedL2=1;
MR.P.VNSAlambdaCorrection=1;
MR.ReconCS
R11=MR.P.Recon;

MR.P.WeightedL2=0;
MR.P.VNSAlambdaCorrection=1;
MR.ReconCS
R01=MR.P.Recon;

MR.P.WeightedL2=1;
MR.P.VNSAlambdaCorrection=0;
MR.ReconCS
R10=MR.P.Recon;

MR.P.WeightedL2=0;
MR.P.VNSAlambdaCorrection=0;
MR.ReconCS
R00=MR.P.Recon;

figure;
imshow(abs(cat(2,squeeze(R00),squeeze(R01),squeeze(R10),squeeze(R11))),[])
title('W/lambda, 00/ 01 / 10 / 11')
