MR.P.reconslices=180;
MR.P.WeightedL2=1;
MR.P.VNSAlambdaCorrection=1;

%%
MR.P.Itnlim=0

for iter=1:5
    
MR.P.Itnlim=MR.P.Itnlim+3
MR.P.VNorm=1
MR.ReconCS
RV1=MR.P.Recon;

MR.P.VNorm=0;
MR.ReconCS
RV0=MR.P.Recon;

MR.P.VNorm=2;
MR.ReconCS
RV2=MR.P.Recon;

MR.P.VNorm=1/2;
MR.ReconCS
RV12=MR.P.Recon;

I{iter}=abs(cat(2,squeeze(RV0),squeeze(RV12),squeeze(RV1),squeeze(RV2)));

end
figure(11);
imshow(cat(1,I{1},I{2},I{3},I{4},I{5}),[])
