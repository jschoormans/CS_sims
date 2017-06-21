MR=Recon_varNSA_CS(strcat(folder,files(iter).name));
MR.Parameter.Recon.ArrayCompression='Yes'; %coil compression
MR.Parameter.Recon.ACNrVirtualChannels=4;

MR.Perform1;

MR.P.resize=true;
MR.P.fixsliceintensities=false;
MR.P.reconslices=[41:100];

MR.P.noNSAcorr=1;
MR.P.TVWeight=1e-5;
MR.P.xfmWeight=2e-5;
MR.P.Itnlim=25;
MR.P.outeriter=1;
MR.P.TGVfactor=0;
MR.ReconCS
MR.fixsliceintensities
MR.resizerecon 

imagine(abs(MR.P.Recon))