MR=Recon_varNSA_CS(strcat(folder,files(iter).name));
MR.Parameter.Recon.ArrayCompression='Yes' %coil compression
MR.Parameter.Recon.ACNrVirtualChannels=2;
MR.P.fixsliceintensities=true;
MR.P.reconslices=[300:350];

MR.Perform1;
MR.P.noNSAcorr=0;
MR.P.TVWeight=1e-5;
MR.P.xfmWeight=2e-5;
MR.P.Itnlim=25;
MR.P.outeriter=1;
MR.P.TGVfactor=0;
MR.ReconCS
MR.fixsliceintensities
MR.resizerecon 

