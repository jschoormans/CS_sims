
clear all; close all
varsVNSA
if ispc()
    folder=('L:\basic\divi\Ima\parrec\Jasper\VNSA\ICVW_RFSPOILED_20_6\2017_06_20\20_8041\' )
else
   folder='/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/VNSA/ICVW_RFSPOILED_20_6/2017_06_20/20_8041/'  
end
files=dir([folder,'/*.raw'])

iter=5
MRC=Recon_varNSA_CS(strcat(folder,files(iter).name));
%% prepare both MR onbjects
MR=MRC.Copy;
MR.P.prewhiten=false   ;
MR.Parameter.Recon.ArrayCompression='No'; %coil compression
MR.Parameter.Recon.ACNrVirtualChannels=4;
MR.Perform1;

MR.P.WeightedL2=1;
MR.P.resize=true;
MR.P.fixsliceintensities=false;
MR.P.parallel=0;
MR.P.Itnlim=25;
MR.P.outeriter=1;
MR.P.TGVfactor=0;

MRW=MRC.Copy;
MRW.P.prewhiten=true   ;
MRW.Parameter.Recon.ArrayCompression='No'; %coil compression
MRW.Parameter.Recon.ACNrVirtualChannels=4;
MRW.Perform1;
MRW.P.WeightedL2=1;
MRW.P.resize=true;
MRW.P.fixsliceintensities=false;
MRW.P.parallel=0;
MRW.P.Itnlim=25;
MRW.P.outeriter=1;
MRW.P.TGVfactor=0;


MRW.P.TVWeight=1e-5*8;
MRW.P.xfmWeight=2e-5*8;
MRW.P.reconslices=[100:105]
MRW.ReconCS

MR.P.TVWeight=MRW.P.TVWeight;
MR.P.xfmWeight=MRW.P.xfmWeight;
MR.P.reconslices=MRW.P.reconslices
MR.ReconCS


figure(1); imshow(abs(cat(2,squeeze(MR.P.Recon(1,:,:)),squeeze(MRW.P.Recon(1,:,:)))),[])
