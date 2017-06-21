% addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations'));
% cd('/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/Reconstructions/20170512_VWI_VSNA')
% file='20_12052017_1311152_8_2_wip3dicvwrfspoiledvsna2senseV4.raw'

clear all; close all
varsVNSA
if ispc()
    folder=('L:\basic\divi\Ima\parrec\Jasper\VNSA\VNSA_71\2017_06_16\VN_5469\' )
else
    folder=('/home/jschoormans/lood_storage/divi/ima/parrec/jasper/VNSA/VNSA_71/2017_06_16/VN_5469/' )
end
files=dir([folder,'/*.raw'])

%% RUN OVER THE WEEKEND 

for iter=4:9
MR=Recon_varNSA_CS(strcat(folder,files(iter).name));
MR.Parameter.Recon.ArrayCompression='Yes'; %coil compression
MR.Parameter.Recon.ACNrVirtualChannels=6;
MR.P.fixsliceintensities=true;

MR.Perform1;
% MR.P.reconslices=[250:280];
MR.P.noNSAcorr=0;
MR.P.TVWeight=1e-5;
MR.P.xfmWeight=2e-5;
MR.P.Itnlim=25;
MR.P.outeriter=2;
MR.P.TGVfactor=0;
MR.ReconCS;

% imagine(abs(MR.P.Recon))
MR.fixsliceintensities
% MR.resizerecon 
% imagine(abs(MR.P.Recon))

nii=make_nii(abs(MR.P.Recon))
filename=['recon',num2str(iter),'.nii']
cd(folder)
save_nii(nii,filename)

end