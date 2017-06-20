%% after loading and perform1 

% addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations'));
% cd('/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/Reconstructions/20170512_VWI_VSNA')
% file='20_12052017_1311152_8_2_wip3dicvwrfspoiledvsna2senseV4.raw'

if ispc()
cd('L:\basic\divi\Projects\cosart\scans\VNSA\VNSA60\VN_150604');
vars
else
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/scans/VNSA/VNSA60/VN_150604');
vars
end
file='vn_18042017_1923287_4_2_wip3dir3fv1V4.raw'
% file='vn_18042017_1941042_7_2_wip3dir6fv1V4.raw'
%%
MR=Recon_varNSA_CS(file);
MR.Parameter.Recon.ArrayCompression='Yes' %coil compression
MR.Parameter.Recon.ACNrVirtualChannels=4;

% MR.Parameter.Parameter2Read.dyn=[0,1].';
% MR.Parameter.Parameter2Read.ky=[-50:50].'
MR.Perform1;
%% ACTUAL CODE 
 lambdas=[1e-4 2e-4 5e-4 1e-3 2e-3 5e-3 1e-2];
for ii=1:2
    for jj=1:length(lambdas)
        
        MR.P.fixsliceintensities=true;
        MR.P.reconslices=[250];
        MR.P.noNSAcorr=ii-1;
        MR.P.TVWeight=0;%1e-5;
        MR.P.xfmWeight=lambdas(jj);%e-5;
        MR.P.Itnlim=15;
        MR.P.outeriter=3;
        MR.P.TGVfactor=0;
        MR.ReconCS
        
        img=abs(squeeze(MR.P.Recon(:,:,250)))
        noisepiece=img(1:15,1:15);
        R(ii,jj)=std(noisepiece(:))
    end
end
figure(1);
semilogx(lambdas,R.','-sq')
xlabel('lambda')
ylabel('stand deviation of noise')

figure(2);
imshow(abs(img),[0 1])

figure(3); 
imshow(abs(noisepiece),[0 1])
