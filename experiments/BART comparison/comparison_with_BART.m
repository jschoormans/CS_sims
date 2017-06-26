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
TVWeight=1e-5
xfmWeight=0;%2*TVWeight;
MR.P.TVWeight=TVWeight
MR.P.xfmWeight=xfmWeight
MR.P.debug_nlcg=true;
MR.P.visualize_nlcg=false

MR.P.reconslices=150

MR.P.WeightedL2=0;
MR.P.VNSAlambdaCorrection=0;
MR.ReconCS
R00=MR.P.Recon;

RBart=bart('pics -l1',MR.Data(slices(2),:,:,:),conj(MR.P.sensemaps(slices(2),:,:,:)));

close all
figure(1);
imshow(abs(cat(2,squeeze(R00)./max(R00(:)),squeeze(abs(RBart./max(RBart(:)))))),[])
