% NEW SIMULATIONS FOR CS 
clear all; close all; clc; 
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/exportfig'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/tightfig/'))
addpath(genpath('/opt/amc/matlab/toolbox/MRecon-3.0.519'))
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/'))
addpath(genpath('/opt/amc/bart')); vars
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab'))
CC=clock; CCC=[num2str(CC(2)),'-',num2str(CC(3)),'-',num2str(CC(4)),'-',num2str(CC(5))];
figfolder=['/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Figures/',CCC]
dir=mkdir(figfolder);cd(figfolder)

%% load k-space
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/Code')
% folder='/home/jschoormans/lood_storage/divi/Projects/cosart/scans/kiwi/'
folder='/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/Reconstructions/kiwi/' 
% file='ph_25072016_2030563_13_2_wipt1wffekiwi12810dynclearV4.raw' 
% file='ph_25072016_2034581_15_2_wipt1wffekiwi128100dynclearV4.raw' %actually 256 10
file='ph_25072016_2036025_16_2_wipt1wffekiwi256100dynclearV4.raw'
% file='ph_25072016_2133129_28_2_wipt1wffekiwi51210dynclearV4.raw' %512 10 dyns
% file='ph_25072016_2125045_27_2_wipt1wffekiwi256100dynclearV4.raw'  
% file='ph_25072016_2137350_29_2_wipt1wffekiwi51230dynclearV4.raw'
% file='ph_26072016_1752392_7_2_wipt1wffe256thinsliceclearV4.raw' 
MR=MRecon([folder,file]); 
MR.Parameter.Parameter2Read.typ=1;
MR.Parameter.Parameter2Read.chan=2;
MR.Parameter.Recon.ImmediateAveraging='no'
MR.ReadData;
MR.RandomPhaseCorrection;
MR.PDACorrection;
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;
MR.SortData;
% MR.Average;
MR.RemoveOversampling;
MR.RingingFilter;
MR.ZeroFill;
% MR.K2I;
% MR.ShowData
K=MR.Data;
%%
KV=std((K),[],12);
figure(101);imshow(KV,[]); title('std of noise: does not look to be iid!')

%%
S=mean(K,12);
SNR=S./KV;
%%
figure(102); hold on; plot(abs(SNR(129,:))); plot(abs(S(129,:))); plot(abs(KV(129,:))); hold off; legend('SNR','S','N')
%%
% plot all values
figure(103);
hold on; for i=1:100
plot(real(K(129,:,i)),'.')
end; hold off
%%
KVr=std(real(K),[],12);
KVi=std(imag(K),[],12);

figure(104); hold on; plot((KVi(129,:))); plot((KVr(129,:))); plot((KV(129,:))); hold off; legend('real','imag','abs')
%%












