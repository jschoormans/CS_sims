
%%GRAPEFRUIT
folder='/home/jschoormans/lood_storage/divi/Projects/cosart/scans/fruit287/'
file='fr_28072016_1937549_9_2_wipffegrapefruit256100clearV4.raw' 

MR=MRecon([folder,file]); 
MR.Parameter.Parameter2Read.typ=1;
MR.Parameter.Parameter2Read.chan=7; % (2 for wrist coil)
MR.Parameter.Recon.ImmediateAveraging='yes'
MR.ReadData;
MR.RandomPhaseCorrection;
MR.PDACorrection;
MR.DcOffsetCorrection;
MR.MeasPhaseCorrection;
MR.SortData;
% MR.Average;
MR.RemoveOversampling;
% MR.RingingFilter;
MR.ZeroFill;
% MR.K2I;
% MR.ShowData
K=MR.Data;

K=K./max(K(:)); %normalize kspace
if true
%   K=K(97:160,97:160);
    K=K(65:190,65:190);
end


sens=ones(size(K,1),size(K,2));
reg=0.02;
Kref=mean(K,12); %still use coils though
ImRef=fftshift(squeeze(bart(['pics -RW:7:0:',num2str(reg),' -S -e -i500'],Kref,sens)),2); %to compare image
ImRef=abs(ImRef);
