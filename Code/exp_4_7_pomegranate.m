

folder='/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/pomegranate/'
folder1='/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/pomegranate/scan1/'
file7_1='ph_04072016_1724533_7_2_t1wffeclearV4.raw'; 
file7='ph_04072016_1745477_7_2_wipt1wffeclearV4.raw'  %partial echo :(
file12='ph_04072016_1759344_12_2_wipt1wffeclearV4.raw' 
MR=MRecon([folder1,file7_1])
%%
MR.Parameter.Recon.ImmediateAveraging='No' 

MR.ReadData;
MR.DcOffsetCorrection;
MR.PDACorrection;
MR.RandomPhaseCorrection;
MR.MeasPhaseCorrection;
MR.SortData;
MR.RemoveOversampling;
%%
%% simulate k-space
K=MR.Data{1};
clear MR;




